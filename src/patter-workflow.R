# Run the particle algorithms for a selected individual
# This workflow assumes:
# * Select packages and functions are available & loaded
# * Julia is connected
# * the map & vmap exist in the Julia session

patter_workflow <- function(id, 
                            map = NULL, map_len,
                            timelines, paths, moorings, detections, metadata, 
                            parameters,
                            interactive = TRUE) {
  
  # id   <- 1
  # interactive <- TRUE
  
  ###########################
  #### Prepare algorithms
  
  #### Define individual datasets
  sim_id <- NULL
  path   <- paths[sim_id == id, ]
  dets   <- detections[sim_id == id, ]
  meta   <- metadata[sim_id == id, ]
  stopifnot(nrow(meta) == 1L)
  
  #### Define timeline (individual-specific)
  timeline <- timelines[[id]]
  
  #### (optional) Examine individual-specific data
  if (interactive) {
    
    # Examine simulated step lengths 
    dist <- terra::distance(cbind(path$x, path$y), 
                            lonlat = FALSE, sequential = TRUE)
    max(dist)
    hist(dist)
    
    # Examine turning angles
    # * (optional) TO DO
  }
  
  #### Define movement model 
  
  # Parameters
  mobility <- parameters$model_move$mobility
  sshape   <- parameters$model_move$step$shape
  sscale   <- parameters$model_move$step$scale
  amean    <- parameters$model_move$angle$mean
  asd      <- parameters$model_move$angle$sd
  
  # Random walk
  # > This is associated with convergence issues
  # state      <- "StateXY"
  # model_move <- move_xy(mobility = mobility, 
  #                       dbn_length = glue("truncated(Gamma({sshape}, {sscale}), lower = 0.0, upper = {mobility})"), 
  #                       dbn_angle = glue("truncated(Normal(0, 1), lower = -pi, upper = pi)"))
  
  # Correlated random walk
  # > This has better convergence properties
  # > The turning angle concentration is turned to limit performance cost
  # > Overly high concentration leads to slow algorithm runs
  # > ... as many trials are required to avoid jumps onto land
  state      <- "StateXYD"
  model_move <- move_xyd(mobility = mobility, 
                         dbn_length = glue("truncated(Gamma({sshape}, {sscale}), lower = 0.0, upper = {mobility})"), 
                         dbn_angle_delta = glue("Normal({0.0}, {0.4})"))

  #### Visualise movement model 
  if (interactive) {
    
    # Visualise model components 
    pp <- par(mfrow = c(1, 2))
    curve(dtrunc(x, spec = "gamma", a = 0, b = mobility, shape = sshape, scale = sscale), 
          from = 0, to = mobility, n = 1e3L,
          xlab = "Step length (m)", ylab = "Density")
    abline(v= max(dist), col = "red")
    curve(dnorm(x, 0, 0.25),
          from = -pi - 0.1, to = pi + 0.1, n = 1e3L,
          xlab = "Turning angle (rad)", ylab = "Density")
    par(pp)
    
    # Visualise realisations of the movement model
    length(timeline)
    pos <- 1:50000
    if (!os_linux() & !is.null(map)) {
      sim_path_walk(.map = map, 
                    .timeline = timeline[pos],
                    .state = state, 
                    .model_move = model_move, 
                    .n_path = 6L, .one_page = TRUE)
    }
    
    # Compare to simulated paths
    if (!os_linux() & !is.null(map)) {
      pp <- par(mfrow = c(2, 2))
      lapply(1:4, function(id) {
        path <- paths[sim_id == id, ]
        path <- path[pos, ]
        terra::plot(map)
        patter:::add_sp_path(path$x, path$y)
      })
      par(pp)
    }
    
  }
  
  #### Define observation models
  # Assemble acoustics (0, 1)
  moorings[, receiver_gamma := parameters$model_obs$receiver_gamma]
  acoustics   <- assemble_acoustics(.timeline   = timeline, 
                                    .detections = dets, 
                                    .moorings   = moorings)
  # Assemble containers
  containers  <- assemble_acoustics_containers(.timeline  = timeline, 
                                               .acoustics = acoustics, 
                                               .mobility  = mobility, 
                                               .threshold = map_len)
  # Define yobs (forward & backward)
  yobs_fwd <- list(ModelObsAcousticLogisTrunc = copy(acoustics),
                   ModelObsAcousticContainer  = copy(containers$forward))
  yobs_bwd <- list(ModelObsAcousticLogisTrunc = copy(acoustics),
                   ModelObsAcousticContainer  = copy(containers$backward))
  
  #### Define initial states
  xinit_fwd <- xinit_bwd <- NULL
  
  
  ###########################
  ###########################
  #### Algorithm runs 
  
  #### Define filter arguments 
  args <- list(.timeline = timeline, 
               .state = state, 
               .xinit = xinit_fwd, 
               .yobs = yobs_fwd,
               .model_move = model_move, 
               .n_particle = 2.5e4L,
               .n_record = 500L,
               .n_iter = 5L,
               .direction = "forward"
  )
  
  #### Run forward filter 
  # Setting initial observations is slow (~1 min)
  # Set yobs to missing in `args` to speed up multiple runs 
  t1             <- Sys.time()
  set_seed()
  fwd            <- do.call(pf_filter, args, quote = TRUE)
  t2             <- Sys.time()
  convergence_dt <- data.table(sim_id      = id, 
                               direction   = "forward", 
                               time        = as.numeric(difftime(t2, t1, units = "mins")),
                               convergence = fwd$convergence, 
                               trials      = fwd$trials)
  qs::qsave(convergence_dt, 
            here_output("convergence", glue("convergence-fwd-{id}.qs")))
  if (!fwd$convergence) {
    return(NULL)
  }
  # qs::qsave(fwd, here_output("tmp", "tmp.qs"))
  # file.size(here_output("tmp", "tmp.qs")) * 1e-9 # 4 GB
  
  #### Run backward filter
  args$.xinit     <- xinit_bwd
  args$.yobs      <- yobs_bwd
  args$.direction <- "backward"
  t1              <- Sys.time()
  set_seed()
  bwd             <- do.call(pf_filter, args, quote = TRUE)
  t2              <- Sys.time()
  convergence_dt  <- data.table(sim_id      = id, 
                                direction   = "backward", 
                                time        = as.numeric(difftime(t2, t1, units = "mins")),
                                convergence = bwd$convergence, 
                                trials      = bwd$trials)
  qs::qsave(convergence_dt, 
            here_output("convergence", glue("convergence-bwd-{id}.qs")))
  if (!bwd$convergence) {
    return(NULL)
  }
  
  #### Run smoother 
  # ETA without vmap: XXX mins
  # ETA with vmap:    XXX mins
  t1   <- Sys.time()
  set_seed()
  smo  <- pf_smoother_two_filter()
  t2   <- Sys.time()
  convergence_dt  <- data.table(sim_id      = id, 
                                direction   = "smoothing", 
                                time        = as.numeric(difftime(t2, t1, units = "mins")),
                                convergence = NA_integer_, 
                                trials      = NA_integer_)
  qs::qsave(convergence_dt, here_output(glue("convergence-smo-{id}.qs")))
  qs::qsave(smo, here_output("particles", glue("smo-{id}.qs")))
  # file_size(here_output("particles", glue("smo-{id}.qs")))
  
  NULL
  
}

