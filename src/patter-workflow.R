#' Patter workflow

# Run the particle algorithms for a selected individual
# This workflow assumes:
# * Select packages and functions are available & loaded
# * Julia is connected
# * the map & vmap exist in the Julia session

patter_workflow <- function(id, 
                            map = NULL, map_len,
                            timelines, paths, moorings, detections, metadata, 
                            parameters, model_move_type = c("low", "medium", "high"),
                            interactive = TRUE) {
  
  # id          <- 2L
  # interactive <- TRUE
  model_move_type    <- match.arg(model_move_type)

  
  ###########################
  #### Prepare data
  
  #### Define individual datasets
  sim_id <- NULL
  path   <- paths[sim_id == id, ]
  dets   <- detections[sim_id == id, ]
  meta   <- metadata[sim_id == id, ]
  stopifnot(nrow(meta) == 1L)
  
  #### Define timeline (individual-specific)
  timeline <- timelines[[id]]
  
  #### (optional) Examine individual-specific data
  if (!os_linux() && interactive) {
    
    # Examine simulated step lengths 
    dist <- terra::distance(cbind(path$x, path$y), 
                            lonlat = FALSE, sequential = TRUE)
    max(dist)
    hist(dist)
    
    # Examine turning angles
    # * (optional) TO DO
  }
  
  
  ###########################
  #### Define movement model 

  #### Random walk
  # This is a simple generalisable model
  # But it is associated with convergence issues 
  # state      <- "StateXY"
  # sspec    <- "gamma"
  mobility <- parameters$model_move$mobility
  # sshape   <- parameters$model_move$step$shape
  # sscale   <- parameters$model_move$step$scale
  # amean    <- parameters$model_move$angle$mean
  # asd      <- parameters$model_move$angle$sd
  # model_move <- move_xy(mobility = mobility, 
  #                       dbn_length = glue("truncated(Gamma({sshape}, {sscale}), lower = 0.0, upper = {mobility})"), 
  #                       dbn_angle = glue("truncated(Normal(0, 1), lower = -pi, upper = pi)"))

  #### Correlated random walk (generalised)
  # A CRW has better convergence properties
  # However, it is difficult to capture the step length distribution adequately across all individuals
  # The turning angle distribution is acceptible:
  # * This is inflated slightly to limit performance cost
  # * Overly high concentration leads to slow algorithm runs
  # ... as many trials are required to avoid jumps onto land
  # * But a sd of <= 0.3 appears to be required for convergence
  # state      <- "StateXYD"
  # sspec      <- "gamma"
  # sshape     <- sshape
  # sscale     <- 60
  # amean      <- 0.0
  # asd        <- 0.3
  # state      <- "StateXYD"
  # model_move <- move_xyd(mobility = mobility, 
  #                        dbn_length = glue("truncated(Gamma({sshape}, {sscale}), lower = 0.0, upper = {mobility})"), 
  #                        dbn_angle_delta = glue("Normal({amean}, {asd})"))

  #### Correlated random walk (low activity)
  # TO DO
  if (model_move_type == "low") {
    state      <- "StateXYD"
    sspec      <- "gamma"
    sshape     <- 3
    sscale     <- 6
    amean      <- 0.0
    asd        <- 0.3
    model_move <- move_xyd(mobility = mobility, 
                           dbn_length = glue("truncated(Normal({sshape}, {sscale}), lower = 0.0, upper = {mobility})"), 
                           dbn_angle_delta = glue("Normal({amean}, {asd})"))
    
    
    curve(dtrunc(x, spec = "gamma", a = 0, b = mobility, shape = sshape, scale = sscale), 
          from = 0, to = mobility, n = 1e3L,
          xlab = "Step length (m)", ylab = "Density")
    abline(v = 12.7)
  }
  
  #### Correlated random walk (medium activity)
  if (model_move_type == "medium") {
    state      <- "StateXYD"
    sspec       <- "norm"
    sshape     <- 63.5
    sscale     <- 5
    amean      <- 0.0
    asd        <- 0.3
    model_move <- move_xyd(mobility = mobility, 
                           dbn_length = glue("truncated(Normal({sshape}, {sscale}), lower = 0.0, upper = {mobility})"), 
                           dbn_angle_delta = glue("Normal({amean}, {asd})"))
  }

  #### Correlated random walk (high activity)
  if (model_move_type == "high") {
    state      <- "StateXYD"
    sspec      <- "gamma"
    sshape     <- 20
    sscale     <- 10
    amean      <- 0.0
    asd        <- 0.3
    model_move <- move_xyd(mobility = mobility, 
                           dbn_length = glue("truncated(Gamma({sshape}, {sscale}), lower = 0.0, upper = {mobility})"), 
                           dbn_angle_delta = glue("Normal({amean}, {asd})"))
  }

  #### Visualise movement model 
  if (interactive) {
    
    # Visualise model components 
    pp <- par(mfrow = c(1, 2))
    if (sspec == "gamma") {
      curve(dtrunc(x, spec = "gamma", a = 0, b = mobility, shape = sshape, scale = sscale), 
            from = 0, to = mobility, n = 1e3L,
            xlab = "Step length (m)", ylab = "Density")
    } else if (sspec == "norm") {
      curve(dtrunc(x, spec = "norm", a = 0, b = mobility, mean = sshape, sd = sscale), 
            from = 0, to = mobility, n = 1e3L,
            xlab = "Step length (m)", ylab = "Density")
    }
    abline(v = max(dist), col = "red")
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
  
  
  ###########################
  #### Define observation models
  
  #### Assemble acoustics (0, 1)
  moorings[, receiver_gamma := parameters$model_obs$receiver_gamma]
  acoustics   <- assemble_acoustics(.timeline   = timeline, 
                                    .detections = dets, 
                                    .moorings   = moorings)
  
  #### Assemble containers
  containers  <- assemble_acoustics_containers(.timeline  = timeline, 
                                               .acoustics = acoustics, 
                                               .mobility  = mobility, 
                                               .threshold = map_len)
  
  #### Define yobs (forward & backward)
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
  # Tune .n_move, .n_particle & .n_record for improved speed during filtering/smoothing
  args <- list(.timeline = timeline, 
               .state = state, 
               .xinit = xinit_fwd, 
               .yobs = yobs_fwd,
               .model_move = model_move, 
               .n_move = 10000L,
               .n_particle = 2.5e4L,
               .n_record = 500L,
               .n_iter = 1L,
               .direction = "forward")
  
  #### Run forward filter 
  # Setting initial observations is slow (~1 min)
  # Set yobs to missing in `args` to speed up multiple runs 
  t1             <- Sys.time()
  set_seed()
  fwd            <- do.call(pf_filter, args, quote = TRUE)
  t2             <- Sys.time()
  convergence_dt <- data.table(sim_id      = id, 
                               direction   = "forward", 
                               model_move  = model_move_type,
                               np          = args$.n_particle, 
                               nt          = length(timeline),
                               time        = as.numeric(difftime(t2, t1, units = "mins")),
                               convergence = fwd$convergence, 
                               trials      = fwd$trials)
  qs::qsave(convergence_dt, 
            here_output_sim("convergence", model_move_type, glue("convergence-fwd-{id}.qs")))
  if (!fwd$convergence) {
    return(NULL)
  }
  # qs::qsave(fwd, here_output_sim("tmp", "tmp.qs"))
  # file.size(here_output_sim("tmp", "tmp.qs")) * 1e-9 # 4 GB
  
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
                                model_move  = model_move_type,
                                np          = args$.n_particle, 
                                nt          = length(timeline),
                                time        = as.numeric(difftime(t2, t1, units = "mins")),
                                convergence = bwd$convergence, 
                                trials      = bwd$trials)
  qs::qsave(convergence_dt, 
            here_output_sim("convergence", model_move_type, glue("convergence-bwd-{id}.qs")))
  if (!bwd$convergence) {
    return(NULL)
  }
  
  #### Run smoother 
  # (TO DO) ETA without vmap: XXX mins
  # (TO DO) ETA with vmap:    XXX mins
  t1   <- Sys.time()
  set_seed()
  smo  <- pf_smoother_two_filter()
  t2   <- Sys.time()
  convergence_dt  <- data.table(sim_id      = id, 
                                model_move  = model_move_type,
                                direction   = "smoothing", 
                                np          = args$.n_particle, 
                                nt          = length(timeline),
                                time        = as.numeric(difftime(t2, t1, units = "mins")),
                                convergence = NA_integer_, 
                                trials      = NA_integer_)
  qs::qsave(convergence_dt, 
            here_output_sim("convergence", model_move_type, glue("convergence-smo-{id}.qs")))
  qs::qsave(smo, here_output_sim("particles", model_move_type, glue("smo-{id}.qs")))
  # file_size(here_output_sim("particles", model_move_type, glue("smo-{id}.qs"))) # 757 MB
  
  NULL
  
}

