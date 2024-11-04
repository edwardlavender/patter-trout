#' Patter workflow

# Run the particle algorithms for a selected individual
# This workflow assumes:
# * Select packages and functions are available & loaded
# * Julia is connected
# * the map & vmap exist in the Julia session

patter_workflow <- function(id, 
                            timelines = NULL, 
                            moorings, detections, cthresholds,
                            model_move_type = c("sim-low", "sim-medium", "sim-high", "real"),
                            n_particle = 1000L,
                            trial = FALSE,
                            here_output) {
  
  #### Print progress
  cat("\n\n\n\n")
  print(paste0(rep("-", 50), collapse = ""))
  print(Sys.time())
  print(paste(which(unique(detections$individual_id) == id), "/", length(unique(detections$individual_id))))
  print(id)
  cat("\n")
  
  #### Function set up
  # id            <- 1L
  model_move_type <- match.arg(model_move_type)
  
  #### Build directories
  folders <- c(here_output("convergence", model_move_type),
               here_output("particles", model_move_type))
  sapply(folders, function(folder) {
    if (!dir.exists(folder)) {
      warning(glue("Directory '{folder}' created."))
      dir.create(folder, recursive = TRUE)
    }
    NULL
  })

  
  ###########################
  #### Prepare data
  
  #### Define individual datasets
  dets   <- detections[individual_id == id, ]
  
  #### Define timeline (individual-specific)
  # Use pre-prepared timeline for simulations
  # Otherwise, assemble timeline 
  if (!is.null(timelines)) {
    timeline <- timelines[[id]]
  } else {
    timeline <- assemble_timeline(.datasets = list(dets), .step = "2 mins")
  }
  
  
  ###########################
  #### Define movement model 

  #### (1) Simulation models
  # Step lengths:
  # * It is difficult to capture step lengths adequately across all individuals
  # * Three models (low, medium, high activity) are required
  # Turning angle:
  # * A RW is associated with convergence challenges
  # * A CRW has better convergence properties
  # * A generic model for turning angle is acceptable: 
  # - This is inflated slightly to limit performance cost
  # - Overly high concentration leads to slow algorithm runs
  # ... as many trials are required to avoid jumps onto land
  # - But an sd of <= 0.3 appears to be required for convergence

  #### (A) Simulation low-activity model
  mobility <- 115
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
  }
  
  #### (B) Simulation medium-activity model 
  if (model_move_type == "medium") {
    state      <- "StateXYD"
    sspec      <- "norm"
    sshape     <- 63.5
    sscale     <- 5
    amean      <- 0.0
    asd        <- 0.3
    model_move <- move_xyd(mobility = mobility, 
                           dbn_length = glue("truncated(Normal({sshape}, {sscale}), lower = 0.0, upper = {mobility})"), 
                           dbn_angle_delta = glue("Normal({amean}, {asd})"))
  }

  #### (C) Simulation high-activity model 
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
  
  #### (2) Real-world model
  if (model_move_type == "real") {
    state      <- "StateXYD"
    # Trial max mobility: 1 m/s for 3 min
    mobility   <- 180 
    model_move <- move_xyd(mobility = mobility, 
                           dbn_length = glue("truncated(Gamma({3.0}, {1/0.1}), lower = 0.0, upper = {mobility})"), 
                           dbn_angle_delta = glue("Normal({0.0}, {0.3})"))
  }
  
  #### Visualise movement model
  # See supporting code in dev/
  
  
  ###########################
  #### Define observation models
  
  #### Assemble acoustics (0, 1)
  acoustics   <- assemble_acoustics(.timeline   = timeline, 
                                    .detections = dets, 
                                    .moorings   = moorings)
  
  #### Assemble containers
  containers  <- assemble_acoustics_containers(.timeline  = timeline, 
                                               .acoustics = acoustics, 
                                               .mobility  = mobility, 
                                               .threshold = 175000)
  # Filter containers by pre-computed cthresholds
  # * Precomputed thresholds are used as `.map` is not currently supported on linux
  # * (when JULIA_SESSION = "TRUE")
  for (direction in c("forward", "backward")) {
    containers[[direction]] <- 
      containers[[direction]] |>
      lazy_dt(immutable = FALSE) |>
      left_join(cthresholds, by = c("sensor_id" = "receiver_id")) |> 
      filter(radius <= distance) |> 
      mutate(distance = NULL) |> 
      as.data.table()
  }

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
  n_record <- ifelse(trial, 1L, 1000L)
  args     <- list(.timeline    = timeline, 
                   .state       = state, 
                   .xinit       = xinit_fwd, 
                   .yobs        = yobs_fwd,
                   .model_move  = model_move, 
                   .n_move      = 10000L,
                   .n_particle  = n_particle,
                   .n_record    = n_record,
                   .n_iter      = 1L,
                   .direction   = "forward")
  
  #### Run forward filter 
  t1             <- Sys.time()
  set_seed()
  fwd            <- do.call(pf_filter, args, quote = TRUE)
  t2             <- Sys.time()
  if (trial) {
    return(fwd$convergence)
  }
  convergence_dt <- data.table(individual_id = id, 
                               direction    = "forward", 
                               model_move   = model_move_type,
                               np           = args$.n_particle, 
                               nt           = length(timeline),
                               time         = as.numeric(difftime(t2, t1, units = "mins")),
                               convergence  = fwd$convergence, 
                               trials       = fwd$trials)
  qs::qsave(convergence_dt, 
            here_output("convergence", model_move_type, glue("convergence-fwd-{id}.qs")))
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
  convergence_dt  <- data.table(individual_id      = id, 
                                direction   = "backward", 
                                model_move  = model_move_type,
                                np          = args$.n_particle, 
                                nt          = length(timeline),
                                time        = as.numeric(difftime(t2, t1, units = "mins")),
                                convergence = bwd$convergence, 
                                trials      = bwd$trials)
  qs::qsave(convergence_dt, 
            here_output("convergence", model_move_type, glue("convergence-bwd-{id}.qs")))
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
  convergence_dt  <- data.table(individual_id = id, 
                                model_move    = model_move_type,
                                direction     = "smoothing", 
                                np            = args$.n_particle, 
                                nt            = length(timeline),
                                time          = as.numeric(difftime(t2, t1, units = "mins")),
                                convergence   = NA_integer_, 
                                trials        = NA_integer_)
  qs::qsave(convergence_dt, 
            here_output("convergence", model_move_type, glue("convergence-smo-{id}.qs")))
  qs::qsave(smo, here_output("particles", model_move_type, glue("smo-{id}.qs")))
  # file_size(here_output("particles", model_move_type, glue("smo-{id}.qs"))) # 757 MB
  
  NULL
  
}

