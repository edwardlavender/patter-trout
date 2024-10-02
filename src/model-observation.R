# (Temporary) ModelObsAcousticContainer R wrappers from patter-flapper

#' Identify the receiver that recorded the next detection
# This function is modified from .acs_setup_obs_receiver_id_next()
# Source:
# https://github.com/edwardlavender/patter/blob/a6c9468ae3c1cc043c41d01e5de48d06910c17e3/R/acs-internals.R
list_sensor_id_next <- function(.sensor_id) {
  rlang::check_installed("zoo")
  dt <- data.table(sensor_id = .sensor_id)
  dt$sensor_id[sapply(dt$sensor_id, is.null)] <- list(NA_integer_)
  out <-
    dt |>
    mutate(sensor_id_next = dplyr::lead(.data$sensor_id),
           sensor_id_next = zoo::na.locf(.data$sensor_id_next,
                                         fromLast = TRUE,
                                         na.rm = FALSE)) |>
    as.data.table()
  out$sensor_id_next[nrow(out)][[1]] <- NA_integer_
  out$sensor_id_next
}

#' Assemble the acoustic containers dataset
# * The default threshold is the length of of the map (metres)
assemble_acoustics_containers <- function(.acoustics,
                                          .direction,
                                          .mobility,
                                          .threshold = 10000) {
  
  # Define test arguments:
  # .acoustics = copy(tmp)
  # .direction = "backward"
  # .mobility = 750
  # .threshold = 10000
  
  # Handle inputs
  # .direction <- match.arg(.direction)
  stopifnot(.direction %in% c("forward", "backward"))
  acoustics <- copy(.acoustics)
  
  # Define the timeline
  # * If .direction == "forward", time_index runs from 1:1440
  # * If .direction == "backward", time_index runs from 1440:1
  # * Using time index to sort the time series sorts the data forwards or backwards in time
  # * We can then use the same code to define acoustic containers
  # * ... for a forward or backward filter run
  # * (On the forward run, the 'next' detection should reflect the next detection)
  # * (On the backward run, the 'next' detection should be the previous detection, to be useful)
  step       <- patter:::diffstep(sort(unique(acoustics$timestamp)))
  timestamps <- seq(min(.acoustics$timestamp),
                    max(.acoustics$timestamp),
                    by = step)
  timeline   <- data.table(timestamp = timestamps)
  if (.direction == "forward") {
    timeline <-
      timeline |>
      mutate(time_index = dense_rank(timestamp)) |>
      as.data.table()
  } else if (.direction == "backward") {
    timeline <-
      timeline |>
      mutate(time_index = dense_rank(desc(timestamp))) |>
      as.data.table()
  }
  
  # Define moorings (e.g., receiver coordinates)
  moorings <-
    acoustics |>
    group_by(.data$sensor_id) |>
    slice(1L) |>
    ungroup() |>
    as.data.table()
  
  # Define detections
  detections <-
    acoustics |>
    filter(.data$obs == 1L) |>
    group_by(.data$timestamp) |>
    # List receivers with detections at each time step
    summarise(sensor_id = list(unique(.data$sensor_id))) |>
    ungroup() |>
    arrange(.data$timestamp) |>
    # Define detection_id(s)
    mutate(detection_id = as.integer(dplyr::row_number())) |>
    select("timestamp", "sensor_id", "detection_id") |>
    as.data.table()
  
  # Define a regular timeline of the sensors that recorded detections
  # * For each time stamp, we have:
  # - A list of sensor_id that recorded detection(s) at that time stamp
  # - A list of sensor_id_next that recorded the next detection(s)
  # - max_dist_mobility, which defines the max moveable distance of the individual from those receivers
  acoustics <-
    # Add the list of sensor_id(s) and detection_id(s)
    timeline |>
    merge(detections, all.x = TRUE, by = "timestamp") |>
    # Carry the last detection_id forward, so if:
    # ... there are no detections at a given time step (detection_id = NA)
    # ... we carry forward the detection_id from the last time step
    # ... (i.e., that group of time steps belongs to the same detection)
    # ... This is required to define max_dist_mobility, below.
    arrange(.data$time_index) |>
    mutate(detection_id = as.integer(data.table::nafill(.data$detection_id, type = "locf"))) |>
    # The individual may be within receiver_gamma + (.mobility * time steps) from the next sensor
    group_by(.data$detection_id) |>
    arrange(.data$time_index, .by_group = TRUE) |>
    mutate(max_dist_mobility = .mobility * rev(dplyr::row_number())) |>
    ungroup() |>
    # At each time step, list the receivers that recorded the 'next' detection
    # * For forward filter runs, the 'next' detection is the next detection
    # * For backward filter runs, the 'next' detection is the previous detection
    arrange(.data$time_index) |>
    mutate(sensor_id_next = list_sensor_id_next(.data$sensor_id)) |>
    as.data.table()
  
  # View(acoustics[, .(timestamp, detection_id, sensor_id, sensor_id_next, max_dist_mobility)])
  
  # Define acoustic containers dataset
  # * The contains the following columns:
  # - timestamp (required)
  # - sensor_id
  # - obs (nominally 1, unused)
  # - receiver_x, receiver_y, max_dist (the additionally required ModelObs structure parameters)
  containers <-
    acoustics |>
    mutate(obs = 1L) |>
    select("timestamp", "obs", sensor_id = "sensor_id_next", "max_dist_mobility") |>
    tidyr::unnest(cols = "sensor_id") |>
    arrange(.data$timestamp, .data$sensor_id) |>
    mutate(receiver_x = moorings$receiver_x[match(.data$sensor_id, moorings$sensor_id)],
           receiver_y = moorings$receiver_y[match(.data$sensor_id, moorings$sensor_id)],
           # We assume that receiver_gamma is constant in time for each receiver for convenience
           receiver_gamma = moorings$receiver_gamma[match(.data$sensor_id, moorings$sensor_id)],
           max_dist = .data$receiver_gamma + .data$max_dist_mobility) |>
    select("timestamp", "obs", "sensor_id", "receiver_x", "receiver_y", "max_dist") |>
    filter(!is.na(sensor_id) & !is.na(max_dist)) |>
    # For speed, we only implement acoustic containers
    # ... when the distance an individual must be from a receiver is < .threshold
    # ... (e.g., threshold may be the size of the study area)
    filter(max_dist < .threshold) |>
    as.data.table()
  
  # Checks
  if (FALSE) {
    if (.direction == "forward") {
      
      # On the forward run, the first detection is:
      detections[1, ]
      # As we approach the first detection, max_dist should shrink
      # ... & at the time stamp immediately preceding the detection,
      # ... max_dist should be receiver_gamma + mobility
      containers[timestamp < detections$timestamp[1], ]
      
    } else if (.direction == "backward") {
      
      # On the backward run, the first detection is effectively:
      detections[.N, ]
      
      # As we approach this detection, moving backwards in time, max_detection should shrink:
      # ... & at the time stamp immediately preceding (after!) the detection,
      # ... max_dist should be receiver_gamma + mobility
      containers[timestamp > detections$timestamp[nrow(detections)], ]
      
    }
  }
  
  
  # Visualise the detections & containers data.table
  # Pick an example detection
  # In containers, work towards (forward or backward) to that time & check values
  # View(detections); View(containers)
  
  # Return the dataset
  containers
  
}
