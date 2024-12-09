###########################
###########################
#### patter-setup-sim.R

#### Aims
# 1) Sets up data for analysis with patter

#### Prerequisites
# 1) Create recs_short -> simulated_moorings.rds in model_occupancy_estimates.R


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())
# try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
Sys.setenv("JULIA_SESSION" = FALSE)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(patter)
library(tictoc)
dv::src()

#### Load data 
champlain  <- terra::vect(here_data("ChamplainRegionsGrouped/ChamplainRegionsGrouped.shp"))
map        <- terra::rast(here_input("map.tif"))
map_bbox   <- qs::qread(here_input("map-bbox.qs"))
moorings   <- readRDS(here_data("simulated_moorings.rds"))
detections <- readRDS(here_data("simulated_detections.rds"))
metadata   <- readRDS(here_data("simulations_metadata.rds"))


###########################
###########################
#### Simulated paths (~40 s)

start <- as.POSIXct("2022-01-01 00:00:00" , tz = "UTC")

file_paths <- here_input_sim("paths.qs")
if (!file.exists(file_paths)) {
  
  tic()
  
  #### Read paths 
  # * `complete_simulated_transmissions_regions.rds` contains the full trajectory
  # > I.e., the position every 127 s
  # * `simulated_positions_July2023.rds` contains the simulated positions only
  # > These positions can be more sparsely spaced depending on velocity
  # > (glatos imposes a timeline afterwards)
  # paths <- readRDS("data/simulated_positions_July2023.rds")
  paths   <- readRDS("data/complete_simulated_transmissions_regions.rds")
  
  #### Define path coordinates (UTM)
  sxy <- 
    paths$geometry |>
    sf::st_coordinates() |>
    terra::vect(crs = terra::crs(champlain)) |> 
    terra::project(terra::crs(map)) |> 
    terra::crds()
  
  #### Define paths data.table
  paths <- 
    paths |> 
    as.data.table() |>
    mutate(timestamp = start + time, 
           x = sxy[, 1], y = sxy[, 2]) |>
    group_by(virt_fish) |> 
    mutate(timestep = row_number()) |> 
    select(individual_id = virt_fish, timestep, timestamp, x, y, region = regions) |> 
    as.data.table()
  
  #### Save 
  qs::qsave(paths, file_paths)
  toc()
  
} else {
  paths <- qs::qread(file_paths)
}


###########################
###########################
#### Simulated observations

#### Define detection data (detections & receiver coordinates)
# Define receiver coordinates (UTM)
rxy <- 
  cbind(moorings$deploy_lon, moorings$deploy_lat) |> 
  terra::vect(crs = terra::crs(champlain)) |> 
  terra::project("EPSG:3175") |>
  terra::crds()
# points(rxy[, 1], rxy[, 2], col = "red")
# stopifnot(all(!is.na(terra::extract(map, rxy)$map_value)))

#### Define moorings data.table
range(detections$detection_timestamp_utc)
moorings <- 
  moorings |>
  rename(receiver_lon = deploy_lon, 
         receiver_lat = deploy_lat) |>
  mutate(receiver_id = row_number(), 
         receiver_start = as.POSIXct("2022-01-01 00:00:00", tz = "UTC"),
         receiver_end = as.POSIXct("2022-10-18 00:00:00 UTC", tz = "UTC"),
         receiver_x = rxy[, 1], 
         receiver_y = rxy[, 2], 
         receiver_alpha = 1.8, 
         receiver_beta = -1/200,
         # Set receiver_gamma to some high value
         # * (Untruncated model used for simulation)
         # * This value is based on explore-sim.R
         receiver_gamma = 3240, 
         receiver_key = paste(receiver_lon, receiver_lat)) |> 
  select(receiver_id, 
         receiver_start, receiver_end,
         receiver_x, receiver_y,
         receiver_alpha, receiver_beta, receiver_gamma,
         receiver_lon, receiver_lat, receiver_key) |> 
  as.data.table()

#### Update detections 
detections <- 
  detections |> 
  select(individual_id = sim_id, 
         timestamp = "detection_timestamp_utc",
         receiver_lon = deploy_lon, 
         receiver_lat = deploy_lat) |> 
  # Define individual_id as integer
  mutate(individual_id = stringr::str_replace(individual_id, "sim_", ""), 
         individual_id = as.integer(individual_id)) |>
  mutate(receiver_key = paste(receiver_lon, receiver_lat)) |> 
  mutate(receiver_id = moorings$receiver_id[match(receiver_key, moorings$receiver_key)]) |>
  select(individual_id, timestamp, receiver_id) |> 
  as.data.table()

#### Clean up moorings
moorings <- 
  moorings |> 
  select("receiver_id", 
         "receiver_start", "receiver_end", 
         "receiver_x", "receiver_y", 
         "receiver_alpha", "receiver_beta", "receiver_gamma") |> 
  as.data.table()

#### Checks
stopifnot(all(!is.na(detections$receiver_id)))

#### Compute container thresholds
# See explanation in setup-real.R
dist <- terra::distance(cbind(moorings$receiver_x, moorings$receiver_y), map_bbox, lonlat = FALSE)
dist <- apply(dist, 1, max)
cthresholds <- data.table(receiver_id = moorings$receiver_id, 
                          distance = dist)


###########################
###########################
#### Metadata

#### Define metadata
metadata <-
  metadata |> 
  # Define individual_id as integer
  mutate(individual_id = stringr::str_replace(sim_id, "sim_", ""), 
         individual_id = as.integer(individual_id), 
         sim_id = NULL) |>
  as.data.table()

#### Collect individual IDs (1:100)
ids <- sort(unique(metadata$individual_id))


###########################
###########################
#### Timeline 

#### Simulation methods
# The timeline is individual specific
# Start time: "2022-01-01 00:00:00"
# Simulated tracks comprised 5000 steps
# Each step comprised 500 m
# But velocity was set to different values in glatos::transmit_along_path()
# I.e., For each individual, the duration of a step length is different 
# Transmissions were generated along the paths every 120 + 7 s
# To assemble the timeline, we can use the simulated paths dataset

#### Define study start time
start <- as.POSIXct("2022-01-01 00:00:00" , tz = "UTC")

#### Define timelines (a list, with one element for each individual)
timelines <- 
  lapply(ids, function(id) {
    # print(id)
    path <- paths[individual_id == id, ]
    timeline <- assemble_timeline(list(data.table(timestamp = start), path),
                                  .step = "2 mins", 
                                  .trim = FALSE)
    timeline
}) 
names(timelines) <- ids
# Examine timeline ranges
timestats <-
  data.table(individual_id = 1:100L, 
             duration = sapply(timelines, 
                               \(x)  as.numeric(difftime(max(x), min(x), units = "weeks"))))
utils.add::basic_stats(timestats$duration)


###########################
###########################
#### Save datasets

#### Write datasets
qs::qsave(start, here_input_sim("start.qs"))
qs::qsave(timelines, here_input_sim("timelines.qs"))
qs::qsave(moorings, here_input_sim("moorings.qs"))
qs::qsave(detections, here_input_sim("detections.qs"))
qs::qsave(cthresholds, here_input_sim("cthresholds.qs"))
qs::qsave(metadata, here_input_sim("metadata.qs"))


#### End of code. 
###########################
###########################