###########################
###########################
#### patter-setup.R

#### Aims
# 1) Sets up data for analysis with patter

#### Prerequisites
# 1) Create recs_short -> moorings.rds in model_occupancy_estimates.R


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())
# try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tictoc)
dv::src()

#### Load data 
champlain  <- terra::vect("data/ChamplainRegions.shp")
moorings   <- readRDS("data/moorings.rds")
detections <- readRDS("data/simulated_detections.rds")
metadata   <- readRDS("data/simulations_metadata.rds")


###########################
###########################
#### Study area (~2 s)

#### Build map 
# Use a coarse map for speed sampling initial locations 
tic()
epsg_utm         <- "EPSG:3175"
champlain$region <- champlain$GNIS_NAME
champlain$land   <- as.numeric(1)
regions <- as_SpatRaster(champlain, .simplify = 0.001, .utm = epsg_utm,
                         .field = "region", .res = 250, .plot = TRUE)
maps    <- as_SpatRaster(champlain, .simplify = 0.001, .utm = epsg_utm,
                         .field = "land", .res = 250, .plot = TRUE)
map     <- maps$SpatRaster
map_len <- terra::ymax(map) - terra::ymin(map)
toc()

#### Examine map properties
# Visualise map simplification
terra::plot(map, col = "blue")
terra::lines(champlain |> terra::project(epsg_utm))
# Check ncell & compare to dat_gebco() for reference
terra::ncell(map)
terra::ncell(dat_gebco())

#### Validity map
# TO DO


###########################
###########################
#### Simulated paths (~40 s)

file_paths <- "data/patter/input/paths.qs"
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
    select(sim_id = virt_fish, timestep, timestamp, x, y, region = regions) |> 
    as.data.table()
  
  #### Save 
  qs::qsave(paths, file_paths)
  toc()
  
} else {
  paths <- qs::qread(file_paths)
}


###########################
###########################
#### Observations

#### Define detection data (detections & receiver coordinates)
# Define receiver coordinates (UTM)
rxy <- 
  cbind(moorings$deploy_lon, moorings$deploy_lat) |> 
  terra::vect(crs = terra::crs(champlain)) |> 
  terra::project(epsg_utm) |>
  terra::crds()
points(rxy[, 1], rxy[, 2], col = "red")
stopifnot(all(!is.na(terra::extract(map, rxy)$map_value)))

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
         # * This is refined later
         receiver_gamma = 5000.0, 
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
  select(sim_id, 
         timestamp = "detection_timestamp_utc",
         receiver_lon = deploy_lon, 
         receiver_lat = deploy_lat) |> 
  # Define sim_id as integer
  mutate(sim_id = stringr::str_replace(sim_id, "sim_", ""), 
         sim_id = as.integer(sim_id)) |>
  mutate(receiver_key = paste(receiver_lon, receiver_lat)) |> 
  mutate(receiver_id = moorings$receiver_id[match(receiver_key, moorings$receiver_key)]) |>
  select(sim_id, timestamp, receiver_id) |> 
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

###########################
###########################
#### Metadata

#### Define metadata
metadata <-
  metadata |> 
  # Define sim_id as integer
  mutate(sim_id = stringr::str_replace(sim_id, "sim_", ""), 
         sim_id = as.integer(sim_id)) |>
  as.data.table()

#### Collect individual IDs (1:100)
ids <- sort(unique(metadata$sim_id))


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
    path <- paths[sim_id == id, ]
    timeline <- assemble_timeline(list(data.table(timestamp = start), path),
                                  .step = "2 mins", 
                                  .trim = FALSE)
    timeline
}) 
names(timelines) <- ids


###########################
###########################
#### Save datasets

#### Write datasets
terra::writeRaster(map, here_input("map.tif"), overwrite = TRUE)
qs::qsave(map_len, here_input("map_len.qs"))
qs::qsave(timelines, here_input("timelines.qs"))
qs::qsave(moorings, here_input("moorings.qs"))
qs::qsave(detections, here_input("detections.qs"))
qs::qsave(metadata, here_input("metadata.qs"))

#### Check sizes 
file.size(here_input("map.tif")) / 1e6 # MB

#### End of code. 
###########################
###########################