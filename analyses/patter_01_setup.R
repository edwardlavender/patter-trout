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

#### Load data 
champlain  <- terra::vect("data/ChamplainRegions.shp")
moorings   <- readRDS("data/moorings.rds")
detections <- readRDS("data/simulated_detections.rds")
metadata   <- readRDS("data/simulations_metadata.rds")
dv::src()


###########################
###########################
#### Prepare datasets

#### Define study start time
start <- as.POSIXct("2022-01-01 00:00:00" , tz = "UTC")

#### Define study area (~2 s)
tic()
epsg_utm         <- "EPSG:3175"
champlain$region <- champlain$GNIS_NAME
champlain$land   <- as.numeric(1)
regions <- as_SpatRaster(champlain, .simplify = 0.001, .utm = epsg_utm,
                         .field = "region", .res = 10, .plot = TRUE)
maps    <- as_SpatRaster(champlain, .simplify = 0.001, .utm = epsg_utm,
                         .field = "land", .res = 10, .plot = TRUE)
map     <- maps$SpatRaster
map_len <- terra::ymax(map) - terra::ymin(map)
toc()

#### Define validity map
# TO DO

#### Define simulated paths (~40 s)
file_paths <- "data/patter/input/paths.qs"
if (!file.exists(file_paths)) {
  
  tic()
  
  # Read paths 
  # * `complete_simulated_transmissions_regions.rds` contains the full trajectory
  # > I.e., the position every 127 s
  # * `simulated_positions_July2023.rds` contains the simulated positions only
  # > These positions can be more sparsely spaced depending on velocity
  # > (glatos imposes a timeline afterwards)
  # paths <- readRDS("data/simulated_positions_July2023.rds")
  paths   <- readRDS("data/complete_simulated_transmissions_regions.rds")
  
  # Define path coordinates (UTM)
  sxy <- 
    paths$geometry |>
    sf::st_coordinates() |>
    terra::vect(crs = terra::crs(champlain)) |> 
    terra::project(terra::crs(map)) |> 
    terra::crds()
  
  # Define paths data.table
  paths <- 
    paths |> 
    as.data.table() |>
    mutate(timestamp = start + time, 
           x = sxy[, 1], y = sxy[, 2]) |>
    group_by(virt_fish) |> 
    mutate(timestep = row_number()) |> 
    select(sim_id = virt_fish, timestep, timestamp, x, y, region = regions) |> 
    as.data.table()
  
  # Save 
  qs::qsave(paths, file_paths)
  toc()
  
} else {
  paths <- qs::qread(file_paths)
}

#### Define detection data (detections & receiver coordinates)
# Define receiver coordinates (UTM)
rxy <- 
  cbind(moorings$deploy_lon, moorings$deploy_lat) |> 
  terra::vect(crs = terra::crs(champlain)) |> 
  terra::project(epsg_utm) |>
  terra::crds()
points(rxy[, 1], rxy[, 2], col = "red")
stopifnot(all(!is.na(terra::extract(map, rxy)$map_value)))
# Define moorings data.table
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
         # (untruncated model used for simulation)
         receiver_gamma = 5000.0, 
         receiver_key = paste(receiver_lon, receiver_lat)) |> 
  select(receiver_id, 
         receiver_start, receiver_end,
         receiver_x, receiver_y,
         receiver_alpha, receiver_beta, receiver_gamma,
         receiver_lon, receiver_lat, receiver_key) |> 
  as.data.table()
# Define detections 
detections <- 
  detections |> 
  select(sim_id, 
         timestamp = "detection_timestamp_utc",
         receiver_lon = deploy_lon, 
         receiver_lat = deploy_lat) |> 
  mutate(receiver_key = paste(receiver_lon, receiver_lat)) |> 
  mutate(receiver_id = moorings$receiver_id[match(receiver_key, moorings$receiver_key)]) |>
  as.data.table()
# Clean up 
moorings <- 
  moorings |> 
  select("receiver_id", 
         "receiver_start", "receiver_end", 
         "receiver_x", "receiver_y", 
         "receiver_alpha", "receiver_beta", "receiver_gamma") |> 
  as.data.table()
detections <- 
  detections |> 
  select(sim_id, timestamp, receiver_id) |> 
  as.data.table()
# Checks
stopifnot(all(!is.na(detections$receiver_id)))

#### Define metadata
metadata <- as.data.table(metadata)

#### Save datasets
terra::writeRaster(map, here_input("map.tif"))
qs::qsave(map_len, here_input("map_len.qs"))
qs::qsave(start, here_input("start.qs"))
qs::qsave(moorings, here_input("moorings.qs"))
qs::qsave(detections, here_input("detections.qs"))
qs::qsave(metadata, here_input("metadata.qs"))


#### End of code. 
###########################
###########################