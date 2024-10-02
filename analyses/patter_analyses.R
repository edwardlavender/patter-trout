###########################
###########################
#### patter-analyses.R

#### Aims
# 1) Implement patter algorithms 

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
library(JuliaCall)
library(lubridate)
library(patter)
library(tictoc)
# TO DO
# (Temporary)
dv::src()

#### Load data 
champlain  <- terra::vect("data/ChamplainRegions.shp")
moorings   <- readRDS("data/moorings.rds")
detections <- readRDS("data/simulated_detections.rds")
metadata   <- readRDS("data/simulations_metadata.rds")

#### Julia connection
julia_connect()
set_seed()
julia_command('include("Julia/src/model-movement.jl");')
julia_command('include("Julia/src/model-observation.jl")')


###########################
###########################
#### Prepare data (generic)

#### Define study area
champlain$region <- champlain$GNIS_NAME
champlain$land   <- as.numeric(1)
regions <- as_SpatRaster(champlain, .simplify = 0.001, .utm = "EPSG:32618",
                         .field = "region", .res = 10, .plot = TRUE)
maps    <- as_SpatRaster(champlain, .simplify = 0.001, .utm = "EPSG:32618",
                         .field = "land", .res = 10, .plot = TRUE)
map <- maps$SpatRaster
ydist <- terra::ymax(map) - terra::ymin(map)
set_map(map)

#### Define simulated paths 
file_paths <- "data/patter/paths.qs"
if (!file.exists(file_paths)) {
  # paths      <- readRDS("data/simulated_positions_July2023.rds")
  paths      <- readRDS("data/complete_simulated_transmissions_regions.rds")
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
    mutate(x = sxy[, 1], y = sxy[, 2]) |>
    select(sim_id = virt_fish, x, y, region = regions) |> 
    group_by(sim_id) |> 
    mutate(timestep = row_number()) |> 
    as.data.table()
  qs::qsave(paths, file_paths)
} else {
  paths <- qs::qread(file_paths)
}

#### Define detection data (detections & receiver coordinates)
# Define receiver coordinates (UTM)
rxy <- 
  cbind(moorings$deploy_lon, moorings$deploy_lat) |> 
  terra::vect(crs = "WGS84") |> 
  terra::project("EPSG:32618") |>
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


###########################
###########################
#### Run algorithms for a selected individual

#### Define individual & associated datasets
id   <- 1
sid  <- paste0("sim_", id)
path <- paths[sim_id == id, ]
dets <- detections[sim_id == sid, ]
meta <- metadata[sim_id == sid, ]
stopifnot(nrow(meta) == 1L)

#### Define study timeline

# The timeline is individual specific
# Start time: "2022-01-01 00:00:00"
# Simulated tracks comprised 5000 steps
# Each step comprised 500 m
# But velocity set to different values in transmit_along_path()
# I.e., step lengths differed in duration
# Transmissions were generated along the paths every 120 s

# Calculate the step duration
speed    <- meta$velocity      # speed @ each time step (m/s)
distance <- 500                # distance moved at each time step (m)
time     <- distance / speed   # duration of each time step (s)

# Define simulation start and end times
start    <- as.POSIXct("2022-01-01 00:00:00" , tz = "UTC")
end      <- max(seq(start, length.out = 5000, by = paste(time, "secs")))

# Define timeline 
# * The timeline is based on a transmission interval of 2 mins
# * The step length (m) in 2 mins given {speed} is speed * 120
timeline <- seq.POSIXt(start, end, by = "2 mins")

# Validate that the timeline spans the detection time series for the ID
stopifnot(interval(min(dets$timestamp), max(dets$timestamp)) %within% interval(min(timeline), max(timeline)))

#### Define movement model
# Paths were simulated via glatos::crw_in_polygon() which calls glatos::crw()
# * https://github.com/ocean-tracking-network/glatos/blob/main/R/sim-crw_in_polygon.r
# * https://github.com/ocean-tracking-network/glatos/blob/main/R/simutil-crw.r
# There were 500 steps
# Step length in 2 mins given {speed} is speed (m/s) * 120 (s)
# heading_{t = 1} = Uniform(0, 360) 
# heading_{t > 1} = Normal(0, meta$theta) + cumsum(heading{t = 1, ..., t})
# But internally glatos enlarges the SD if the simulation steps outside the polygon! 
state             <- "StateXYD"
model_move_length <- speed * 120
model_move        <- move_xyd(length = model_move_length, theta = meta$theta)

#### Define observation model
model_obs   <- c("ModelObsAcousticLogisTrunc", "ModelObsAcousticContainer")
acoustics   <- assemble_acoustics(.timeline = timeline, .acoustics = dets, .moorings = moorings)
containers  <- assemble_acoustics_containers(.acoustics = acoustics, 
                                            .direction = "forward", 
                                            .mobility = distance, 
                                            .threshold = ydist)
yobs        <- list(ModelObsAcousticLogisTrunc = acoustics, 
                    ModelObsAcousticContainer = containers)

#### Implement particle filter
# Define filter arguments 
args <- list(.map = map, 
             .timeline = timeline, 
             .state = state, 
             .xinit = NULL, 
             .xinit_pars = list(mobility = model_move_length),
             .yobs = yobs, 
             .model_obs = model_obs, 
             .model_move = model_move, 
             .n_particle = 25000L, 
             .direction = "forward"
             )
# Run forward filter 
# * NB: Setting observations is slow (1 min)
fwd <- do.call(pf_filter, args, quote = TRUE)
# Backward filter
# * Update containers & yobs for .direction = "backward", then run filter
containers      <- assemble_acoustics_containers(.acoustics = acoustics, 
                                                 .direction = "backward", 
                                                 .mobility = distance, 
                                                 .threshold = ydist)
yobs            <- list(ModelObsAcousticLogisTrunc = acoustics, 
                        ModelObsAcousticContainer = containers)
args$.yobs      <- yobs
args$.direction <- "backward"
do.call(pf_filter, args, quote = TRUE)

#### Run smoother 
# Run smoother
out_smo <- pf_smoother_two_filter()

#### Residency analysis
# True residency in each region
path_res <- 
  path |> 
  group_by(region) |> 
  summarise(n = n()) |> 
  mutate(residency = (n / sum(n)) * 100, 
         estimate = "path", 
         sim_id = id) |> 
  as.data.table()
# Particle residency estimates
part_res <- 
  out_smo |> 
  mutate(region = terra::extract(map, cbind(x, y))) |>
  group_by(region) |> 
  mutate(residency = (n / sum(n)) * 100, 
         estimate = "smoother", 
         sim_id = id) |> 
  as.data.table()
# Collect data
res <- rbind(path_res, part_res)
qs::qsave(res, "data/patter/qresidency/", paste0(id, ".qs"))


###########################
###########################
#### Synthesis

# Read residency data for each individual 
residency <- 
  lapply(unique(paths$sim_id), function(id) {
  qs::qread("data/patter/qresidency.qs")
}) |> rbindlist()

# Visualise residency ~ individual, coloured by truth/algorithm
# 
# > TO DO


#### End of code. 
###########################
###########################
