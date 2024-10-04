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
julia_command('include("Julia/src/model-observation.jl")')


###########################
###########################
#### Prepare data (generic)

#### Define study start time
start <- as.POSIXct("2022-01-01 00:00:00" , tz = "UTC")

#### Define study area
epsg_utm         <- "EPSG:3175"
champlain$region <- champlain$GNIS_NAME
champlain$land   <- as.numeric(1)
regions <- as_SpatRaster(champlain, .simplify = 0.001, .utm = epsg_utm,
                         .field = "region", .res = 10, .plot = TRUE)
maps    <- as_SpatRaster(champlain, .simplify = 0.001, .utm = epsg_utm,
                         .field = "land", .res = 10, .plot = TRUE)
map <- maps$SpatRaster
ydist <- terra::ymax(map) - terra::ymin(map)
set_map(map)

#### Define simulated paths (~40 s)
file_paths <- "data/patter/paths.qs"
if (!file.exists(file_paths)) {
  
  tic()
  
  # complete_simulated_transmissions_regions.rds contains the full trajectory
  # > I.e., the position every 127 s
  # simulated_positions_July2023.rds contains the simulated positions only
  # > These positions can be more sparsely spaced depending on velocity7
  # > glatos imposes a timeline afterwards! 

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
    mutate(timestamp = start + time, 
           x = sxy[, 1], y = sxy[, 2]) |>
    group_by(virt_fish) |> 
    mutate(timestep = row_number()) |> 
    select(sim_id = virt_fish, timestep, timestamp, x, y, region = regions) |> 
    as.data.table()
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
# I.e., For each individual, the duration of a step length is different 
# Transmissions were generated along the paths every 120 + 7 s

# Calculate the step duration
trms_time    <- 120 # + 7                    # transmission time (every 127 s)
mvt_speed    <- meta$velocity              # speed @ each _movement_ time step (m/s)
mvt_distance <- 500                        # distance moved at each _movement_ time step (m)
mvt_time     <- mvt_distance / mvt_speed   # duration of each _movement_ time step (s)

# Define simulation start and end times
end      <- max(seq(start, length.out = 5000, by = paste(mvt_time, "secs")))

# Define timeline 
# * The timeline is based on a transmission interval of ~2 mins
# * (The step length (m) in ~2 mins given {mvt_speed} is {mvt_speed} * {trms_time})
# * (The angle is simulated every {mvt_time} s)
timeline <- seq.POSIXt(start, end, by = paste(trms_time, "secs"))

# Validate that the timeline spans the detection time series for the ID
stopifnot(interval(min(dets$timestamp), max(dets$timestamp)) %within% interval(min(timeline), max(timeline)))

#### Movement & timeline validation
# Check simulated step lengths and turning angles at 127 s resolution 
# * The path was simulated with 5000 steps, each of {mvt_time} secs
# * complete_simulated_transmissions_regions.rds contains the 'full' path @ resolution of 127 s
# * Step lengths should match {mvt_speed} (m/s) * {trms_time} (s)
mvt_len <- mvt_speed * trms_time
if (FALSE) {
  movements <- 
    path |> 
    lazy_dt() |> 
    select(sim_id, date = timestamp, x, y) |> 
    as.data.frame() |> 
    bayesmove::prep_data(coord.names = c("x", "y"), id = "sim_id")
  unique(diff(movements$date)) # time stamps evenly spaced
  unique(movements$step)       # step lengths not exactly equal 
  hist(movements$step)         # but most step lengths are correct
  unique(movements$angle)      # angles (radians)
}
# Check distance of individual to receiver @ moment of detection
# TO DO

#### Define movement model
# Paths were simulated via glatos::crw_in_polygon() which calls glatos::crw()
# * https://github.com/ocean-tracking-network/glatos/blob/main/R/sim-crw_in_polygon.r
# * https://github.com/ocean-tracking-network/glatos/blob/main/R/simutil-crw.r
# Step length in {trms_time} given {mvt_speed} is {mvt_speed} (m/s) * {trms_time} (s)
# heading_{t = 1} = Uniform(0, 360) 
# heading_{t > 1} = Normal(0, meta$theta) + cumsum(heading{t = 1, ..., t}) 
# Angles updated every {mvt_time} (s) / {trms_time} (s) time steps
# (Internally glatos enlarges the SD if the simulation steps outside the polygon!)
state              <- "StateXYD"
model_move         <- move_xyd(length = mvt_len, theta = meta$theta)
update_angle_every <- round(mvt_time / trms_time)
julia_assign("update_angle_every", update_angle_every)
julia_command('include("Julia/src/model-movement.jl");')

#### Define observation model
model_obs   <- c("ModelObsAcousticLogisTrunc", "ModelObsAcousticContainer")
acoustics   <- assemble_acoustics(.timeline = timeline, .acoustics = dets, .moorings = moorings)
containers  <- assemble_acoustics_containers(.acoustics = acoustics, 
                                            .direction = "forward", 
                                            .mobility = mvt_len, 
                                            .threshold = ydist)
yobs        <- list(ModelObsAcousticLogisTrunc = acoustics, 
                    ModelObsAcousticContainer = containers)

#### Implement particle filter
# Define filter arguments 
args <- list(.map = map, 
             .timeline = timeline, 
             .state = state, 
             .xinit = NULL, 
             .xinit_pars = list(mobility = mvt_len),
             .yobs = yobs, 
             .model_obs = model_obs, 
             .model_move = model_move, 
             .n_particle = 10000L, 
             .direction = "forward"
             )
# Run forward filter 
# * NB: Setting observations is slow (1 min)
# > Weights from filter (1 -> 208292) are zero at time 62367
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
