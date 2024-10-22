###########################
###########################
#### patter-analyses

#### Aims
# 1) Run patter analyses

#### Prerequisites
# 1) patter-setup.R


###########################
###########################
#### Setup 

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

#### Load data 
champlain  <- terra::vect("data/ChamplainRegions.shp")
moorings   <- readRDS("data/moorings.rds")
detections <- readRDS("data/simulated_detections.rds")
metadata   <- readRDS("data/simulations_metadata.rds")
dv::src()

#### Julia setup 
julia_connect()
set_seed()
set_map()

###########################
###########################
#### Prepare algorithms


###########################
#### Individual datasets

#### Define individual & associated datasets
id   <- 1
sid  <- paste0("sim_", id)
path <- paths[sim_id == id, ]
dets <- detections[sim_id == sid, ]
meta <- metadata[sim_id == sid, ]
stopifnot(nrow(meta) == 1L)


###########################
#### Define study timeline

# The timeline is individual specific
# Start time: "2022-01-01 00:00:00"
# Simulated tracks comprised 5000 steps
# Each step comprised 500 m
# But velocity set to different values in transmit_along_path()
# I.e., For each individual, the duration of a step length is different 
# Transmissions were generated along the paths every 120 + 7 s

# Calculate the step duration
trms_time    <- 120 # + 7                  # transmission time (every 127 s)
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


###########################
#### Movement model

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


###########################
#### Define observation model

model_obs   <- c("ModelObsAcousticLogisTrunc", "ModelObsAcousticContainer")
acoustics   <- assemble_acoustics(.timeline = timeline, .acoustics = dets, .moorings = moorings)
containers  <- assemble_acoustics_containers(.acoustics = acoustics, 
                                            .direction = "forward", 
                                            .mobility = mvt_len, 
                                            .threshold = ydist)
yobs        <- list(ModelObsAcousticLogisTrunc = acoustics, 
                    ModelObsAcousticContainer = containers)


###########################
###########################
#### Algorithm runs 

#### Define filter arguments 
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


#### Run forward filter 
# * NB: Setting observations is slow (1 min)
# > Weights from filter (1 -> 208292) are zero at time 62367
fwd <- do.call(pf_filter, args, quote = TRUE)

#### Run backward filter
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
out_smo <- pf_smoother_two_filter()


###########################
###########################
#### Quick residency analysis

#### True residency in each region
path_res <- 
  path |> 
  group_by(region) |> 
  summarise(n = n()) |> 
  mutate(residency = (n / sum(n)) * 100, 
         estimate = "path", 
         sim_id = id) |> 
  as.data.table()

#### Particle residency estimates
part_res <- 
  out_smo |> 
  mutate(region = terra::extract(map, cbind(x, y))) |>
  group_by(region) |> 
  mutate(residency = (n / sum(n)) * 100, 
         estimate = "smoother", 
         sim_id = id) |> 
  as.data.table()

#### Collect data
res <- rbind(path_res, part_res)
qs::qsave(res, "data/patter/output/residency/qresidency/", paste0(id, ".qs"))


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
