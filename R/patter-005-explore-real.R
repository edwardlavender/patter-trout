###########################
###########################
#### patter-explore-real.R

#### Aims
# 1) Explore real-world datasets

#### Prerequisites
# 1) NA


###########################
###########################
#### Setup 

#### Wipe workspace 
rm(list = ls())
# try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
Sys.setenv("JULIA_SESSION" = FALSE)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tictoc)
dv::src()

#### Load data 
map         <- terra::rast(here_input("map.tif"))
moorings    <- qs::qread(here_input_real("moorings.qs"))
detections  <- qs::qread(here_input_real("detections.qs"))


###########################
###########################
#### Detections

#### Examine timeline duration
# Compute durations
duration <- 
  detections |> 
  group_by(individual_id) |> 
  summarise(duration = as.numeric(difftime(max(timestamp), min(timestamp), units = "days"))) |> 
  arrange(duration) |> 
  as.data.table()
table(table(duration$duration))
# Examine durations
duration
hist(duration$duration)

#### Examine memory requirements
# Compute pf_particles object size for largest detection time series
# > Build example states and diagnostics data.tables
# > Check sizes
tic()
np <- 500L
nt <- ceiling(60/2 * 24 * 1027.806516)  # number of time steps
nr <- nt * np                           # number of rows: 370,010,346
states <- data.table(path_id = rep(1L:np, each = nt),
                     timestep = rep(1L:nt, np), 
                     timestamp = rep(seq.POSIXt(as.POSIXct("2016-01-01"), by = "2 mins", length.out = nt), np), 
                     map_value = runif(nr), 
                     x = runif(nr), 
                     y = runif(nr), 
                     angle = runif(nr))
diagnostics <- data.table(timestep = 1L:nt,
                          timestamp = seq.POSIXt(as.POSIXct("2016-01-01"), by = "2 mins", length.out = nt), 
                          ess = runif(nt), 
                          maxlp = runif(nt))
fwd <- list(states = states, diagnostics = diagnostics, convergence = TRUE, trials = 1L)
pryr::object_size(states)      # 16.28 GB
pryr::object_size(diagnostics) # 17.76 MB
pryr::object_size(fwd)         # 16.30 GB
toc()
# Results
# > A single output from the filter may be ~16 GB (x 2: in R and Julia)
# > We also require outputs from the backward filter and the smoother
# > For the longer time series, batching is required


###########################
###########################
#### Mobility 

#### Estimation method
# We compute three speed metrics & for each we compute the min/mean/max value
# Metrics:
# * speed_min: minimum speed, based on movement between nearest container edges 
# * speed_avg: average speed, based on movement between receivers
# * speed_max: maximum speed, based on movement between farthest container edge
# Interpretation:
# * min speeds are 'inappropriately low' b/c indirect routes between receivers are taken
# * maximum values for each metric can inform us about mobility 

#### Estimate mobility (~4 s)
if (requireNamespace("flapper", quietly = TRUE)) {
  tic()
  detections[, fct := individual_id]
  msp <- sp::SpatialPoints(moorings[, c("receiver_x", "receiver_y")], sp::CRS(terra::crs(map)))
  msp <- sp::SpatialPointsDataFrame(msp, data.frame(receiver_id = moorings$receiver_id))
  mvt <- flapper::get_mvt_mobility_from_acoustics(data = detections, 
                                                  fct = "individual_id", 
                                                  moorings = msp, 
                                                  detection_range = moorings$receiver_gamma[1], 
                                                  calc_distance = "lcp", 
                                                  bathy = raster::raster(map),
                                                  step = 120,
                                                  transmission_interval = 160)
  toc()
}

#### Results
# --------------------------------------
#   Estimates (m/s)-----------------------
#   variable min mean  max
# 1 speed_min_ms   0 0.03 0.40
# 2 speed_avg_ms   0 0.10 0.91
# 3 speed_max_ms   0 0.18 1.65
# --------------------------------------
#   Estimates (m/step)--------------------
#   variable  min  mean    max
# 1 speed_min_mstep 0.00  3.95  48.12
# 2 speed_avg_mstep 0.07 12.54 109.43
# 3 speed_max_mstep 0.12 21.14 198.51
# --------------------------------------

#### Conclusions 
# Even if we assume 'minimum distance' movements, maximum speeds of 48 m/2 min are apparent.
# (But travelled  speeds are likely to be higher than this value 
# ... as this assumes the minimum possible travel distance)
# If we assume larger movements, speeds may be up to 109 m/2 min.
# This is a reasonable middle-of-the-road estimate for maximum movement speeds. 
# In the most extreme case (unlikely), movement speeds up to 198 m/2 min are apparent.
# (This is the upper bound for movement speeds suggested by the data.)
# These results are consistent for different transmission intervals
# (But the movement model needs to consider the effect of random transmission)


#### End of code. 
###########################
###########################