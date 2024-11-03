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
library(tictoc)
dv::src()

#### Load data 
map       <- terra::rast(here_input("map.tif"))
moorings   <- qs::qread(here_input_real("moorings.qs"))
detections <- qs::qread(here_input_real("detections.qs"))


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
# (But the movement model needs to consider the effect of random tranmission)


#### End of code. 
###########################
###########################