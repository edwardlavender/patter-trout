###########################
###########################
#### patter-exploration

#### Aims
# 1) Examine simulated datasets for patter

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
Sys.setenv("JULIA_SESSION" = FALSE)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(dv)
library(JuliaCall)
library(lubridate)
library(patter)
library(tictoc)
library(truncdist)
dv::src()

#### Load data 
map         <- terra::rast(here_input("map.tif"))
map_len     <- qs::qread(here_input("map_len.qs"))
start       <- qs::qread(here_input("start.qs"))
paths       <- qs::qread(here_input("paths.qs"))
moorings    <- qs::qread(here_input("moorings.qs"))
detections  <- qs::qread(here_input("detections.qs"))
metadata    <- qs::qread(here_input("metadata.qs"))


###########################
###########################
#### Examine movements

#### Movement simulation
# Paths were simulated via glatos::crw_in_polygon() which calls glatos::crw()
# * https://github.com/ocean-tracking-network/glatos/blob/main/R/sim-crw_in_polygon.r
# * https://github.com/ocean-tracking-network/glatos/blob/main/R/simutil-crw.r
# Step length in {trms_time} given {mvt_speed} is {mvt_speed} (m/s) * {trms_time} (s)
# heading_{t = 1} = Uniform(0, 360) 
# heading_{t > 1} = Normal(0, meta$theta) + cumsum(heading{t = 1, ..., t}) 
# Angles updated every {mvt_time} (s) / {trms_time} (s) time steps
# (Internally glatos enlarges the SD if the simulation steps outside the polygon!)

#### Compute path metrics
# `complete_simulated_transmissions_regions.rds` contains the 'full' path @ resolution of 127 s
# We can use this to conveniently check the distribution of step lengths & turning angles
# > Each simulated path was 5000 steps
# > Each step was 500 m
# > During a step, velocity was set to different values
# > Hence, steps lasted different durations (500 / velocity)
outfile <- here_input("path-metrics.qs")
if (!file.exists(outfile)) {
  
  # ~1.25 mins with 12 cl
  paths_metrics <- cl_lapply(unique(paths$sim_id), .cl = 12L, .fun = function(id) {
    
    path_metrics <- 
      paths |> 
      lazy_dt() |> 
      filter(sim_id == id) |> 
      select(sim_id, date = timestamp, x, y) |> 
      as.data.frame() |> 
      bayesmove::prep_data(coord.names = c("x", "y"), id = "sim_id") |> 
      select(sim_id, timestamp = date, x, y, step, angle, dt) |> 
      as.data.table()
    
    if (FALSE) {
      unique(diff(path_metrics$date)) # time stamps evenly spaced
      unique(path_metrics$step)       # step lengths not exactly equal 
      hist(path_metrics$step)         # but most step lengths are correct
      unique(path_metrics$angle)      # angles (radians)
    }
    
    path_metrics
    
  }) |> rbindlist()
  
  qs::qsave(path_metrics, outfile)
  
} else {
  paths_metrics <- qs::qread(outfile)
}

#### Examine step lengths
# Examine step lengths by individual
# * id 1: 12.7
# * id 2: 63.5
# * id 7: 114
paths_metrics |> 
  group_by(sim_id) |> 
  summarise(step = max(step, na.rm = TRUE)) |> 
  arrange(as.integer(sim_id))
# Estimate mobility 
steps <- paths_metrics$step[!is.na(paths_metrics$step)]
max(steps)
# 114.355 
# Estimate approximately suitable model across all individuals
spar <- run(file = here_input("step-fitdistr.qs"), 
            expr = {
              # ~2 mins
              MASS::fitdistr(steps, "gamma")
            }, 
            read = qs::qread, 
            write = qs::qsave)
step_shape <- spar$estimate["shape"] |> as.numeric()
step_scale <- 1 / spar$estimate["rate"] |> as.numeric()
mobility   <- 115.0
# shape           rate    
# 1.494954e+00   5.714505e-02 
# (6.233589e-04) (2.822856e-05)
hist(steps, prob = TRUE)
curve(dgamma(x, shape = step_shape, scale = step_scale), 
      lwd = 3, add = TRUE)
# Convergence
# > The averaged gamma model works for individual 1 (12.7 m step length)
# > But fails for individual 2 (63.5 m step length)
# > A more flexible model is required to handle the range of movements across all individuals

#### Compute validity map from mobility
vmap <- patter:::spatVmap(.map = map, .mobility = mobility, .plot = TRUE)
terra::writeRaster(vmap, here_input("vmap.tif"), overwrite = TRUE)

#### Examine turning angles
# This is a quick way of generating an approximately suitable model across all individuals
angles     <- paths_metrics$angle[!is.na(paths_metrics$angle)]
apar       <- MASS::fitdistr(angles, "normal")
angle_mean <- apar$estimate["mean"] |> as.numeric()
angle_sd   <- apar$estimate["sd"] |> as.numeric()
# mean             sd      
# -4.278476e-05    1.362498e-01 
# ( 4.416044e-05) ( 3.122614e-05)
hist(angles, prob = TRUE)
curve(dnorm(x, mean = angle_mean, sd = angle_sd), 
      lwd = 3, add = TRUE)

#### Movement model dimensions
png(here_fig("model-move.png"), height = 5, width = 10, units = "in", res = 600)
pp <- par(mfrow = c(1, 2))
# Step lengths
hist(steps, prob = TRUE, 
     xlab = "Step length (m)", ylab = "Density")
curve(dtrunc(x, spec = "gamma", a = 0, b = mobility, shape = step_shape, scale = 60), 
      from = 0, to = mobility, n = 1e3L, lwd = 3, add = TRUE)
# Turning angles
hist(angles, prob = TRUE, 
     xlab = "Turning angle (rad)", ylab = "Density")
curve(dnorm(x, 0, 0.4),
      from = -pi - 0.1, to = pi + 0.1, n = 1e3L, lwd = 3, add = TRUE)
par(pp)
dev.off()


###########################
###########################
#### Examine observations

#### Examine detection range (~2 s)
# > Compute the detection distances (ddists) between individuals & receivers at moment of detection

ddists <- 
  cl_lapply(split(paths, paths$sim_id), function(path) {
  
  # path <- split(paths, paths$sim_id)[[1]]
    
  # Combine detections with receiver coordinates 
  ddist <- merge(detections, moorings, by = "receiver_id")
  
  # At the moment of detection, add individual location
  ddist <- merge(ddist, path, by = "timestamp")
  
  # Compute distances between individual and receiver
  ddist |>
    select("timestamp", "receiver_x", "receiver_y", "x", "y") |> 
    mutate(dist = terra::distance(cbind(ddist$receiver_x, ddist$receiver_y),
                                  cbind(ddist$x, ddist$y), pairwise = TRUE, lonlat = FALSE)) |>
    as.data.table()
  
}) |> rbindlist()

# The maximum detection range in the simulated study system:
max(ddists$dist)
# 3236.828
receiver_gamma <- 3240

#### Examine movement observations
# What can we learn about movements from the observations? 
# Is there a case to increase information in the prior? 
if (requireNamespace("flapper", quietly = TRUE)) {
  detections[, fct := sim_id]
  msp <- sp::SpatialPoints(moorings[, c("receiver_x", "receiver_y")], sp::CRS(terra::crs(map)))
  msp <- sp::SpatialPointsDataFrame(msp, data.frame(receiver_id = moorings$receiver_id))
  mvt <- flapper::get_mvt_mobility_from_acoustics(data = detections, 
                                                  fct = "sim_id", 
                                                  moorings = msp, 
                                                  detection_range = receiver_gamma, 
                                                  step = 120,
                                                  transmission_interval = 127)
}


###########################
###########################
#### Save outputs

parameters <- 
  list(model_move = list(mobility = mobility, 
                         step = list(shape = step_shape, scale = step_scale), 
                         angle = list(mean = angle_mean, sd = angle_sd)), 
       model_obs = list(receiver_gamma = receiver_gamma))

qs::qsave(parameters, here_input("parameters.qs"))


#### End of code. 
###########################
###########################