###########################
###########################
#### patter-setup-real.R

#### Aims
# 1) Sets up real-world data for analysis with patter

#### Prerequisites
# 1) 


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
library(lubridate)
library(truncdist)
library(tictoc)
dv::src()

#### Load data 
champlain  <- terra::vect(here_data("ChamplainRegionsGrouped/ChamplainRegionsGrouped.shp"))
map        <- terra::rast(here_input("map.tif"))
moorings   <- readRDS(here_data("OriginalReceiverSummary_2013-2017.rds"))
detections <- readRDS(here_data("lkt_detections_2013-2017.rds"))


###########################
###########################
#### Prepare datasets

#### Prepare moorings

head(moorings)

# Define receiver coordinates (UTM)
rxy <- 
  cbind(moorings$deploy_lon, moorings$deploy_lat) |> 
  terra::vect(crs = terra::crs(champlain)) |> 
  terra::project("EPSG:3175") |>
  terra::crds()
stopifnot(all(!is.na(terra::extract(map, rxy)$map_value)))

# Process moorings
moorings <- 
  moorings |> 
  mutate(receiver_id = row_number(),
         receiver_sn = as.integer(as.character(receiver_sn)),
         receiver_start = as.POSIXct(paste0(deploy_date_time, "00:00:00"), tz = "UTC"), 
         receiver_end = as.POSIXct(paste0(recover_date_time, "00:00:00"), tz = "UTC"), 
         receiver_int = lubridate::interval(receiver_start, receiver_end),
         receiver_x = rxy[, 1],
         receiver_y = rxy[, 2]) |> 
  select(receiver_id, receiver_sn, receiver_start, receiver_end, receiver_int, receiver_x, receiver_y) |>
  as.data.frame()
         
# Add detection probability parameters
# * Implemented below

#### Prepare detections

# Examine selected columns
nrow(detections)
head(detections)
table(detections$passFilter)

# Process detections
detections <- 
  detections |> 
  mutate(individual_id = as.integer(as.character(animal_id)), 
         timestamp = as.POSIXct(detection_timestamp_utc, tz = "UTC"),
         receiver_id = NA_integer_, 
         receiver_sn = as.integer(as.character(receiver_sn))) |>
  select(individual_id, 
         timestamp,
         receiver_sn) |> 
  as.data.table()

# Order detections by duration
# * This improves speed during algorithm testing 
# * (NB: the code below works because duration is unique to each individual)
detections <- 
  detections |>
  group_by(individual_id) |> 
  mutate(duration = as.numeric(difftime(max(timestamp), min(timestamp), units = "days"))) |> 
  arrange(duration, timestamp) |> 
  mutate(-duration) |>
  as.data.table()

# Define receiver_id (~3 s)
# * match using receiver_sn and the time stamps
for (i in 1:nrow(moorings)) {
  detections[receiver_sn == moorings$receiver_sn[i] & 
               timestamp %within% moorings$receiver_int[i], receiver_id := moorings$receiver_id[i]]
}
table(is.na(detections$receiver_id))

# Drop 'extra' detections
detections <- detections[!is.na(receiver_id), ]

#### Clean up

# Clean up moorings 

head(moorings)
moorings <- 
  moorings |> 
  select(receiver_id, receiver_start, receiver_end, receiver_x, receiver_y) |> 
  as.data.table()

# Clean up detections

head(detections)

detections <-
  detections |> 
  select(individual_id, timestamp, receiver_id) |> 
  as.data.table()


###########################
###########################
#### Prepare parameters

#### Detection probability model
curve(plogis(1.8 + -0.005 * x), from = 0, to = 6000, ylim = c(0, 1))
moorings[, receiver_alpha := 1.8]
moorings[, receiver_beta := -0.005]
moorings[, receiver_gamma := 6000]

#### Movement model 

# 0.1 m/s -> 12 m/2 min
# 0.2 m/s -> 24 m/2 min
# 0.3 m/s -> 36 m/2 min
# 0.4 m/s -> 48 m/2 min
# 0.9 m/ms -> 108 m/2 min

# Normal 
# * This distribution probably permits overly low step lengths 
# * But otherwise covers a broad range of possible cruising speeds
curve(dtrunc(x, "norm", a = 0, b = 108, 24, 20), from = 0, to = 108)

# Gamma
# * Reduce shape to shift to left (scale-dependent)
# * Reduce rate to widen distribution
curve(dtrunc(x, "gamma", a = 0, b = 108, 5, 0.25), from = 0, to = 108)
curve(dtrunc(x, "gamma", a = 0, b = 108, 3, 0.15), from = 0, to = 108) # best-guess
curve(dtrunc(x, "gamma", a = 0, b = 108, 3, 0.1), from = 0, to = 108)

# Log normal
# * Reduce sdlog to broaden distribution
# * This parameterisation peaks too early
curve(dtrunc(x, "lnorm", a = 0, b = 108, 3, 1), from = 0, to = 108)

# Cauchy
# * This distribution also permits overly low step lengths
# * But has a longer tail to the right
# * Increase scale to widen distribution 
curve(dtrunc(x, "cauchy", a = 0, b = 108, 20, 10), from = 0, to = 108)

#### Validity map
vmap <- patter:::spatVmap(.map = map, .mobility = 108, .plot = TRUE)
terra::writeRaster(vmap, here_input_real("vmap.tif"), overwrite = TRUE)


###########################
###########################
#### Write outputs

qs::qsave(moorings, here_input_real("moorings.qs"))
qs::qsave(detections, here_input_real("detections.qs"))


#### End of code. 
###########################
###########################