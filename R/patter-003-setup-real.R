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
library(ggplot2)
library(lubridate)
library(truncdist)
library(tictoc)
dv::src()

#### Load data 
champlain  <- terra::vect(here_data("ChamplainRegionsGrouped/ChamplainRegionsGrouped.shp"))
map        <- terra::rast(here_input("map.tif"))
map_bbox   <- qs::qread(here_input("map-bbox.qs"))
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

# Define receiver_id (~4 s)
# * match using receiver_sn and the time stamps
tic()
for (i in 1:nrow(moorings)) {
  detections[receiver_sn == moorings$receiver_sn[i] & 
               timestamp %within% moorings$receiver_int[i], receiver_id := moorings$receiver_id[i]]
}
toc()
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
#### Compute container thresholds 

#### Method
# In assemble_acoustics_containers(), the .map argument 
# ... isn't supported on linux if JULIA_SESSION = "TRUE".
# Here, we prepare a data.table and use it to filter the output
# ... from assemble_acoustics_containers() on the fly. 

#### Compute max. dist. btwn receivers and study area edges
# Compute maximum distances
dist <- terra::distance(cbind(moorings$receiver_x, moorings$receiver_y), map_bbox, lonlat = FALSE)
dist <- apply(dist, 1, max)
# Compute container thresholds
cthresholds <- data.table(receiver_id = moorings$receiver_id, 
                          distance = dist)


###########################
###########################
#### Batch datasets

#### Method
# Batch the detection datasets to manage memory (see also explore-real.R)
# This should also mitigate potential convergence issues
# Split at the moment of a detection, roughly into one-month batches 
# Join batches detections (last time step, first time step)
# Then we can treat the time series as if they are from different individuals

#### Define individuals/months
detections <- 
  detections |> 
  group_by(individual_id) |> 
  mutate(unit_id = as.character(cut(timestamp, "months")), 
         unit_id = stringr::str_replace_all(unit_id, "-", "")) |> 
  select(individual_id, unit_id, timestamp, receiver_id) |>
  as.data.table()

#### Update detections 
detections <- 
  lapply(split(detections, detections$individual_id), function(d_id) {
  # d_id <- split(detections, detections$individual_id)[[1]]
  d_batch <- split(d_id, d_id$unit_id)
  n_batch <- length(d_batch)
  # (optional) Merge batches
  # * Merge batches if subsequent batches only contain a few additional observations
  if (n_batch > 1L) {
    for (i in n_batch:2L) {
      # Compute duration of current batch
      duration <- difftime(max(d_batch[[i]]$timestamp), 
                           min(d_batch[[i]]$timestamp), 
                           units = "weeks")
      # If duration is less than threshold (e.g., 1 week), merge with previous batch
      if (duration < 1) {
        # Combine datasets, using unit_id of preceeding batch
        d_batch[[i]][, unit_id := d_batch[[i - 1]]$unit_id[1]]
        d_batch[[i - 1]] <- rbind(d_batch[[i - 1]], d_batch[[i]])
        d_batch[[i]] <- NULL
      }
    }
    d_batch <- plyr::compact(d_batch)
  }
  n_batch <- length(d_batch)
  # Link batches
  # * The next batch should begin with the last observation from the previous batch
  # * This means we always start/stop filter @ a detection
  # * (i.e., in a 'happy' place) 
  if (n_batch > 1L) {
    for (i in 1:(n_batch - 1)) {
      last_row <- d_batch[[i]][.N, ]
      d_batch[[i + 1]] <- rbind(last_row, d_batch[[i + 1]])
    }
  }
  # Update unit_id
  # (use last unit_id as first one for batches 2:N is different)
  for (i in 1:n_batch) {
    d_batch[[i]][, unit_id := paste0(individual_id, "_", unit_id[.N])]
  }
  # Rejoin batches
  rbindlist(d_batch)
}) |> rbindlist()

#### Check code works
det_1 <- detections[individual_id == 26786, ]
det_1 <- split(det_1, det_1$unit_id)
length(det_1)
(det_1a <- det_1[[1]][.N, ])
(det_1b <- det_1[[2]][1, ])
det_1[[2]][1:3, ]
stopifnot(all.equal(det_1a$individual_id, det_1b$individual_id))
stopifnot(all.equal(det_1a$timestamp, det_1b$timestamp))
stopifnot(all.equal(det_1a$receiver_id, det_1b$receiver_id))

#### Check the number of observations & time range per batch
nobs <- 
  detections |> 
  group_by(unit_id) |> 
  summarise(n = n(), 
            duration = as.numeric(difftime(max(timestamp), min(timestamp)),
                                  units = "mins")) |> 
  arrange(n) |>
  as.data.table()
head(sort(nobs$n))
head(sort(nobs$duration)) 

#### Record mapping between individual_id and unit_id
detections_units <- 
  detections |> 
  select(individual_id, unit_id) |> 
  group_by(unit_id) |> 
  slice(1L) |> 
  ungroup() |> 
  group_by(individual_id) |>
  mutate(n_batch = n()) |> 
  ungroup() |>
  arrange(individual_id, unit_id) |>
  as.data.table()
# View(detections_units)

#### Set individual_id = unit_id (backwards compatibility)
# We implement the algorithms for each 'individual'
# To compute residency:
# - We iterate over individuals in detections_units
# - For each individual, we read the data for each batch
# - We process (shrink) the data for that batch
# - For batchs 2:N, we drop the first row to avoid douple counting 
# - see analysis-real.R (TO DO)
detections[, individual_id := unit_id]
detections[, unit_id := NULL]


###########################
###########################
#### Prepare detection parameters

#### Best-guess summer model 
curve(plogis(1.8 + -0.005 * x), from = 0, to = 6000, ylim = c(0, 1))

#### Best-guess winter model
curve(plogis(3 + -0.0025 * x), from = 0, to = 6000, col = "blue", add = TRUE)

#### Compromise model
curve(plogis(2.5 + -0.003 * x), from = 0, to = 6000, col = "darkgreen", add = TRUE)

#### ggplot representation
alpha          <- 2.5    # intercept
beta           <- -0.003 # rate of decline (larger values, nearer 0, increase steepness)
receiver_gamma <- 7000
ggplot(data.frame(x = c(0, receiver_gamma)), aes(x = x)) +
  stat_function(fun = function(x) {
    prob <- plogis(alpha + beta * x)
    prob[x > receiver_gamma] <- 0
    prob
  }) +
  ylim(0, 1) + 
  scale_x_continuous(breaks = seq(0, receiver_gamma, by = 500))  + 
  theme_bw()

#### Update moorings
moorings[, receiver_alpha := 2.25]
moorings[, receiver_beta := -0.0022]
moorings[, receiver_gamma := 7000]


###########################
###########################
#### Prepare movement parameters

#### Speeds
# 0.1 m/s -> 12 m/2 min,    18 m/3 min
# 0.2 m/s -> 24 m/2 min,    36 m/3 min
# 0.3 m/s -> 36 m/2 min,    54 m/3 min
# 0.4 m/s -> 48 m/2 min,    72 m/3 min
# 0.9 m/ms -> 108 m/2 min,  180 m/3 min

#### Define mobility 
mobility <- 200

#### Normal 
# This distribution probably permits overly low step lengths 
# But otherwise covers a broad range of possible cruising speeds
curve(dtrunc(x, "norm", a = 0, b = mobility, 24, 20), from = 0, to = mobility)

#### Gamma
# Reduce shape to shift to left (scale-dependent)
# Reduce rate to widen distribution
curve(dtrunc(x, "gamma", a = 0, b = mobility, 5, 0.25), from = 0, to = mobility)
curve(dtrunc(x, "gamma", a = 0, b = mobility, 3, 0.15), from = 0, to = mobility)

#### Gamma (ggplot2)
# Define example parameters
p1 <- c(2.8, 1 / 0.05)
# Shrink/expand distribution e.g., by 120 % while maintaining the same mode 
p2 <- gamma_rescale(p1[1], p1[2], fact = 1.2)
p3 <- gamma_rescale(p1[1], p1[2], fact = 0.8)
# Visualise Gamma distributions
ggplot(data.frame(x = c(0, mobility)), aes(x = x)) +
  # stat_dtruncgamma(mobility, shape = 4.5, scale = 1/0.1) + 
  # stat_dtruncgamma(mobility, shape = 3.2, scale = 1/0.06, col = "blue") + 
  stat_dtruncgamma(mobility, shape = 4.5, scale = 1/0.10, 
                   col = "grey", linetype = 2, linewidth = 1.5) + 
  stat_dtruncgamma(mobility, shape = p1[1], scale = p1[2], col = "black") + 
  stat_dtruncgamma(mobility, shape = p2[1], scale = p2[2], col = "green") + 
  stat_dtruncgamma(mobility, shape = p3[1], scale = p3[2], col = "red") + 
  scale_x_continuous(breaks = seq(0, mobility, by = 10))  + 
  theme_bw()

#### Log normal
# Reduce sdlog to broaden distribution
# This parameterisation peaks too early
curve(dtrunc(x, "lnorm", a = 0, b = mobility, 3, 1), from = 0, to = mobility)

#### Cauchy
# This distribution also permits overly low step lengths
# But has a longer tail to the right
# Increase scale to widen distribution 
curve(dtrunc(x, "cauchy", a = 0, b = mobility, 20, 10), from = 0, to = mobility)

#### Validity map
vmap <- patter:::spatVmap(.map = map, .mobility = mobility, .plot = TRUE)
terra::writeRaster(vmap, here_input_real("vmap.tif"), overwrite = TRUE)


###########################
###########################
#### Write outputs

qs::qsave(moorings, here_input_real("moorings.qs"))
qs::qsave(detections, here_input_real("detections.qs"))
qs::qsave(detections_units, here_input_real("detections-units.qs"))
qs::qsave(cthresholds, here_input_real("cthresholds.qs"))


#### End of code. 
###########################
###########################