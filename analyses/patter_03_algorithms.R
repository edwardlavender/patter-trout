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
library(glue)
library(JuliaCall)
library(lubridate)
library(patter)
library(tictoc)
dv::src()

#### Load data 
#  map <- terra::rast(here_input("map.tif"))
map_len     <- qs::qread(here_input("map_len.qs"))
timelines   <- qs::qread(here_input("timelines.qs"))
paths       <- qs::qread(here_input("paths.qs"))
metadata    <- qs::qread(here_input("metadata.qs"))
moorings    <- qs::qread(here_input("moorings.qs"))
detections  <- qs::qread(here_input("detections.qs"))
parameters  <- qs::qread(here_input("parameters.qs"))

#### Julia setup 
julia_connect()
set_seed()
set_map("./data/patter/input/map.tif")


###########################
###########################
#### Prepare algorithms

#### Define individual & associated datasets
id   <- 1
path <- paths[sim_id == id, ]
dets <- detections[sim_id == id, ]
meta <- metadata[sim_id == id, ]
stopifnot(nrow(meta) == 1L)

#### Define timeline (individual-specific)
timeline <- timelines[[id]]

#### Define movement model 
# Parameters
mobility <- parameters$model_move$mobility
sshape   <- parameters$model_move$step$shape
sscale   <- parameters$model_move$step$scale
amean    <- parameters$model_move$angle$mean
asd      <- parameters$model_move$angle$sd
# Model
state      <- "StateXY"
model_move <- move_xy(mobility = mobility, 
                      dbn_length = glue("truncated(Gamma({sshape}, {sscale}), lower = 0.0, upper = {mobility})"), 
                      dbn_angle = glue("Normal({amean}, {asd})"))

#### Define observation models
# Assemble acoustics (0, 1)
moorings[, receiver_gamma := parameters$model_obs$receiver_gamma]
acoustics   <- assemble_acoustics(.timeline   = timeline, 
                                  .detections = dets, 
                                  .moorings   = moorings)
# Assemble containers
containers  <- assemble_acoustics_containers(.timeline  = timeline, 
                                             .acoustics = acoustics, 
                                             .mobility  = mobility, 
                                             .threshold = map_len)
# Define yobs (forward & backward)
yobs_fwd <- list(ModelObsAcousticLogisTrunc = acoustics, 
                 ModelObsAcousticContainer  = containers$forward)
yobs_bwd <- list(ModelObsAcousticLogisTrunc = acoustics, 
                 ModelObsAcousticContainer  = containers$backward)


###########################
###########################
#### Algorithm runs 

#### Define filter arguments 
args <- list(.timeline = timeline, 
             .state = state, 
             .xinit = NULL, 
             .yobs = yobs_fwd,
             .model_move = model_move, 
             .n_particle = 1e4L, 
             .direction = "forward"
             )

#### Run forward filter 
# * NB: Setting observations is slow (1 min)
# > Weights from filter (1 -> 208292) are zero at time 62367
fwd <- do.call(pf_filter, args, quote = TRUE)

#### Run backward filter
args$.yobs      <- yobs_bwd
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
