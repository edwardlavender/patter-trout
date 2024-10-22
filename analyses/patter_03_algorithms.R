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
if (!os_linux()) {
  map <- terra::rast(here_input("map.tif"))
}
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
                      dbn_angle = glue("Uniform(-pi, pi)"))

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
             .n_particle = 2e4L, 
             .direction = "forward"
             )

#### Run forward filter 
fwd <- do.call(pf_filter, args, quote = TRUE)

#### Run backward filter
args$.yobs      <- yobs_bwd
args$.direction <- "backward"
do.call(pf_filter, args, quote = TRUE)

#### Run smoother 
smo <- pf_smoother_two_filter()
qs::qsave(smo, here_output("particles", paste0(id, ".qs")))


#### End of code. 
###########################
###########################
