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
library(truncdist)
dv::src()

#### Load data 
if (!os_linux()) {
  map <- terra::rast(here_input("map.tif"))
  terra::plot(map)
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
julia_source("./Julia/src/model-move.jl")
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

#### (optional) Examine individual-specific data
# Simulated step lengths & turning angles
dist <- terra::distance(cbind(path$x, path$y), lonlat = FALSE, sequential = TRUE)
max(dist)
hist(dist)

#### Define movement model 
# Parameters
mobility <- parameters$model_move$mobility
sshape   <- parameters$model_move$step$shape
sscale   <- parameters$model_move$step$scale
amean    <- parameters$model_move$angle$mean
asd      <- parameters$model_move$angle$sd
# Model
# state      <- "StateXY"
# model_move <- move_xy(mobility = mobility, 
#                       dbn_length = glue("truncated(Gamma({sshape}, {sscale}), lower = 0.0, upper = {mobility})"), 
#                       dbn_angle = glue("truncated(Normal(0, 1), lower = -pi, upper = pi)"))
state      <- "StateXYD"
model_move <- move_xyd(mobility = mobility, 
                       dbn_length = glue("truncated(Gamma({sshape}, {sscale}), lower = 0.0, upper = {mobility})"), 
                       dbn_angle_delta = glue("Normal({0.0}, {0.4})"))
# Comments
# * Overly high angle concentration leads to slow algorithm runs (many trials required to avoid jumping on land)

#### Visualise movement model 
# Visualise model components 
pp <- par(mfrow = c(1, 2))
curve(dtrunc(x, spec = "gamma", a = 0, b = mobility, shape = sshape, scale = sscale), 
      from = 0, to = mobility, n = 1e3L,
      xlab = "Step length (m)", ylab = "Density")
abline(v= max(dist), col = "red")
curve(dnorm(x, 0, 0.25),
      from = -pi - 0.1, to = pi + 0.1, n = 1e3L,
      xlab = "Turning angle (rad)", ylab = "Density")
par(pp)
# Visualise realisations of the movement model
length(timeline)
pos <- 1:50000
if (!os_linux()) {
  sim_path_walk(.map = map, 
                .timeline = timeline[pos],
                .state = state, 
                .model_move = model_move, 
                .n_path = 6L, .one_page = TRUE)
}
# Compare to simulated paths
if (!os_linux()) {
  pp <- par(mfrow = c(2, 2))
  lapply(1:4, function(id) {
    path <- paths[sim_id == id, ]
    path <- path[pos, ]
    terra::plot(map)
    patter:::add_sp_path(path$x, path$y)
  })
  par(pp)
}

#### Define observation models
# Assemble acoustics (0, 1)
moorings[, receiver_gamma := parameters$model_obs$receiver_gamma]
acoustics   <- assemble_acoustics(.timeline   = timeline, 
                                  .detections = dets, 
                                  .moorings   = moorings)
# Record detections (1)
# * This is useful for checking convergence failures
detections  <- acoustics[obs == 1L, ]
detections[, timestep := 0L]
for (i in 1:nrow(detections)) {
  detections[i, timestep := which(timeline == timestamp)]
}
# Assemble containers
containers  <- assemble_acoustics_containers(.timeline  = timeline, 
                                             .acoustics = acoustics, 
                                             .mobility  = mobility, 
                                             .threshold = map_len)
# Define yobs (forward & backward)
yobs_fwd <- list(ModelObsAcousticLogisTrunc = copy(acoustics),
                 ModelObsAcousticContainer  = copy(containers$forward))
yobs_bwd <- list(ModelObsAcousticLogisTrunc = copy(acoustics),
                 ModelObsAcousticContainer  = copy(containers$backward))

#### Define initial states
xinit_fwd <- xinit_bwd <- NULL


###########################
###########################
#### Algorithm runs 

#### Define filter arguments 
args <- list(.timeline = timeline, 
             .state = state, 
             .xinit = xinit_fwd, 
             # .yobs = yobs_fwd,
             .model_move = model_move, 
             .n_particle = 5e4L,
             .n_iter = 1L,
             .direction = "forward"
             )

#### Run forward filter 
# Setting initial observations is slow (~1 min)
# Set yobs to missing in `args` to speed up multiple runs 
fwd <- do.call(pf_filter, args, quote = TRUE)
# qs::qsave(fwd, here_output("tmp", "tmp.qs"))

#### Run backward filter
args$.xinit     <- xinit_bwd
args$.yobs      <- yobs_bwd
args$.direction <- "backward"
do.call(pf_filter, args, quote = TRUE)

#### Run smoother 
smo <- pf_smoother_two_filter()
qs::qsave(smo, here_output("particles", paste0(id, ".qs")))


#### End of code. 
###########################
###########################