###########################
###########################
#### patter-run-sim.R

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
Sys.setenv("JULIA_SESSION" = TRUE)
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
timelines   <- qs::qread(here_input_sim("timelines.qs"))
# paths     <- qs::qread(here_input_sim("paths.qs"))
# metadata  <- qs::qread(here_input_sim("metadata.qs"))
moorings    <- qs::qread(here_input_sim("moorings.qs"))
detections  <- qs::qread(here_input_sim("detections.qs"))


###########################
###########################
#### Run algorithms 

#### Julia setup 
julia_connect()
julia_source("./Julia/src/model-move.jl")
set_seed()
set_map(here_input("map.tif"))
set_vmap(.vmap = here_input("vmap.tif"))

#### Testing
test <- FALSE
if (test) {
  # For testing, use a short timeline
  timelines <- lapply(timelines, function(timeline) timeline[1:1000])
}

#### Timing & parallelisation
# Approximate timing for individual 1, forward filter on siam-linux20:
# * 15 threads: ~28 s vs. 12 threads MacBook ~10 s
# * 30 threads: ~25 s
# * 64 threads: 28 s
# * 120 threads: 53 s

#### Cleanup
if (FALSE) {
  old <- list.files(here_output_sim(), pattern = "\\.qs$", recursive = TRUE, full.names = TRUE)
  file.remove(old)
}

#### Run workflow for each individual
# To test convergence, use individuals 1 (12.7 m), 2 (63.5 m) and 7 (114 m)
# Repeat for each movement model ("low", "medium", "high")
tic()
cl_lapply(c(1L, 2L, 7L), function(id) {
  patter_workflow(id = id, 
                  timelines = timelines,
                  moorings = moorings, detections = detections, 
                  model_move = "sim-low")
  
})
toc()


#### End of code. 
###########################
###########################