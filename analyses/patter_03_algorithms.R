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
map_len     <- qs::qread(here_input("map_len.qs"))
timelines   <- qs::qread(here_input("timelines.qs"))
paths       <- qs::qread(here_input("paths.qs"))
metadata    <- qs::qread(here_input("metadata.qs"))
moorings    <- qs::qread(here_input("moorings.qs"))
detections  <- qs::qread(here_input("detections.qs"))
parameters  <- qs::qread(here_input("parameters.qs"))


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

#### Run workflow for each individual
tic()
cl_lapply(1:100L, function(id) {
  patter_workflow(id = id, 
                  map = NULL, map_len = map_len, 
                  timeline = timelines, paths = paths, 
                  moorings = moorings, detections = detections, 
                  metadata = metadata, 
                  parameters = parameters, interactive = FALSE)
  
})
toc()

#### Convergence record
# ID 1: 
# * Random walk with uniform heading: convergence failures @ time ~62172 even with many particles
# * Moderately correlated random walk: success after ~1.11 hours on 12 workers (5e4L)


#### End of code. 
###########################
###########################