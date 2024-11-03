###########################
###########################
#### patter-run-real.R

#### Aims
# 1) Run particle algorithms for real-world datasets

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
moorings    <- qs::qread(here_input("moorings.qs"))
detections  <- qs::qread(here_input("detections.qs"))


###########################
###########################
#### Run algorithms 

#### Julia setup 
julia_connect()
julia_source("./Julia/src/model-move.jl")
set_seed()
set_map(here_input("map.tif"))
set_vmap(.vmap = here_input_real("vmap.tif"))

#### Cleanup
if (FALSE) {
  old <- list.files(here_output_real(), pattern = "\\.qs$", recursive = TRUE, full.names = TRUE)
  file.remove(old)
}

#### Testing
test <- FALSE
if (test) {
  # For testing, use a short timeline
  timelines <- lapply(timelines, function(timeline) timeline[1:1000])
}

#### Timing & parallelisation
# TO DO

#### Run workflow for each individual
tic()
cl_lapply(c(1L, 2L, 7L), function(id) {
  patter_workflow(id = id, 
                  moorings = moorings, detections = detections, 
                  model_move = "real")
  
})
toc()


#### End of code. 
###########################
###########################