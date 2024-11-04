###########################
###########################
#### patter-run-real.R

#### Aims
# 1) Run particle algorithms for real-world datasets

#### Prerequisites
# 1) patter-setup.R

#### TO DO
# Implement batching for big datasets (see setup-real.R)


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
moorings    <- qs::qread(here_input_real("moorings.qs"))
detections  <- qs::qread(here_input_real("detections.qs"))
cthresholds <- qs::qread(here_input_real("cthresholds.qs"))


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
# Pick individuals with shorter time series for testing
# ids <- as.integer(c(26790, 26793, 26795, 26804, 24320, 24390, 24346, 24377, 24326, 26785))
ids <- unique(detections$individual_id)

#### Timing & parallelisation
# TO DO

#### Run workflow for each individual
length(ids)
log.txt <- sink_open(NULL)
# log.txt <- sink_open(here_output_real("log"))
tic()
out <- cl_lapply(ids, function(id) {
  patter_workflow(id = id, 
                  moorings = moorings, detections = detections, cthresholds = cthresholds,
                  model_move_type = "real", 
                  n_particle = 1e4L,
                  trial = TRUE,
                  here_output = here_output_real)
})
toc()
sink_close(log.txt)


#### End of code. 
###########################
###########################