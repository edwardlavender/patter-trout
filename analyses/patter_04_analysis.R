###########################
###########################
#### patter-analyses.R

#### Aims
# 1) Quick residency analysis

#### Prerequisites
# 1) patter-algorithms.R


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())
# try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tictoc)
dv::src()

#### Load data
map_region <- terra::rast(here_input("regions.tif"))
paths      <- qs::qread(here_input("paths.qs"))


###########################
###########################
#### Algorithm properties

# Define runs
runs <- c("forward", "backward", "smoothing")

# Read convergence properties
convergence <- 
  here_output("convergence") |>
  list.files(full.names = TRUE) |> 
  lapply(qs::qread) |>
  rbindlist() |> 
  mutate(direction = factor(direction, levels = runs)) |> 
  arrange(sim_id, direction) |> 
  as.data.table()

# Examine convergence failures
convergence[direction == "forward" & convergence == 0L, ]
convergence[direction == "backward" & convergence == 0L, ]

# Examine computation time
utils.add::basic_stats(convergence$time[convergence$direction != "smoothing"])
utils.add::basic_stats(convergence$time[convergence$direction == "smoothing"])


###########################
###########################
#### Quick residency analysis

#### Compute quick residency statistics
# Compute residency 
tic()
cl_lapply(1:100L, function(id) {
  qresidency(id = id, map = map_region, paths = paths)
})
toc()
# Read residency estimates
residency <- 
  here_output("residency", "qresidency") |>
  list.files(full.names = TRUE) |>
  lapply(qs::qread) |> 
  rbindlist()

#### Visualise residency ~ individual, coloured by truth/algorithm
# 
# > TO DO


#### End of code.
###########################
###########################