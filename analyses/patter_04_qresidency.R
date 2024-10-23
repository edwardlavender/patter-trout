###########################
###########################
#### patter-qresidency.R

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
paths       <- qs::qread(here_input("paths.qs"))


###########################
###########################
#### Quick residency analysis

#### TO DO
# * Iterate over indidividuals
# * Compute true & estimated residency 

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