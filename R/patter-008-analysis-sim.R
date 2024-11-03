###########################
###########################
#### patter-analysis-sim.R

#### Aims
# 1) Quick residency analysis for simulated data

#### Prerequisites
# 1) patter-algorithms.R

#
#
# TO DO
# Update script in line with revised directory structure
# Rename sim_id -> individual_id
# Handle inputs from different models etc.
#
#


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
library(glue)
library(patter)
library(tictoc)
dv::src()

#### Load data
map_region <- terra::rast(here_input("regions.tif"))
paths      <- qs::qread(here_input("paths.qs"))

#### Set options
op <- options(terra.pal = rev(terrain.colors(256L)))


###########################
###########################
#### Algorithm properties

#### Define runs
runs <- c("forward", "backward", "smoothing")

#### Read convergence properties
convergence <- 
  here_output("convergence") |>
  list.files(full.names = TRUE) |> 
  lapply(qs::qread) |>
  rbindlist() |> 
  mutate(direction = factor(direction, levels = runs)) |> 
  arrange(sim_id, direction) |> 
  as.data.table()

#### Examine convergence failures
convergence[direction == "forward" & convergence == 0L, ]
convergence[direction == "backward" & convergence == 0L, ]

#### Examine computation time
if (requireNamespace("utils.add", quietly = TRUE)) {
  # On siam-linux20:
  # * ~ 1.5 hours for each filter (x 2)
  # * ~ 2.0 hours for smoothing
  utils.add::basic_stats(convergence$time[convergence$direction != "smoothing"])
  utils.add::basic_stats(convergence$time[convergence$direction == "smoothing"])
}


###########################
###########################
#### Example maps

#### Define example individual (path & states)
id     <- 1L
path   <- paths[sim_id == id, ]
states <- qs::qread(here_output("particles", glue("smo-{id}.qs")))$states

#### Visualise regions, path & states
# The reconstructed map should look reasonable
pp <- par(mfrow = c(1, 3))
tic()
terra::plot(map_region)
map_pou(.map = map_region, .coord = path)
map_pou(.map = map_region, .coord = states)
toc()
par(pp)


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
  rbindlist() |> 
  arrange(sim_id, region, estimate) |> 
  as.data.table()
# Examine estimates
head(residency, 10)
residency[is.na(region), ]

#### Visualise residency ~ individual, coloured by truth/algorithm
# NB: residency includes an NA category 
# This is probably due to the simplification of the regions polygon 
png(here_fig("residency.png"), 
    height = 5, width = 10, units = "in", res = 600)
residency |>
  ggplot(aes(x = region, y = perc, fill = estimate)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ sim_id) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Region", y = "Percentage", fill = "Estimate")
dev.off()


#### End of code.
###########################
###########################