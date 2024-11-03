###########################
###########################
#### patter-setup-map.R

#### Aims
# 1) Sets up the map (used in both simulation & real-world analyses)

#### Prerequisites
# 1) Use ChamplainRegionsGrouped shapefile from M. Futia. 


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())
# try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
Sys.setenv("JULIA_SESSION" = FALSE)
library(tictoc)
dv::src()

#### Load data 
champlain  <- terra::vect(here_data("ChamplainRegionsGrouped/ChamplainRegionsGrouped.shp"))


###########################
###########################
#### Study area (~2 s)

#### Build map 
# Use a coarse map for speed sampling initial locations 
tic()
epsg_utm         <- "EPSG:3175"
champlain$land   <- as.numeric(1)
regions <- as_SpatRaster(champlain, .simplify = 0.001, .utm = epsg_utm,
                         .field = "region", .res = 200, .plot = TRUE)
maps    <- as_SpatRaster(champlain, .simplify = 0.001, .utm = epsg_utm,
                         .field = "land", .res = 200, .plot = TRUE)
map     <- maps$SpatRaster
map_len <- terra::ymax(map) - terra::ymin(map)
toc()

#### Examine map properties
# Visualise map simplification
terra::plot(map, col = "blue")
terra::lines(champlain |> terra::project(epsg_utm))
# Zoom-in to check resolution
map_zoom <- terra::crop(map, 
                        cbind(1787707, 977397.3) |>
                          terra::vect() |>
                          terra::buffer(width = 10000) |>
                          terra::ext())
terra::plot(map_zoom)
terra::lines(champlain |> terra::project(epsg_utm))
# Check ncell & compare to dat_gebco() for reference
terra::ncell(map)                  # 107625
terra::ncell(patter::dat_gebco())  # 50160
# Check map size
terra::ncell(map) * 8 / 1e6        # 0.81 MB

#### Map dimensions
# max dimension
terra::ymax(map) - terra::ymin(map) # 175000 m
terra::xmax(map) - terra::xmin(map) # 24600 m
# bbox
map_bbox <- patter:::map_bbox(map)
stopifnot(nrow(map_bbox) == 4L)
terra::plot(map)
points(map_bbox, pch = 21, bg = "red", cex = 5)

#### Write maps
terra::writeRaster(map, here_input("map.tif"), overwrite = TRUE)
terra::writeRaster(regions$SpatRaster, here_input("regions.tif"), overwrite = TRUE)
qs::qsave(map_bbox, here_input("map-bbox.qs"))
file.size(here_input("map.tif")) / 1e6 # MB


#### End of code. 
###########################
###########################