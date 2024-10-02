# (TEMPORARY) as_SpatRaster() function

# TO DO
# * Rename to assemble_maps()
# * If SpatVector supplied, implement as shown
# * If SpatRaster supplied, create SpatVector (?)
# * Add to patter

as_SpatRaster <- function(.x, .simplify = NULL, .utm = NULL,
                          .res, .field, .touches = TRUE, 
                          .plot = FALSE, ...) {
  
  # Translate sf objects to SpatVectors
  if (inherits(.x, "sf")) {
    .x <- terra::vect(.x)
  }
  
  # (optional) Simplify SpatVector 
  # > This improves speed
  if (!is.null(.simplify)) {
    .x <- terra::simplifyGeom(.x, tolerance = .simplify)
  }
  
  # Translate to UTM
  if (!is.null(.utm)) {
    .x <- terra::project(.x, .utm)
  }
  
  # Define map_value
  .x$map_value <- .x[[.field]]
  
  # Define blank raster for rasterization 
  r <- terra::rast(terra::ext(.x), 
                   crs = terra::crs(.x),
                   res = .res)
  
  # Rasterise the SpatVector 
  map <- terra::rasterize(.x, r, field = "map_value", touches = .touches, ...)
  
  # (optional) Plot 
  if (.plot) {
    terra::plot(map)
    terra::lines(.x)
  }
  
  # Return list
  list(SpatVector = .x, SpatRaster = map)
  
}