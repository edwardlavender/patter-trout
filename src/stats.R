# Add truncated gamma curves to ggplot2
stat_dtruncgamma <- function(mobility, shape, scale, col = "black", ...) {
  ggplot2::stat_function(fun = function(x) {
    truncdist::dtrunc(x, "gamma", a = 0, b = mobility, shape = shape, scale = scale)
  }, colour = col, ...)
}

# Get parameters of shrunk/widened gamma distribution
gamma_rescale <- function(shape, scale, fact = 1.2) {
  mode      <- (shape - 1) * scale
  new_shape <- shape / fact
  new_scale <- mode / (new_shape - 1)
  c(new_shape, new_scale)
}
