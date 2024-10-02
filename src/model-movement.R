states_init.StateXYD <- function(.state, .coords) {
  z <- map_value <- angle <- NULL
  # Add angle (rad)
  .coords[, angle := runif(.N, 0, 2 * pi)]
  .coords
}

move_xyd <- function(length, theta) {
  length <- as.numeric(length)
  theta  <- as.numeric(theta)
  julia_assign("model_move_length", length)
  julia_assign("model_move_theta", theta)
  patter:::julia_check_exists("env")
  'ModelMoveXYD(env, model_move_length, Normal(0.0, model_move_theta));'
}
