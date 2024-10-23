move_xyd <- function(mobility = "750.0",
                     dbn_length = "truncated(Gamma(1.0, 750.0), upper = 750.0)",
                     dbn_angle_delta = "Normal(0, 0.5)") {
  patter:::julia_check_exists("env")
  glue('ModelMoveXYD(env, {mobility}, {dbn_length}, {dbn_angle_delta});')
}