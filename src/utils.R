#' here::here() wrappers

here_fig <- function(...) {
  here::here("fig", ...)
}

here_fig_sim <- function(...) {
  here_fig("sim", ...)
}

here_fig_real <- function(...) {
  here_fig("real", ...)
}

here_data <- function(...) {
  here::here("data", ...)
}

here_input <- function(...) {
  here::here("data", "patter", "input", ...)
}

here_input_sim <- function(...) {
  here_input("sim", ...)
}

here_input_real <- function(...) {
  here_input("real", ...)
}

here_output <- function(...) {
  here::here("data", "patter", "output", ...)
}

here_output_sim <- function(...) {
  here_output("sim", ...)
}

here_output_real <- function(...) {
  here_output("real", ...)
}

#' os_*() wrappers

os_linux <- function() {
  grepl("linux", tolower(Sys.info()["sysname"]))
}
