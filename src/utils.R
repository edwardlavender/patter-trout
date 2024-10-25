#' here::here() wrappers

here_fig <- function(...) {
  here::here("fig", ...)
}

here_input <- function(...) {
  here::here("data", "patter", "input", ...)
}

here_output <- function(...) {
  here::here("data", "patter", "output", ...)
}

here_output_sim <- function(...) {
  here::here("data", "patter", "output", "simulation", ...)
}

#' os_*() wrappers

os_linux <- function() {
  grepl("linux", tolower(Sys.info()["sysname"]))
}
