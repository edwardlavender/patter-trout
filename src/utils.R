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

#' sink_*() wrappers

sink_open <- function(log.folder = NULL) {
  log.txt <- NULL
  if (!is.null(log.folder)) {
    # Define file name
    log.txt <- paste0("log-", as.numeric(Sys.time()), ".txt")
    # Define full file path & validate 
    log.txt <- file.path(log.folder, log.txt)
    stopifnot(!file.exists(log.txt))
    # Define connection 
    log.txt <- file(log.txt, open = "wt")
    # Open sink
    sink(log.txt, append = TRUE)
    sink(log.txt, type = "message", append = TRUE)
    # Print start time
    print(Sys.time())
  }
  invisible(log.txt)
}

sink_close <- function(log.txt = NULL) {
  if (!is.null(log.txt)) {
    print(Sys.time())
    sink()
    sink(type = "message")
  }
  invisible(NULL)
}