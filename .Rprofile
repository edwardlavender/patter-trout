source("renv/activate.R")

repos <- c(CRAN = "https://cloud.r-project.org")
options(repos = repos)

if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}

Sys.setenv(JULIA_PROJ = here::here("Julia"))
Sys.setenv(JULIA_NUM_THREADS = parallel::detectCores())
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
