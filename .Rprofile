source("renv/activate.R")

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}

Sys.setenv(JULIA_PROJ = here::here("Julia"))
Sys.setenv(JULIA_NUM_THREADS = parallel::detectCores())
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
