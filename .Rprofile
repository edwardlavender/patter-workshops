source("renv/activate.R")

if (!requireNamespace("here", quietly = TRUE)) {
  renv::install("here", prompt = FALSE)
}

Sys.setenv(JULIA_PROJ = here::here("Julia"))

Sys.setenv(JULIA_NUM_THREADS = parallel::detectCores())

Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
