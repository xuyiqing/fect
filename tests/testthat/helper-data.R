# Helpers to ensure datasets are available
if (!exists("simdata", inherits = FALSE) || !exists("simgsynth", inherits = FALSE)) {
  f <- file.path("data", "fect.RData")
  if (file.exists(f)) {
    try(load(f, envir = .GlobalEnv), silent = TRUE)
  }
}
