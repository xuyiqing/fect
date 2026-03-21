## Ensure datasets are available in tests across installed/build contexts
if (!exists("simdata", inherits = TRUE)) {
  suppressWarnings(try(utils::data("simdata", package = "fect"), silent = TRUE))
}
if (!exists("simdata", inherits = TRUE)) {
  f <- system.file("data", "fect.RData", package = "fect")
  if (nzchar(f) && file.exists(f)) {
    try(load(f, envir = environment()), silent = TRUE)
  }
}
if (!exists("simdata", inherits = TRUE)) {
  f <- file.path("data", "fect.RData")
  if (file.exists(f)) {
    try(load(f, envir = environment()), silent = TRUE)
  }
}
if (!exists("simgsynth", inherits = TRUE)) {
  suppressWarnings(try(utils::data("simgsynth", package = "fect"), silent = TRUE))
}
