# Shared setup for all vignette chapters
# When rendering locally from the source tree, devtools::load_all()
# ensures the LATEST R functions (including new features like cm) and
# all datasets are available, even if the installed package is outdated.
# During R CMD check the package is freshly installed, so library() suffices.
if (file.exists("../DESCRIPTION") &&
    requireNamespace("devtools", quietly = TRUE)) {
  tryCatch(
    suppressMessages(devtools::load_all("..", quiet = TRUE)),
    error = function(e) {
      message("devtools::load_all() failed: ", conditionMessage(e))
      message("Falling back to library(fect)")
      suppressMessages(library(fect))
    }
  )
} else {
  suppressMessages(library(fect))
}
