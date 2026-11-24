skip_unless_honest_solvers <- function(extra = character(0)) {
  for (p in c("tibble", "purrr", "TruncatedNormal", "matrixStats", "Rglpk", "pracma", "lpSolveAPI", "Matrix", extra))
    testthat::skip_if_not_installed(p)
}

# Compare two numeric vectors position by position. Where both are finite the absolute difference must be
# at most `tol`; elsewhere the two must agree in kind (NA, NaN, Inf, -Inf). Returns the largest finite
# absolute difference invisibly (0 if none).
expect_close <- function(a, b, tol = 1e-8, info = NULL) {
  kind <- function(x) ifelse(is.na(x) & !is.nan(x), "NA",
                      ifelse(is.nan(x), "NaN",
                      ifelse(is.infinite(x), ifelse(x > 0, "Inf", "-Inf"), "finite")))
  if (length(a) != length(b)) {
    testthat::fail(paste0(if (is.null(info)) "" else paste0(info, ": "),
                          "lengths differ: ", length(a), " vs ", length(b)))
    return(invisible(0))
  }
  a <- as.numeric(a)
  b <- as.numeric(b)
  ka <- kind(a)
  kb <- kind(b)
  both <- ka == "finite" & kb == "finite"
  d <- abs(a - b)
  worst <- if (any(both)) max(d[both]) else 0
  bad <- which((ka != kb) | (both & d > tol))
  if (length(bad) > 0L) {
    i <- bad[1L]
    testthat::fail(paste0(if (is.null(info)) "" else paste0(info, ": "),
                          "first mismatch at position ", i, " (", format(a[i]), " vs ", format(b[i]),
                          "); largest absolute difference ", format(worst), " (tolerance ", format(tol), ")"))
  } else {
    testthat::succeed()
  }
  invisible(worst)
}
