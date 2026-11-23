## Minimal CRAN smoke tests for fect.
##
## This is the ONLY test file that runs on CRAN: tests/testthat.R filters the
## suite to this file when NOT_CRAN is unset (see that driver). Keep it fast and
## dependency-light -- a single core fit and the print method. The full
## regression suite (every other test-*.R) runs locally and in CI, where
## NOT_CRAN=true.

library(testthat)
data(simdata, package = "fect")

test_that("core fect() fit runs and returns the expected structure", {
  out <- fect::fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                    method = "fe", force = "two-way", se = FALSE, parallel = FALSE)
  expect_s3_class(out, "fect")
  expect_true(is.matrix(out$eff))
  expect_true(is.numeric(out$att.avg))
  expect_true(length(out$time) >= 1L)
})

test_that("print.fect() runs without error", {
  out <- fect::fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                    method = "fe", force = "two-way", se = FALSE, parallel = FALSE)
  printed <- capture.output(print(out))
  expect_type(printed, "character")
})

test_that("fect_sens engine reproduces the stored bounds", {
  fx <- readRDS(testthat::test_path("fixtures", "honest-sens-ref.rds"))
  cases <- Filter(function(z) isTRUE(z$cran), fx$cases)
  expect_true(length(cases) >= 1L)
  for (p in unique(unlist(lapply(cases, `[[`, "needs")))) testthat::skip_if_not_installed(p)
  worst <- 0
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    res <- as.data.frame(suppressWarnings(suppressMessages(
      do.call(get(cs$fn, envir = asNamespace("fect")), cs$args))))
    expect_identical(names(res), names(cs$expected), info = nm)
    expect_identical(nrow(res), nrow(cs$expected), info = nm)
    for (v in c("lb", "ub")) {
      d <- abs(res[[v]] - cs$expected[[v]])
      fin <- is.finite(res[[v]]) & is.finite(cs$expected[[v]])
      expect_identical(is.finite(res[[v]]), is.finite(cs$expected[[v]]), info = paste(nm, v))
      expect_true(all(d[fin] <= fx$tolerance), info = paste(nm, v))
      if (any(fin)) worst <- max(worst, d[fin])
    }
  }
  cat("\nhonest-fixture check:", length(cases), "case(s) compared, largest difference", format(worst), "\n")
})
