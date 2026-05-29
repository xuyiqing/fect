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
