# Code adapted from the 'HonestDiD' package, version 0.2.8, by Ashesh Rambachan and Jonathan Roth,
# through 'HonestDiDFEct' 1.0.0 by Ziyi Liu (github.com/lzy318/HonestDiDFEct, commit 29c831d).
# Original file: tests/testthat/test-honestdid.R
# Copyright (c) 2023 HonestDiD authors. MIT licence; the full notice is in inst/COPYRIGHTS.
# Copied into fect 2.5.1 without change to any numerical line; the edits are listed in the fect repository history.
#
# The bounds equal HonestDiD on the differenced, recentred path plus the shift. These tests catch a merge
# from upstream that changes the engine in a way that does not carry over.
base_version <- "0.2.8"

compare <- function(fn, ...) {
  h <- by_hand(path_b, path_V, 3L, 4L, path_l)
  a <- quiet(get(fn, envir = asNamespace("fect"))(path_b, path_V, 3, 4, l_vec = path_l, ...))
  b <- quiet(getExportedValue("HonestDiD", fn)(h$betahat, h$sigma, h$numPrePeriods, 4L, l_vec = path_l, ...))
  expect_equal(a$lb, b$lb + h$shift, tolerance = 1e-10)
  expect_equal(a$ub, b$ub + h$shift, tolerance = 1e-10)
  expect_identical(a$method, b$method)
  expect_identical(a$Delta, b$Delta)
}

skip_unless_base <- function() {
  skip_on_cran()
  skip_if_not_installed("HonestDiD")
  skip_if(as.character(utils::packageVersion("HonestDiD")) != base_version)
}

test_that("SD and its variants match HonestDiD", {
  skip_unless_honest_solvers(c("CVXR", "ECOSolveR"))
  skip_unless_base()
  compare("createSensitivityResults", method = "FLCI", Mvec = 0.02)
  compare("createSensitivityResults", method = "Conditional", Mvec = 0.02)
  compare("createSensitivityResults", method = "C-LF", Mvec = 0.02, biasDirection = "positive")
  compare("createSensitivityResults", method = "C-LF", Mvec = 0.02, monotonicityDirection = "increasing")
})

test_that("RM and SDRM match HonestDiD", {
  skip_unless_honest_solvers()
  skip_unless_base()
  compare("createSensitivityResults_relativeMagnitudes", Mbarvec = 0.5, gridPoints = 50)
  compare("createSensitivityResults_relativeMagnitudes", Mbarvec = 0.5, gridPoints = 50,
          bound = "deviation from linear trend")
})
