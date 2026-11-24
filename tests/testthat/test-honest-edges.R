# Code adapted from the 'HonestDiD' package, version 0.2.8, by Ashesh Rambachan and Jonathan Roth,
# through 'HonestDiDFEct' 1.0.0 by Ziyi Liu (github.com/lzy318/HonestDiDFEct, commit 29c831d).
# Original file: tests/testthat/test-edges.R
# Copyright (c) 2023 HonestDiD authors. MIT licence; the full notice is in inst/COPYRIGHTS.
# Copied into fect 2.5.1 without change to any numerical line; the edits are listed in the fect repository history.
#
# The two edits to R/arp-nuisance.R and the one-post-period path.

test_that("a test statistic with zero variance gives a result, not an error", {
  skip_unless_honest_solvers()
  S <- rbind(c(1, -1, 0), c(-1, 1, 0), c(0, 0, 1))
  X <- matrix(c(1, -1, 0), ncol = 1)
  go <- function(y) {
    fect:::.lp_conditional_test_fn(theta = 0, y_T = y, X_T = X, sigma = S, alpha = 0.05, hybrid_flag = "ARP",
                                            hybrid_list = list(hybrid_kappa = 0.005), rowsForARP = 1:3)
  }
  above <- go(c(0.1, 0.2, -5))
  expect_equal(above$eta, 0.15)
  expect_equal(above$reject, 1)
  below <- go(c(-0.1, -0.2, -5))
  expect_equal(below$eta, -0.15)
  expect_equal(below$reject, 0)
})

test_that("a 1 x 1 matrix at a grid end gives no warning in the nuisance branch", {
  skip_unless_honest_solvers()
  h <- by_hand(path_b, path_V, 3L, 4L, path_l)
  expect_silent(suppressMessages(
    fect:::computeConditionalCS_DeltaSDB(h$betahat, h$sigma, 2L, 4L, l_vec = path_l, M = 0.02,
                                                  hybrid_flag = "LF", gridPoints = 20L)
  ))
})

test_that("one post period works with l_vec = 1", {
  skip_unless_honest_solvers()
  b <- path_b[1:4]
  V <- path_V[1:4, 1:4]
  h <- by_hand(b, V, 3L, 1L, 1)
  r <- quiet(fect:::createSensitivityResults(b, V, 3, 1, method = "C-LF", Mvec = 0.02))
  e <- quiet(fect:::computeConditionalCS_DeltaSD(h$betahat, h$sigma, 2L, 1L, l_vec = matrix(1), M = 0.02,
                                                          hybrid_flag = "LF"))
  expect_equal(r$lb, min(e$grid[e$accept == 1]) + h$shift)
  expect_equal(r$ub, max(e$grid[e$accept == 1]) + h$shift)
})

test_that("the default Mvec is worked out on the differenced path", {
  skip_unless_honest_solvers()
  b <- path_b[2:4]
  V <- path_V[2:4, 2:4]
  r <- quiet(fect:::createSensitivityResults(b, V, 2, 1, method = "C-LF"))
  expect_equal(max(r$M), sqrt(V[1, 1] + V[2, 2] - 2 * V[1, 2]))

  b <- path_b[1:4]
  V <- path_V[1:4, 1:4]
  h <- by_hand(b, V, 3L, 1L, 1)
  r <- quiet(fect:::createSensitivityResults(b, V, 3, 1, method = "C-LF"))
  expect_equal(max(r$M), fect:::DeltaSD_upperBound_Mpre(h$betahat, h$sigma, h$numPrePeriods))
})
