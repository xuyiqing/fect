## A failed singular value decomposition stops with a clear message instead of
## an Armadillo index error (panel_factor: IFE; panel_FE: matrix completion).
## The divide-and-conquer failure itself was seen on a 183 x 327 matrix from a
## very unbalanced panel and is not reproduced here; these tests cover the
## guard and check that a regular matrix still decomposes as before.

test_that("panel_factor and panel_FE stop with a clear message on non-finite input", {
  E <- matrix(stats::rnorm(60), 10, 6)
  E[3, 2] <- NaN
  expect_error(fect:::panel_factor(E, 2L), "non-finite")
  expect_error(fect:::panel_FE(E, 0.1, 0L), "non-finite")
  E[3, 2] <- Inf
  expect_error(fect:::panel_factor(E, 2L), "non-finite")
})

test_that("panel_factor and panel_FE still decompose a regular matrix", {
  set.seed(1)
  f <- matrix(stats::rnorm(20), 10, 2)
  l <- matrix(stats::rnorm(12), 6, 2)
  E <- f %*% t(l) + matrix(stats::rnorm(60, sd = 0.01), 10, 6)
  pf <- fect:::panel_factor(E, 2L)
  expect_equal(dim(pf$factor), c(10L, 2L))
  expect_equal(dim(pf$lambda), c(6L, 2L))
  expect_true(max(abs(pf$FE - E)) < 0.05)
  # wide matrix (T < N) takes the other branch
  pf_t <- fect:::panel_factor(t(E), 2L)
  expect_true(max(abs(pf_t$FE - t(E))) < 0.05)
  fe <- fect:::panel_FE(E, 0, 0L)  # no shrinkage: the full reconstruction
  expect_equal(fe, E, tolerance = 1e-8)
})
