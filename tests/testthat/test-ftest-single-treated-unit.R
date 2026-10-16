test_that("the pre-trend F test with one treated unit returns NA instead of stopping the fit", {
  set.seed(1)
  n <- 20; tt <- 12
  d <- expand.grid(id = 1:n, time = 1:tt)
  d$D <- as.integer(d$id == 1 & d$time == tt)        # one treated unit, treated in the last period
  d$Y <- rnorm(nrow(d)) + d$id / 10 + d$time / 5 + d$D
  fit <- suppressMessages(fect(Y ~ D, data = d, index = c("id", "time"), method = "fe",
                               se = TRUE, nboots = 20, parallel = FALSE))
  expect_true(is.finite(fit$att.avg))
  expect_true(is.na(fit$test.out$f.stat))
  expect_true(is.na(fit$test.out$f.p))
})

## A one-period window for the pre-trend F test, with one treated unit. This is
## the case in which the F test's small-sample factor is 0/0 (NaN); the default
## window above has several periods and is not enough to expose the old error.
## The fit must still return its estimates, with the F statistic and p-value
## reported as NA and a message that says why.
test_that("a one-period pre-trend window with one treated unit gives NA, not an error", {
  set.seed(1)
  n <- 20; tt <- 12
  d <- expand.grid(id = 1:n, time = 1:tt)
  d$D <- as.integer(d$id == 1 & d$time == tt)
  d$Y <- rnorm(nrow(d)) + d$id / 10 + d$time / 5 + d$D
  for (pp in list(0, -1, c(0, 0))) {
    expect_message(
      fit <- fect(Y ~ D, data = d, index = c("id", "time"), method = "fe",
                  se = TRUE, nboots = 20, parallel = FALSE, pre.periods = pp),
      "F statistic")
    expect_true(is.finite(fit$att.avg))
    expect_true(is.na(fit$test.out$f.stat))
    expect_true(is.na(fit$test.out$f.p))
  }
})
