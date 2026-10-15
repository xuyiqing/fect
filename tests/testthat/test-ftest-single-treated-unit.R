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
