## did before 2.5.1 reported a CSDID cell it could not estimate on an
## unbalanced panel as 0 and averaged it in; did_wrapper() now needs 2.5.1.
## Here one cohort (adopting in period 3) has no rows after period 3, so its
## cells for periods 4 and 5 cannot be estimated. The true effect is 1.

test_that("did_wrapper CSDID leaves unestimable cells out of the average on an unbalanced panel", {
  skip_on_cran()
  skip_if_not_installed("did", "2.5.1")
  set.seed(1)
  d <- expand.grid(id = 1:200, period = 1:5)
  g <- ifelse(d$id <= 80, 3, Inf)
  d$D <- as.integer(d$period >= g)
  d$y <- stats::rnorm(nrow(d)) + d$D
  d <- d[!(g < Inf & d$period > 3), ]
  d <- d[!(d$id %% 7 == 0 & d$period == 2), ]
  out <- suppressWarnings(suppressMessages(
    did_wrapper(d, Y = "y", D = "D", index = c("id", "period"), method = "cs_notyet", se = "default")))
  att <- as.data.frame(out$est.avg)$ATT.avg[1]
  expect_true(is.finite(att))
  expect_gt(att, 0.8)   # 0.331 with did 2.3.0, which averaged two placeholder zeros in
})
