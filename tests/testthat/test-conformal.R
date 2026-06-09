## Tests for the estimator-agnostic conformal core (R/conformal.R).
## The full vartype = "conformal" pipeline (leave-one-control-out calibration via
## impute_Y0) runs end to end but its output-slot integration and coverage
## validation are in progress; those tests are added with that wiring.

test_that("conformal scores are finite and ordered as expected", {
  set.seed(1)
  e.post <- rnorm(12); e.pre <- rnorm(24)
  for (s in c("meanabs", "rmse", "ratio", "studentized")) {
    expect_true(is.finite(fect:::.conformal_score(e.post, e.pre, s)))
  }
  ## a larger post deviation gives a larger score (meanabs)
  expect_gt(fect:::.conformal_score(e.post + 5, e.pre, "meanabs"),
            fect:::.conformal_score(e.post,     e.pre, "meanabs"))
  ## ratio with zero pre-variation is NA, not Inf-by-accident
  expect_true(is.na(fect:::.conformal_score(e.post, rep(0, 24), "ratio")))
})

test_that("meanabs closed-form interval covers near nominal", {
  set.seed(20260608)
  reps <- 3000; alpha <- 0.10; Nco <- 39; att <- 3
  cov <- replicate(reps, {
    g.co <- rnorm(Nco)                 # no-effect donor mean gaps
    g.tr <- att + rnorm(1)             # treated mean gap = att + exchangeable noise
    ci <- fect:::.conformal_ci_meanabs(g.tr, g.co, alpha = alpha)
    (att >= ci[1]) && (att <= ci[2])
  })
  expect_gt(mean(cov), 0.86)
  expect_lt(mean(cov), 0.94)           # ~ 1 - alpha
})

test_that("conformal p-value lies in (0,1] and respects weights", {
  s.co <- rnorm(20)
  p <- fect:::.conformal_pval(0, s.co)
  expect_gt(p, 0); expect_lte(p, 1)
  ## a very large treated score is the most extreme -> smallest p = 1/(n+1)
  expect_equal(fect:::.conformal_pval(100, s.co), 1 / (length(s.co) + 1))
  ## zero weight on a donor drops it from the count
  w <- rep(1, 20); w[which.max(s.co)] <- 0
  expect_true(is.finite(fect:::.conformal_pval(0, s.co, w.co = w)))
})

test_that("resolution guard: alpha below 1/(Nco+1) yields an unbounded interval", {
  ci <- fect:::.conformal_ci_meanabs(0, rnorm(20), alpha = 1 / 25)  # 1/25 < 1/21
  expect_true(all(is.infinite(ci)))
})

test_that("effective sample size is sensible", {
  expect_equal(fect:::.conformal_neff(rep(1, 10)), 10)
  expect_lt(fect:::.conformal_neff(c(rep(1, 9), 100)), 10)  # one big weight -> n_eff < n
})
