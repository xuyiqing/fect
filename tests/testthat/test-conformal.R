## Tests for the estimator-agnostic conformal core (R/conformal.R) AND the full
## leave-one-control-out calibration (conformal_calibrate, via impute_Y0).
## Coverage is validated here through the calibration itself; the fect() output-slot
## wiring (eff.calendar / est.* from vartype = "conformal") is still in progress and
## its print/plot tests are added with that wiring.

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

## ---- grid inversion (dispersion scores) -------------------------------------

test_that("grid denominator is outcome-scaled per score", {
  e.pre <- rnorm(30)
  expect_equal(fect:::.conformal_denom(e.pre, "studentized"), stats::sd(e.pre))
  expect_equal(fect:::.conformal_denom(e.pre, "ratio"), sqrt(mean(e.pre^2)))
  expect_equal(fect:::.conformal_denom(e.pre, "rmse"), 1)
})

test_that("grid interval brackets a constant effect", {
  set.seed(2)
  s.co     <- abs(rnorm(40, mean = 1, sd = 0.3))   # control scores
  gap.post <- 3 + rnorm(8, sd = 0.3)               # true constant effect = 3
  gap.pre  <- rnorm(15, sd = 1)
  ci <- fect:::.conformal_ci_grid(gap.post, gap.pre, s.co, "rmse", alpha = 0.10)
  expect_true(all(is.finite(ci)))
  expect_lt(ci[1], 3); expect_gt(ci[2], 3)
})

test_that("grid returns an explicit empty set when every tau is rejected", {
  set.seed(3)
  gap.post <- c(-10, 10, -10, 10, -10, 10)         # huge dispersion: no tau fits
  gap.pre  <- rnorm(15, sd = 1)
  s.co     <- abs(rnorm(40, sd = 0.5))             # tight control scores
  ci <- fect:::.conformal_ci_grid(gap.post, gap.pre, s.co, "rmse", alpha = 0.10)
  expect_true(all(is.na(ci)))
  expect_identical(attr(ci, "empty"), "rejected_all_tau")
})

test_that("grid is unbounded when alpha is below the resolution floor", {
  set.seed(4)
  gap.post <- rnorm(5, sd = 1); gap.pre <- rnorm(15, sd = 1)
  s.co <- abs(rnorm(20))
  ci <- fect:::.conformal_ci_grid(gap.post, gap.pre, s.co, "rmse", alpha = 1 / 25) # < 1/21
  expect_true(any(is.infinite(ci)))
})

test_that("grid expands when the initial span is too tight", {
  set.seed(7)
  gap.post <- 5 + rnorm(6, sd = 0.3); gap.pre <- rnorm(15, sd = 1)
  s.co <- abs(rnorm(40, mean = 1, sd = 0.3))
  ## start far inside the true center (~5); adaptive doubling must still recover it
  ci <- fect:::.conformal_ci_grid(gap.post, gap.pre, s.co, "rmse", alpha = 0.10, span = 0.5)
  expect_true(all(is.finite(ci)))
  expect_lt(ci[1], 5); expect_gt(ci[2], 5)
})

## ---- full leave-one-control-out calibration (coverage) ----------------------
## A block DGP with a known factor structure and zero true effect. Each rep fits
## gsynth, calibrates against the held-out donor scores, and records coverage.
## These are regression guards (no unbounded intervals; empties always flagged;
## coverage near nominal) -- precise calibration is the coverage study's job.

.conf_dgp <- function(seed, N = 25, T = 20, T0 = 15, r = 2) {
  set.seed(seed)
  Fm <- matrix(rnorm(T * r), T, r); L <- matrix(rnorm(N * r), N, r)
  D  <- matrix(0, T, N); D[(T0 + 1):T, 1] <- 1
  list(mu = Fm %*% t(L), D = D, N = N, T = T)
}

test_that("meanabs calibration covers near nominal and never empties/unbounds", {
  skip_on_cran()
  g <- .conf_dgp(20260608); reps <- 50; alpha <- 0.10
  cov <- logical(reps); bad <- 0L
  for (b in seq_len(reps)) {
    Y <- g$mu + matrix(rnorm(g$N * g$T), g$T, g$N)
    dat <- data.frame(id = rep(1:g$N, each = g$T), time = rep(1:g$T, g$N),
                      Y = as.vector(Y), D = as.vector(g$D))
    fit <- suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
              method = "gsynth", force = 3, CV = FALSE, r = 2, se = FALSE))
    cc <- conformal_calibrate(Y = Y, D = g$D, I = fit$I, II = fit$II, T.on = fit$T.on,
              r.cv = fit$r.cv, eff = fit$eff, method = "gsynth",
              score = "meanabs", alpha = alpha)
    if (cc$status != "ok") bad <- bad + 1L
    cov[b] <- (0 >= cc$ci[1]) && (0 <= cc$ci[2])
  }
  expect_equal(bad, 0L)              # meanabs: no empty, no unbounded
  expect_gt(mean(cov), 0.82)         # >= nominal, allowing MC noise over 50 reps
  expect_lte(mean(cov), 1.0)
})

test_that("studentized calibration is calibrated and flags every empty", {
  skip_on_cran()
  ## 80 reps so the coverage bound is robust to Monte-Carlo noise (SE ~ 0.033 at
  ## nominal 0.90); the bound is a gross-miscalibration guard, not a precise check.
  g <- .conf_dgp(424242); reps <- 80; alpha <- 0.10
  cov.tot <- logical(reps); na.unflagged <- 0L; empties <- 0L
  for (b in seq_len(reps)) {
    Y <- g$mu + matrix(rnorm(g$N * g$T), g$T, g$N)
    dat <- data.frame(id = rep(1:g$N, each = g$T), time = rep(1:g$T, g$N),
                      Y = as.vector(Y), D = as.vector(g$D))
    fit <- suppressMessages(fect(Y ~ D, data = dat, index = c("id", "time"),
              method = "gsynth", force = 3, CV = FALSE, r = 2, se = FALSE))
    cc <- suppressWarnings(conformal_calibrate(Y = Y, D = g$D, I = fit$I, II = fit$II,
              T.on = fit$T.on, r.cv = fit$r.cv, eff = fit$eff, method = "gsynth",
              score = "studentized", alpha = alpha))
    if (any(is.na(cc$ci))) {
      empties <- empties + 1L
      if (cc$status != "empty") na.unflagged <- na.unflagged + 1L  # never silent
    }
    ## total coverage: an empty set does NOT contain the truth -> FALSE
    cov.tot[b] <- if (any(is.na(cc$ci))) FALSE else (0 >= cc$ci[1]) && (0 <= cc$ci[2])
  }
  expect_equal(na.unflagged, 0L)            # every empty interval is flagged
  expect_lt(empties / reps, 0.30)           # empties are the exception, not the rule
  expect_gt(mean(cov.tot), 0.78)            # total coverage (empties = miss) ~ nominal
})
