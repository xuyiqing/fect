# Regression test: starting values for the covariate coefficients when feols()
# drops a collinear covariate (found 2026-10-27, gsynth-revisit Carter 2025).
#
# align_beta0() read names(beta_raw) after as.numeric(beta_raw), which drops
# the names, so it always placed coefficients by position. When feols() drops
# x.k, every later coefficient moved one slot left. The main fit is protected
# (fect.default drops such covariates first), but the parametric bootstrap's
# refits are not: in a resampled set of controls two covariates can coincide.
# The EM solver then started from, e.g., a per-capita coefficient on total GDP
# and stopped at max.iteration far from the solution. In Carter dem-filter one
# such draw (Hungary as the pseudo-treated unit) had prediction errors up to
# 2,557, inflated the AR error covariance, and gave SE 35.8 instead of 0.08.

test_that("align_beta0 places coefficients by name when one is missing", {
  b <- c(x.1 = 1, x.2 = 2, x.4 = 4)
  expect_equal(as.numeric(fect:::align_beta0(b, 4)), c(1, 2, 0, 4))
  b10 <- c(x.1 = 1, x.3 = 3, x.10 = 10)
  expect_equal(as.numeric(fect:::align_beta0(b10, 10)), c(1, 0, 3, rep(0, 6), 10))
  ## unnamed input (lm() path keeps NA in place): positional, NA -> 0
  expect_equal(as.numeric(fect:::align_beta0(c(1, NA, 3), 3)), c(1, 0, 3))
})

## A control-only panel shaped like the Carter (2025) pseudo sample (what
## fect_nevertreated() fits, and what each parametric-bootstrap draw refits):
## x1 is a rare event and x2 equals it on every cell (a coup attempt that also
## succeeded), so feols() drops x.2; x3 is on the scale of total GDP, x4 of GDP
## per capita, x5 an electrification rate stuck at 100 for most units; 15% of
## cells are missing, so the EM solver inter_fe_ub() is used.
.ab0_controls <- function() {
  set.seed(1)
  N <- 40; TT <- 20
  id <- rep(seq_len(N), each = TT); tt <- rep(seq_len(TT), N)
  x1 <- as.numeric(stats::runif(N * TT) < 0.01); x2 <- x1
  x3 <- stats::runif(N, 1e9, 1e12)[id] * (1 + 0.04 * tt) * exp(stats::rnorm(N * TT, 0, 0.02))
  x4 <- stats::runif(N, 500, 30000)[id] * (1 + 0.03 * tt) + stats::rnorm(N * TT, 0, 50)
  x5 <- pmin(100, stats::runif(N, 60, 130)[id] + 0.5 * tt)
  x6 <- stats::rnorm(N * TT)
  y <- stats::rnorm(N)[id] + 0.2 * x1 + 2e-13 * x3 + 1e-5 * x4 + 0.01 * x5 + 0.1 * x6 +
    stats::rnorm(N * TT, 0, 0.3)
  obs <- stats::runif(N * TT) > 0.15
  list(data = cbind(y, id, tt, x1, x2, x3, x4, x5, x6), obs = obs, N = N, TT = TT)
}
.ab0_ref <- function(s) {
  fixest::feols(y ~ x1 + x3 + x4 + x5 + x6 | id, data = as.data.frame(s$data[s$obs, ]),
                notes = FALSE)
}

test_that("initialFit's starting coefficients follow feols's names when it drops a covariate", {
  skip_if_not_installed("fixest")
  s <- .ab0_controls()
  ini <- suppressMessages(fect:::initialFit(data = s$data, force = 1, oci = which(s$obs)))
  b <- coef(.ab0_ref(s))
  ## before the fix: c(b1, b3, b4, b5, b6, 0), so x3 (about 1e11) got x4's coefficient
  expect_equal(as.numeric(ini$beta0), unname(c(b["x1"], 0, b["x3"], b["x4"], b["x5"], b["x6"])),
               tolerance = 1e-8)
})

test_that("the EM solver started from initialFit reaches the fixed-effects solution", {
  ## before the fix: niter 5001 and a coefficient of about 24,840 on x5 (truth 0.01)
  skip_if_not_installed("fixest")
  s <- .ab0_controls()
  ini <- suppressMessages(fect:::initialFit(data = s$data, force = 1, oci = which(s$obs)))
  Y <- matrix(s$data[, "y"], s$TT, s$N); I <- matrix(as.numeric(s$obs), s$TT, s$N)
  X <- array(s$data[, 4:9], dim = c(s$TT, s$N, 6))
  Y[I == 0] <- 0
  b0 <- ini$beta0; b0[is.na(b0)] <- 0
  out <- fect:::inter_fe_ub(Y, ini$Y0, X, I, as.matrix(0), b0, 0L, 1L, 1e-5, 5000L, NULL)
  b <- coef(.ab0_ref(s))
  expect_lt(out$niter, 5000)
  ## x1 and x2 are equal on every cell, so only their sum is identified
  expect_equal(sum(out$beta[1:2]), unname(b["x1"]), tolerance = 1e-4)
  expect_equal(as.numeric(out$beta[3:6]), unname(b[c("x3", "x4", "x5", "x6")]), tolerance = 1e-4)
})
