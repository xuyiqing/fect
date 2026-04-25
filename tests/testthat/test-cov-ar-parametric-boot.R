# Test: parametric bootstrap on unbalanced AR(1) panels now preserves within-unit
# serial correlation (fix at R/boot.R:900, 904; cov.ar = 0 -> cov.ar = TT - 1L).
#
# Strategy: run fect twice on the same DGP, once with the fix (branch HEAD)
# and once with res.vcov monkey-patched to force cov.ar = 0 (simulating the
# pre-fix behavior). Compare the bootstrap sampling SD of the post-period
# time-averaged ATT — this is the quantity affected by within-unit serial
# correlation, because cross-period correlation in residuals shows up only
# when the estimator aggregates across time.
#
# Observed magnitudes (empirical, at nboots=200, N_co=30, N_tr=6, T=40, T_0=20):
#   rho=0.0 : ratio ~= 1.0  (no serial correlation to miss)
#   rho=0.8 : ratio ~= 1.57 (fix restores cross-period cov; dampened by IFE refit)
#   rho=0.9 : ratio ~= 2.17
# The naive asymptotic VIF sqrt((1+rho)/(1-rho)) = 3 / 4.36 is an upper bound;
# the IFE factor projection absorbs part of the serial pattern.

# ---- helpers ----

.gen_ar1 <- function(T, rho, n, seed) {
  set.seed(seed)
  out <- matrix(0, T, n)
  e <- matrix(rnorm(T * n), T, n)
  if (abs(rho) < 1) {
    out[1, ] <- e[1, ] / sqrt(1 - rho^2)
  } else {
    out[1, ] <- e[1, ]
  }
  if (T > 1) {
    for (t in 2:T) out[t, ] <- rho * out[t - 1, ] + e[t, ]
  }
  out
}

.make_dgp <- function(N_co, N_tr, T, T_0, rho, miss_rate, seed) {
  N <- N_co + N_tr
  D_mat <- matrix(0, T, N)
  D_mat[(T_0 + 1):T, (N_co + 1):N] <- 1

  set.seed(seed + 1)
  alpha <- matrix(rep(rnorm(N), each = T), T, N)
  lambda <- matrix(rep(rnorm(T), times = N), T, N)
  eps <- .gen_ar1(T, rho, N, seed + 2)
  Y_mat <- alpha + lambda + eps + 0.5 * D_mat

  if (miss_rate > 0) {
    set.seed(seed + 3)
    miss <- matrix(runif(T * N) < miss_rate, T, N)
    # never blank treated post-period
    miss[(T_0 + 1):T, (N_co + 1):N] <- FALSE
    Y_mat[miss] <- NA
  }

  data.frame(
    id   = rep(1:N, each = T),
    time = rep(1:T, times = N),
    Y    = as.vector(Y_mat),
    D    = as.vector(D_mat)
  )
}

.fit_one <- function(df, nboots, seed) {
  fect(
    Y ~ D, data = df, index = c("id", "time"),
    method = "ife", force = "two-way",
    CV = FALSE, r = 1, se = TRUE,
    vartype = "parametric", nboots = nboots, parallel = FALSE, seed = seed,
    time.component.from = "nevertreated"
  )
}

.with_forced_cov_ar_zero <- function(expr) {
  orig <- fect:::res.vcov
  # Capture `orig` in a closure so the patched function can find it when
  # dispatched from within the fect namespace (the raw `orig` binding is not
  # visible from that scope).
  .mk <- local({
    .orig <- orig
    function() function(res, cov.ar = 1) .orig(res = res, cov.ar = 0)
  })
  patched <- .mk()
  assignInNamespace("res.vcov", patched, ns = "fect")
  on.exit(assignInNamespace("res.vcov", orig, ns = "fect"), add = TRUE)
  force(expr)
}

.se_post_avg <- function(fit) {
  # SD of the time-averaged bootstrap ATT — this is the quantity where
  # within-unit serial correlation in residuals matters.
  b <- as.vector(fit$att.avg.boot)
  sd(b, na.rm = TRUE)
}

# ---- tests ----

test_that("T2a: unbalanced AR(1), rho=0.8, fix inflates att.avg bootstrap SD", {
  skip_on_cran()
  df <- .make_dgp(N_co = 30, N_tr = 6, T = 40, T_0 = 20, rho = 0.8,
                  miss_rate = 0.20, seed = 42)
  fit_new <- .fit_one(df, nboots = 200, seed = 42)
  fit_old <- .with_forced_cov_ar_zero(.fit_one(df, nboots = 200, seed = 42))
  se_new <- .se_post_avg(fit_new)
  se_old <- .se_post_avg(fit_old)
  ratio <- se_new / se_old
  message(sprintf("[T2a] att.avg SD  new=%.4f  old=%.4f  ratio=%.3f", se_new, se_old, ratio))
  # Expect ratio clearly above 1: fix restores the cross-period covariance
  # component that the old code zeroed. IFE refit damps the naive
  # sqrt((1+rho)/(1-rho))=3 asymptote to ~1.5-1.7 in finite samples.
  expect_gt(ratio, 1.25)
  expect_lt(ratio, 3.5)   # guard against runaway inflation too
})

test_that("T2b: balanced-panel parametric-boot SE is unchanged by the fix", {
  skip_on_cran()
  df <- .make_dgp(N_co = 30, N_tr = 6, T = 20, T_0 = 12, rho = 0.8,
                  miss_rate = 0.0, seed = 42)
  expect_true(!anyNA(df$Y))  # truly balanced
  fit_new <- .fit_one(df, nboots = 200, seed = 42)
  fit_old <- .with_forced_cov_ar_zero(.fit_one(df, nboots = 200, seed = 42))
  # Balanced-panel branch is the `else` side of `if (0 %in% I)` at
  # R/boot.R:906; uses direct vector resampling and never calls res.vcov.
  # Must match bit-for-bit.
  se_new <- .se_post_avg(fit_new)
  se_old <- .se_post_avg(fit_old)
  message(sprintf("[T2b] balanced att.avg SD new=%.6f old=%.6f diff=%.3e",
                  se_new, se_old, se_new - se_old))
  expect_equal(se_new, se_old, tolerance = 1e-12)

  # Also per-period SE bitwise-identical
  post_rows <- as.numeric(rownames(fit_new$est.att)) >= 1
  expect_equal(fit_new$est.att[post_rows, "S.E."],
               fit_old$est.att[post_rows, "S.E."],
               tolerance = 1e-12)
})

test_that("T2c: unbalanced rho=0, fix leaves att.avg bootstrap SD ~unchanged", {
  skip_on_cran()
  df <- .make_dgp(N_co = 30, N_tr = 6, T = 40, T_0 = 20, rho = 0.0,
                  miss_rate = 0.20, seed = 42)
  fit_new <- .fit_one(df, nboots = 200, seed = 42)
  fit_old <- .with_forced_cov_ar_zero(.fit_one(df, nboots = 200, seed = 42))
  se_new <- .se_post_avg(fit_new)
  se_old <- .se_post_avg(fit_old)
  ratio <- se_new / se_old
  message(sprintf("[T2c] rho=0 att.avg SD new=%.4f old=%.4f ratio=%.3f",
                  se_new, se_old, ratio))
  # No serial correlation -> no cross-period covariance to restore. Any
  # observed ratio is MC noise between two bootstrap runs with the same
  # seed but slightly different covariance-estimation code paths.
  expect_gt(ratio, 0.80)
  expect_lt(ratio, 1.25)
})
