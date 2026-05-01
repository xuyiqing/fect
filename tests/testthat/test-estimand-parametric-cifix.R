## ---------------------------------------------------------------------------
## v2.4.2 post-hoc fix: parametric bootstrap CI calibration
##
## The pre-fix code passed eff.boot (H0-centered) directly to .compute_ci()
## for all ci.methods except "normal".  This caused percentile/basic/bc/bca
## CIs to be anchored near zero instead of near the point estimate, giving
## 0% coverage for non-zero ATTs.
##
## The fix (R/po-estimands.R): for vartype="parametric", shift att_b by
##   att_b <- att_b - mean(att_b) + estimate
## before calling .compute_ci().  This re-centers the distribution while
## preserving sd(), so "normal" CIs are unaffected.
##
## Tests in this file (from test-spec.md §3-§8):
##   P-INV-1  normal CI satisfies Wald structure exactly
##   P-INV-2  point estimate unchanged by fix
##   P-INV-3  se is identical for normal and basic (shift preserves sd)
##   P-INV-4  .compute_ci() with mock distribution: shifted vs unshifted
##   P-INV-5  bca jackknife path produces non-degenerate CI
##   P-EDGE-2 vartype="none" still returns NA for SE/CI
##   P-EDGE-3 bca on very small nboots does not hard-error
##   P-EDGE-4 shift does not affect vartype="bootstrap" fits
##   P-REG-1  pre-fix degenerate coverage is eliminated (smoke test, 5 reps)
##   P-COV-1  coverage >=0.85 for {basic,percentile,bca,normal}, DGP-A, fe
##   P-COV-2  coverage >=0.85 for bca, DGP-G, gsynth
##   P-COV-3  coverage >=0.75 for bca on aptt, DGP-B-positive
##   P-WIDTH-1 CI widths within 15% of bootstrap widths
## ---------------------------------------------------------------------------

## ---- DGP helpers ------------------------------------------------------------

## DGP-A: two-way FE, additive, true ATT = 3.0
make_panel_A <- function(seed) {
  set.seed(seed)
  N <- 40; TT <- 20; T0 <- 12; Ntr <- 12
  alpha_i <- rnorm(N, 0, 2)
  xi_t    <- rnorm(TT, 0, 1)
  D <- matrix(0L, TT, N); D[(T0 + 1):TT, 1:Ntr] <- 1L
  eps <- matrix(rnorm(N * TT, 0, 1), TT, N)
  Y <- outer(xi_t, rep(1, N)) + outer(rep(1, TT), alpha_i) + 3.0 * D + eps
  data.frame(id = rep(1:N, each = TT), time = rep(1:TT, N),
             Y = as.vector(Y), D = as.vector(D))
}

## DGP-G: interactive fixed effects with factor structure, true ATT = 3.0
make_panel_G <- function(seed, r = 2) {
  set.seed(seed)
  N <- 40; TT <- 20; T0 <- 12; Ntr <- 12
  Lambda <- matrix(rnorm(N * r, 0, 1), N, r)
  F_t    <- matrix(rnorm(TT * r, 0, 1), TT, r)
  eps    <- matrix(rnorm(N * TT, 0, 0.5), TT, N)
  Y0     <- F_t %*% t(Lambda) + eps
  D      <- matrix(0L, TT, N); D[(T0 + 1):TT, 1:Ntr] <- 1L
  Y      <- Y0 + 3.0 * D
  data.frame(id = rep(1:N, each = TT), time = rep(1:TT, N),
             Y = as.vector(Y), D = as.vector(D))
}

## DGP-B-positive: multiplicative effect, true APTT ~ 0.3, true log.ATT ~ log(1.3)
make_panel_B_pos <- function(seed) {
  set.seed(seed)
  N <- 40; TT <- 20; T0 <- 12; Ntr <- 12
  alpha_i <- rnorm(N, 0, 0.5)
  xi_t    <- rnorm(TT, 0, 0.3)
  eps     <- matrix(rnorm(N * TT, 0, 0.5), TT, N)
  Y0      <- 20 + outer(xi_t, rep(1, N)) + outer(rep(1, TT), alpha_i) + eps
  D       <- matrix(0L, TT, N); D[(T0 + 1):TT, 1:Ntr] <- 1L
  Y1      <- Y0; Y1[D == 1L] <- Y0[D == 1L] * 1.3
  Y_obs   <- pmax(ifelse(D == 1L, Y1, Y0), 1)
  treated <- which(D == 1L)
  list(
    df = data.frame(id = rep(1:N, each = TT), time = rep(1:TT, N),
                    Y = as.vector(Y_obs), D = as.vector(D)),
    true_aptt   = mean((Y1[treated] - Y0[treated]) / Y0[treated]),
    true_logatt = mean(log(Y1[treated] / Y0[treated]))
  )
}

## ---- Shared parametric fit helper -------------------------------------------
## Uses same fixture as test-estimand-parametric.R (.make_parametric_fit())
## but we define our own per-DGP helper so we can control the seed cleanly.

.make_parafix_fit <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    skip_on_cran()
    e <- new.env()
    data(sim_linear, package = "fect", envir = e)
    set.seed(42)
    suppressMessages(
      fect(Y ~ D, data = e$sim_linear, index = c("id", "time"),
           method = "ife", force = "two-way", se = TRUE,
           nboots = 30, r = 2, CV = FALSE, keep.sims = TRUE,
           vartype = "parametric",
           time.component.from = "nevertreated",
           parallel = FALSE)
    )
  }
})

## ---- P-INV-1: normal CI satisfies Wald structure ----------------------------

test_that("P-INV-1: normal CI is byte-stable (Wald structure: ci = est +/- z*se)", {
  skip_on_cran()
  fit <- .make_parafix_fit()
  res <- fect::estimand(fit, "att", "overall", window = c(1, 5),
                        ci.method = "normal")
  z <- stats::qnorm(0.975)
  expect_equal(res$ci.lo, res$estimate - z * res$se, tolerance = 1e-10,
               label = "ci.lo = estimate - z*se")
  expect_equal(res$ci.hi, res$estimate + z * res$se, tolerance = 1e-10,
               label = "ci.hi = estimate + z*se")
  ## Width = 2*z*se
  expect_equal(res$ci.hi - res$ci.lo, 2 * z * res$se, tolerance = 1e-10,
               label = "CI width = 2*z*se")
  ## Non-degenerate
  expect_gt(res$ci.hi - res$ci.lo, 1e-6, label = "CI width > 1e-6")
})

## ---- P-INV-2: point estimate is unchanged by fix ----------------------------

test_that("P-INV-2: point estimate matches mean(fit$eff) at treated cells (fix-invariant)", {
  skip_on_cran()
  fit <- .make_parafix_fit()
  res <- fect::estimand(fit, "att", "overall", window = c(1, 5),
                        ci.method = "bca")
  ## Manual point estimate: mean over treated post-treatment cells with T.on in 1..5
  mask <- fit$D.dat == 1 &
          !is.na(fit$T.on) &
          fit$T.on >= 1 &
          fit$T.on <= 5
  manual_est <- mean(fit$eff[mask], na.rm = TRUE)
  expect_equal(res$estimate, manual_est, tolerance = 1e-10,
               label = "estimate == mean(eff[treated cells, T.on 1..5])")
})

## ---- P-INV-3: se is identical for normal and basic (shift preserves sd) -----

test_that("P-INV-3: se is identical for ci.method=normal and ci.method=basic on parametric fit", {
  skip_on_cran()
  fit <- .make_parafix_fit()
  res_normal <- fect::estimand(fit, "att", "overall", window = c(1, 5),
                               ci.method = "normal")
  res_basic  <- fect::estimand(fit, "att", "overall", window = c(1, 5),
                               ci.method = "basic")
  expect_equal(res_normal$se, res_basic$se, tolerance = 1e-10,
               label = "se(normal) == se(basic): shift preserves sd()")
})

## ---- P-INV-4: .compute_ci() mock test: shifted vs unshifted distribution ----

test_that("P-INV-4: .compute_ci() shifted distribution contains estimate; unshifted does not", {
  skip_on_cran()
  set.seed(100)
  boot_h0 <- rnorm(500, mean = 0, sd = 0.15)   ## H0-centered
  estimate <- 2.5
  boot_centered <- boot_h0 - mean(boot_h0) + estimate   ## shifted to point estimate

  ci_bc_shifted   <- fect:::.compute_ci(estimate, boot_centered, "bc", 0.95)
  ci_bc_unshifted <- fect:::.compute_ci(estimate, boot_h0,       "bc", 0.95)

  ## Shifted CI contains the estimate
  expect_lte(ci_bc_shifted$ci.lo, estimate,
             label = "shifted ci.lo <= estimate")
  expect_gte(ci_bc_shifted$ci.hi, estimate,
             label = "shifted ci.hi >= estimate")

  ## Unshifted (H0-centered) CI does NOT contain estimate=2.5
  expect_false(
    ci_bc_unshifted$ci.lo <= estimate && estimate <= ci_bc_unshifted$ci.hi,
    label = "H0-centered CI should not contain estimate=2.5"
  )

  ## Shifted CI lower bound is close to estimate - 1.96*sd (Wald-like)
  expected_lo_approx <- estimate - 1.96 * sd(boot_h0)
  expect_equal(ci_bc_shifted$ci.lo, expected_lo_approx,
               tolerance = 0.15,  ## 15% tolerance — bc differs slightly from Wald
               label = "shifted ci.lo ≈ estimate - 1.96*sd(boot_h0)")

  ## Shifted CI is non-degenerate (width > 1e-6)
  expect_gt(ci_bc_shifted$ci.hi - ci_bc_shifted$ci.lo, 1e-6,
            label = "shifted CI width > 1e-6")

  ## Unshifted CI is degenerate (width << 0.01): confirms the pre-fix pathology
  expect_lt(ci_bc_unshifted$ci.hi - ci_bc_unshifted$ci.lo, 0.01,
            label = "H0-centered bc CI is near-degenerate (pre-fix pathology)")
})

## ---- P-INV-5: bca jackknife path is not broken by the shift -----------------

test_that("P-INV-5: bca on parametric fit produces non-degenerate non-NA CI", {
  skip_on_cran()
  fit <- .make_parafix_fit()
  res <- fect::estimand(fit, "att", "overall", window = c(1, 5),
                        ci.method = "bca")
  expect_s3_class(res, "data.frame")
  expect_true(!is.na(res$ci.lo),  label = "bca ci.lo is not NA")
  expect_true(!is.na(res$ci.hi),  label = "bca ci.hi is not NA")
  expect_true(!is.na(res$se),     label = "bca se is not NA")
  ## Non-degenerate CI
  expect_gt(res$ci.hi - res$ci.lo, 0.01,
            label = "bca CI width > 0.01")
})

## ---- P-EDGE-2: vartype="none" still returns NA for SE/CI --------------------

test_that("P-EDGE-2: vartype='none' returns NA se/ci.lo/ci.hi (unchanged by fix)", {
  skip_on_cran()
  fit <- .make_parafix_fit()
  res <- fect::estimand(fit, "att", "overall", window = c(1, 5),
                        vartype = "none")
  expect_true(!is.na(res$estimate), label = "estimate is not NA with vartype='none'")
  expect_true(is.na(res$se),        label = "se is NA with vartype='none'")
  expect_true(is.na(res$ci.lo),     label = "ci.lo is NA with vartype='none'")
  expect_true(is.na(res$ci.hi),     label = "ci.hi is NA with vartype='none'")
})

## ---- P-EDGE-3: bca on tiny nboots does not hard-error -----------------------

test_that("P-EDGE-3: bca with nboots=5 does not hard-error (returns NA CI, not crash)", {
  skip_on_cran()
  ## Fit with very small nboots; bca guard handles sum(valid)==0 or <2 jackknife points
  d <- make_panel_A(seed = 77)
  set.seed(77)
  fit_tiny <- suppressMessages(
    fect(Y ~ D, data = d, index = c("id", "time"),
         method = "fe", force = "two-way", se = TRUE,
         nboots = 5, CV = FALSE, keep.sims = TRUE,
         vartype = "parametric",
         time.component.from = "nevertreated",
         parallel = FALSE)
  )
  ## Should NOT throw an error — guard returns NA CI if needed
  expect_no_error(
    suppressWarnings(
      fect::estimand(fit_tiny, "att", "overall", window = c(1, 8),
                     ci.method = "bca")
    )
  )
  res <- suppressWarnings(
    fect::estimand(fit_tiny, "att", "overall", window = c(1, 8),
                   ci.method = "bca")
  )
  expect_s3_class(res, "data.frame")
  ## estimate is always computable even with tiny nboots
  expect_true(!is.na(res$estimate), label = "estimate non-NA even with tiny nboots")
})

## ---- P-EDGE-4: shift does NOT affect vartype="bootstrap" fits ---------------

test_that("P-EDGE-4: shift guard is inactive for vartype='bootstrap' (fits unaffected)", {
  skip_on_cran()
  d <- make_panel_A(seed = 42)
  set.seed(42)
  fit_boot <- suppressMessages(
    fect(Y ~ D, data = d, index = c("id", "time"),
         method = "fe", force = "two-way", se = TRUE,
         nboots = 50, CV = FALSE, keep.sims = TRUE,
         vartype = "bootstrap",
         parallel = FALSE)
  )
  expect_equal(fit_boot$vartype, "bootstrap",
               label = "fit$vartype is 'bootstrap'")
  ## bc on bootstrap: should produce a sensible CI (no shift applied)
  res <- fect::estimand(fit_boot, "att", "overall", window = c(1, 8),
                        ci.method = "bc")
  expect_s3_class(res, "data.frame")
  expect_true(!is.na(res$ci.lo), label = "bootstrap bc ci.lo non-NA")
  expect_true(!is.na(res$ci.hi), label = "bootstrap bc ci.hi non-NA")
  expect_gt(res$ci.hi - res$ci.lo, 1e-6, label = "bootstrap bc CI non-degenerate")
})

## ---- P-REG-1: pre-fix degenerate coverage is eliminated (5-rep smoke test) --

test_that("P-REG-1: bca CI contains true ATT=3 on parametric fit (smoke: 5 reps, DGP-A)", {
  skip_on_cran()
  true_att <- 3.0
  n_reps   <- 5L
  in_ci    <- logical(n_reps)
  for (r in seq_len(n_reps)) {
    d <- make_panel_A(seed = 42 + r - 1L)
    set.seed(42 + r - 1L + 5000L)
    fit <- suppressMessages(
      fect(Y ~ D, data = d, index = c("id", "time"),
           method = "fe", force = "two-way", se = TRUE,
           nboots = 200, CV = FALSE, keep.sims = TRUE,
           vartype = "parametric",
           time.component.from = "nevertreated",
           parallel = FALSE)
    )
    res <- fect::estimand(fit, "att", "overall", window = c(1, 8),
                          ci.method = "bca")
    in_ci[r] <- !is.na(res$ci.lo) &&
                res$ci.lo <= true_att &&
                true_att  <= res$ci.hi
  }
  ## At least 3 of 5 reps must cover true ATT (catastrophic failure was 0/5 pre-fix)
  expect_gte(sum(in_ci), 3L,
             label = paste0("bca coverage >= 3/5 reps (got ", sum(in_ci), "/5)"))
})

test_that("P-REG-1b: basic CI width > 0.1 on parametric fit (non-degenerate, DGP-A seed=42)", {
  skip_on_cran()
  d <- make_panel_A(seed = 42)
  set.seed(5042)
  fit <- suppressMessages(
    fect(Y ~ D, data = d, index = c("id", "time"),
         method = "fe", force = "two-way", se = TRUE,
         nboots = 200, CV = FALSE, keep.sims = TRUE,
         vartype = "parametric",
         time.component.from = "nevertreated",
         parallel = FALSE)
  )
  res_basic <- fect::estimand(fit, "att", "overall", window = c(1, 8),
                              ci.method = "basic")
  res_bca   <- fect::estimand(fit, "att", "overall", window = c(1, 8),
                              ci.method = "bca")
  ## Pre-fix: basic CI was ~[5.7, 6.3] for true ATT=3 → 0% coverage
  ## Post-fix: basic CI width must be > 0.1
  expect_gt(res_basic$ci.hi - res_basic$ci.lo, 0.1,
            label = "basic CI width > 0.1 (not degenerate)")
  ## bca CI must be non-degenerate
  expect_gt(res_bca$ci.hi - res_bca$ci.lo, 0.1,
            label = "bca CI width > 0.1 (not degenerate)")
  ## basic CI should contain or be near true ATT=3
  true_att <- 3.0
  half_width <- (res_basic$ci.hi - res_basic$ci.lo) / 2
  dist_to_ci <- max(0, res_basic$ci.lo - true_att, true_att - res_basic$ci.hi)
  expect_lte(dist_to_ci, half_width,
             label = "basic CI is within one half-width of true ATT=3")
})

## ---- All 5 ci.methods produce non-NA, non-degenerate CIs --------------------

test_that("All 5 ci.methods produce non-NA, non-degenerate CIs on parametric fit (DGP-A)", {
  skip_on_cran()
  d <- make_panel_A(seed = 42)
  set.seed(5042)
  fit <- suppressMessages(
    fect(Y ~ D, data = d, index = c("id", "time"),
         method = "fe", force = "two-way", se = TRUE,
         nboots = 200, CV = FALSE, keep.sims = TRUE,
         vartype = "parametric",
         time.component.from = "nevertreated",
         parallel = FALSE)
  )
  for (m in c("basic", "percentile", "bc", "bca", "normal")) {
    res <- fect::estimand(fit, "att", "overall", window = c(1, 8),
                          ci.method = m)
    expect_s3_class(res, "data.frame")
    expect_true(!is.na(res$ci.lo),
                label = paste(m, "ci.lo is not NA"))
    expect_true(!is.na(res$ci.hi),
                label = paste(m, "ci.hi is not NA"))
    expect_gt(res$ci.hi - res$ci.lo, 1e-6,
              label = paste(m, "CI width > 1e-6 (not degenerate)"))
  }
})

## ---- gsynth method: all 5 ci.methods produce reasonable CIs -----------------

test_that("gsynth + parametric: all 5 ci.methods produce non-NA, non-degenerate CIs (DGP-G)", {
  skip_on_cran()
  d <- make_panel_G(seed = 42, r = 2)
  set.seed(5042)
  fit_g <- suppressMessages(
    fect(Y ~ D, data = d, index = c("id", "time"),
         method = "gsynth", r = 2, se = TRUE,
         nboots = 200, CV = FALSE, keep.sims = TRUE,
         vartype = "parametric",
         parallel = FALSE)
  )
  expect_equal(fit_g$vartype, "parametric",
               label = "gsynth fit$vartype is 'parametric'")
  for (m in c("basic", "percentile", "bc", "bca", "normal")) {
    res <- fect::estimand(fit_g, "att", "overall", window = c(1, 8),
                          ci.method = m)
    expect_s3_class(res, "data.frame")
    expect_true(!is.na(res$ci.lo),
                label = paste("gsynth", m, "ci.lo non-NA"))
    expect_true(!is.na(res$ci.hi),
                label = paste("gsynth", m, "ci.hi non-NA"))
    expect_gt(res$ci.hi - res$ci.lo, 1e-6,
              label = paste("gsynth", m, "CI non-degenerate"))
    ## Estimate shifted to CI: CI should bracket estimate (or near)
    ci_width <- res$ci.hi - res$ci.lo
    dist_from_est <- max(0, res$ci.lo - res$estimate,
                            res$estimate - res$ci.hi)
    expect_lte(dist_from_est, ci_width,
               label = paste("gsynth", m, "CI within one width of estimate"))
  }
})

## ---- P-COV-1: ATT coverage >= 0.85 for {basic,percentile,bca,normal}, DGP-A --

test_that("P-COV-1: ATT coverage >= 0.85 for basic/percentile/bca/normal (DGP-A, 100 reps)", {
  skip_on_cran()
  true_att <- 3.0
  n_reps   <- 100L
  methods  <- c("basic", "percentile", "bca", "normal")

  coverage <- setNames(numeric(length(methods)), methods)

  for (r in seq_len(n_reps)) {
    d    <- make_panel_A(seed = 1000L + r)
    set.seed(1000L + r + 5000L)
    fit  <- suppressMessages(
      fect(Y ~ D, data = d, index = c("id", "time"),
           method = "fe", force = "two-way", se = TRUE,
           nboots = 200, CV = FALSE, keep.sims = TRUE,
           vartype = "parametric",
           time.component.from = "nevertreated",
           parallel = FALSE)
    )
    for (m in methods) {
      res <- fect::estimand(fit, "att", "overall", window = c(1, 8),
                            ci.method = m)
      if (!is.na(res$ci.lo) && res$ci.lo <= true_att && true_att <= res$ci.hi) {
        coverage[m] <- coverage[m] + 1L
      }
    }
  }
  coverage <- coverage / n_reps

  for (m in methods) {
    expect_gte(coverage[m], 0.85,
               label = paste0("P-COV-1: ", m, " coverage >= 0.85 (got ",
                              round(coverage[m], 3), ")"))
  }

  ## Anti-regression: any method with < 0.50 coverage is a catastrophic failure
  for (m in methods) {
    expect_gte(coverage[m], 0.50,
               label = paste0("P-COV-1 catastrophic: ", m, " coverage < 0.50"))
  }
})

## ---- P-COV-2: bca coverage >= 0.85, DGP-G, gsynth --------------------------

test_that("P-COV-2: ATT bca coverage >= 0.85 for gsynth+parametric (DGP-G, 100 reps)", {
  skip_on_cran()
  true_att <- 3.0
  n_reps   <- 100L
  in_ci    <- logical(n_reps)

  for (r in seq_len(n_reps)) {
    d    <- make_panel_G(seed = 1000L + r, r = 2)
    set.seed(1000L + r + 5000L)
    fit  <- suppressMessages(
      fect(Y ~ D, data = d, index = c("id", "time"),
           method = "gsynth", r = 2, se = TRUE,
           nboots = 200, CV = FALSE, keep.sims = TRUE,
           vartype = "parametric",
           parallel = FALSE)
    )
    res <- fect::estimand(fit, "att", "overall", window = c(1, 8),
                          ci.method = "bca")
    in_ci[r] <- !is.na(res$ci.lo) &&
                res$ci.lo <= true_att &&
                true_att  <= res$ci.hi
  }
  cov_bca <- mean(in_ci)
  expect_gte(cov_bca, 0.85,
             label = paste0("P-COV-2: gsynth bca coverage >= 0.85 (got ",
                            round(cov_bca, 3), ")"))
})

## ---- P-COV-3: APTT bca coverage >= 0.75, DGP-B-positive --------------------

test_that("P-COV-3: APTT bca coverage >= 0.75 (DGP-B-positive, 100 reps)", {
  skip_on_cran()
  n_reps <- 100L
  in_ci  <- logical(n_reps)

  for (r in seq_len(n_reps)) {
    dgp  <- make_panel_B_pos(seed = 1000L + r)
    true_aptt <- dgp$true_aptt
    set.seed(1000L + r + 5000L)
    fit  <- suppressMessages(
      fect(Y ~ D, data = dgp$df, index = c("id", "time"),
           method = "fe", force = "two-way", se = TRUE,
           nboots = 200, CV = FALSE, keep.sims = TRUE,
           vartype = "parametric",
           time.component.from = "nevertreated",
           parallel = FALSE)
    )
    res_aptt <- suppressMessages(
      fect::estimand(fit, "aptt", "event.time", ci.method = "bca")
    )
    ## Check coverage at event.time = 1
    et1 <- res_aptt[res_aptt$event.time == 1, ]
    if (nrow(et1) == 1L && !is.na(et1$ci.lo)) {
      in_ci[r] <- et1$ci.lo <= true_aptt && true_aptt <= et1$ci.hi
    }
  }
  cov_aptt <- mean(in_ci)
  expect_gte(cov_aptt, 0.75,
             label = paste0("P-COV-3: APTT bca coverage >= 0.75 (got ",
                            round(cov_aptt, 3), ")"))
})

## ---- P-WIDTH-1: parametric CI widths within 15% of bootstrap widths ---------

test_that("P-WIDTH-1: parametric CI widths within 15% of bootstrap widths (DGP-A, seed=42)", {
  skip_on_cran()
  d <- make_panel_A(seed = 42)

  set.seed(42)
  fit_par  <- suppressMessages(
    fect(Y ~ D, data = d, index = c("id", "time"),
         method = "fe", force = "two-way", se = TRUE,
         nboots = 200, CV = FALSE, keep.sims = TRUE,
         vartype = "parametric",
         time.component.from = "nevertreated",
         parallel = FALSE)
  )
  set.seed(42)
  fit_boot <- suppressMessages(
    fect(Y ~ D, data = d, index = c("id", "time"),
         method = "fe", force = "two-way", se = TRUE,
         nboots = 200, CV = FALSE, keep.sims = TRUE,
         vartype = "bootstrap",
         parallel = FALSE)
  )

  for (m in c("basic", "percentile", "normal")) {
    res_par  <- fect::estimand(fit_par,  "att", "overall", window = c(1, 8),
                               ci.method = m)
    res_boot <- fect::estimand(fit_boot, "att", "overall", window = c(1, 8),
                               ci.method = m)
    w_par  <- res_par$ci.hi  - res_par$ci.lo
    w_boot <- res_boot$ci.hi - res_boot$ci.lo

    ## Width within 15% of bootstrap width (both > 0)
    expect_gt(w_par, 0, label = paste(m, "parametric CI width > 0"))
    expect_gt(w_boot, 0, label = paste(m, "bootstrap CI width > 0"))
    rel_diff <- abs(w_par - w_boot) / w_boot
    expect_lte(rel_diff, 0.15,
               label = paste0("P-WIDTH-1: ", m,
                              " parametric vs bootstrap width within 15% (rel_diff = ",
                              round(rel_diff, 3), ")"))
  }
})
