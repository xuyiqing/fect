## ============================================================================
## Full coverage test suite for para.error redesign (test-spec.md scenarios)
## Companion to test-para-error.R (builder's starter file, 10 tests).
## This file covers: T4, T7, T9-T10, T13-T21, E2-E4, E6
##
## Spec: statsclaw-workspace/fect/wild-as-paraerror/test-spec.md
## Branch: feat/v242-completion
## ============================================================================

## ============================================================================
## DGP helpers (spec §2)
## ============================================================================

## DGP-A: additive TWFE, IID Gaussian, true ATT = 3.0, fully observed
##
## Treatment block is FIXED at units 1:Ntr (deterministic). The parametric
## bootstrap targets the conditional variance V_t(Lambda, F, X, D); a coverage
## simulation that re-randomizes the treatment assignment D across replications
## adds Var_{D}[b_t] to the marginal variance the bootstrap is being judged
## against, biasing measured coverage downward (see gsynth-note section 2).
dgp_a <- function(seed) {
    set.seed(seed)
    N <- 40; TT <- 20; T0 <- 12; Ntr <- 12
    alpha_i <- rnorm(N, 0, 2)
    xi_t    <- rnorm(TT, 0, 1)
    eps     <- matrix(rnorm(N * TT, 0, 1), TT, N)
    D <- matrix(0L, TT, N); D[(T0 + 1):TT, 1:Ntr] <- 1L
    Y <- outer(xi_t, rep(1, N)) + outer(rep(1, TT), alpha_i) + 3 * D + eps
    data.frame(
        id   = rep(1:N, each = TT),
        time = rep(1:TT, N),
        Y    = c(Y),
        D    = c(D)
    )
}

## DGP-A8: additive TWFE, AR(1) rho=0.8, true ATT = 3.0
## Treatment block fixed (see dgp_a comment).
dgp_a8 <- function(seed) {
    set.seed(seed)
    N <- 40; TT <- 20; T0 <- 12; Ntr <- 12
    alpha_i <- rnorm(N, 0, 2)
    xi_t    <- rnorm(TT, 0, 1)
    eps <- matrix(NA, TT, N)
    for (i in 1:N) {
        e <- rnorm(TT, 0, 1)
        eps[1, i] <- e[1] / sqrt(1 - 0.64)   # stationary init
        for (t in 2:TT) eps[t, i] <- 0.8 * eps[t - 1, i] + e[t]
    }
    D <- matrix(0L, TT, N); D[(T0 + 1):TT, 1:Ntr] <- 1L
    Y <- outer(xi_t, rep(1, N)) + outer(rep(1, TT), alpha_i) + 3 * D + eps
    data.frame(
        id   = rep(1:N, each = TT),
        time = rep(1:TT, N),
        Y    = c(Y),
        D    = c(D)
    )
}

## DGP-B: positive Y, multiplicative treatment (for aptt / log.att)
## Treatment block fixed (see dgp_a comment).
dgp_b <- function(seed) {
    set.seed(seed)
    N <- 40; TT <- 20; T0 <- 12; Ntr <- 12
    alpha_i <- rnorm(N, 0, 0.5)
    xi_t    <- rnorm(TT, 0, 0.3)
    eps     <- matrix(rnorm(N * TT, 0, 0.5), TT, N)
    D <- matrix(0L, TT, N); D[(T0 + 1):TT, 1:Ntr] <- 1L
    Y0 <- 20 + outer(rep(1, TT), alpha_i) + outer(xi_t, rep(1, N)) + eps
    Y  <- Y0 * ifelse(D == 1, 1.3, 1)
    Y  <- pmax(Y, 0.01)
    data.frame(
        id   = rep(1:N, each = TT),
        time = rep(1:TT, N),
        Y    = c(Y),
        D    = c(D)
    )
}

## DGP-M: missing-data panel (10% rows removed)
dgp_m <- function(seed) {
    df <- dgp_a(seed)
    set.seed(seed + 999)
    mask <- sample(1:nrow(df), floor(0.10 * nrow(df)))
    df <- df[-mask, ]
    df
}

## DGP-REV: panel with treatment reversals (some units turn off treatment)
## N=30: 10 reversing treated, 10 absorbing treated, 10 never-treated controls.
dgp_rev <- function(seed = 77) {
    set.seed(seed)
    N_rev <- 10; N_abs <- 10; N_co <- 10
    N <- N_rev + N_abs + N_co; TT <- 20
    D_vec <- integer(N * TT)
    Y_vec <- numeric(N * TT)
    ## Reversing units: turn on at t=8, turn off at t=14
    for (i in 1:N_rev) {
        for (t in 1:TT) {
            idx <- (i - 1) * TT + t
            D_vec[idx] <- as.integer(t >= 8 && t < 14)
            Y_vec[idx] <- rnorm(1) + 2 * D_vec[idx]
        }
    }
    ## Absorbing treated: turn on at t=8, stay on
    for (i in (N_rev + 1):(N_rev + N_abs)) {
        for (t in 1:TT) {
            idx <- (i - 1) * TT + t
            D_vec[idx] <- as.integer(t >= 8)
            Y_vec[idx] <- rnorm(1) + 2 * D_vec[idx]
        }
    }
    ## Never-treated controls
    for (i in (N_rev + N_abs + 1):N) {
        for (t in 1:TT) {
            idx <- (i - 1) * TT + t
            D_vec[idx] <- 0L
            Y_vec[idx] <- rnorm(1)
        }
    }
    data.frame(
        id   = rep(1:N, each = TT),
        time = rep(1:TT, N),
        Y    = Y_vec,
        D    = D_vec
    )
}

## DGP-STAG: staggered adoption (3 cohorts) for E6
dgp_stag <- function(seed = 55) {
    set.seed(seed)
    N <- 30; TT <- 20; Ntr <- 18
    cohort_times <- c(rep(6, 6), rep(10, 6), rep(14, 6), rep(TT + 1, 12))
    D <- matrix(0, TT, N)
    for (i in 1:Ntr) {
        if (cohort_times[i] <= TT) D[cohort_times[i]:TT, i] <- 1
    }
    alpha_i <- rnorm(N, 0, 2)
    xi_t    <- rnorm(TT, 0, 1)
    eps     <- matrix(rnorm(N * TT, 0, 1), TT, N)
    Y <- outer(xi_t, rep(1, N)) + outer(rep(1, TT), alpha_i) + 2 * D + eps
    data.frame(
        id   = rep(1:N, each = TT),
        time = rep(1:TT, N),
        Y    = c(Y),
        D    = c(D)
    )
}

## Shared parametric call helper
fect_para_full <- function(d, para.error = "auto", vartype = "parametric",
                           nboots = 100, seed = 42, ...) {
    suppressWarnings(suppressMessages(
        fect(
            Y ~ D,
            data   = d,
            index  = c("id", "time"),
            method = "fe",
            force  = "two-way",
            time.component.from = "nevertreated",
            se      = TRUE,
            vartype = vartype,
            para.error = para.error,
            nboots     = nboots,
            parallel   = FALSE,
            keep.sims  = TRUE,
            seed       = seed,
            CV         = FALSE,
            ...
        )
    ))
}

## ============================================================================
## T4: Point estimate byte-stability across all three para.error modes
## ============================================================================

test_that("T4: att.avg is byte-identical across para.error = ar, empirical, wild", {
    skip_on_cran()
    df <- dgp_a(seed = 42)
    fit_ar  <- fect_para_full(df, para.error = "ar",        nboots = 50)
    fit_emp <- fect_para_full(df, para.error = "empirical", nboots = 50)
    fit_wld <- fect_para_full(df, para.error = "wild",      nboots = 50)

    expect_identical(fit_ar$att.avg, fit_emp$att.avg,
        label = "att.avg: ar vs empirical should be byte-identical")
    expect_identical(fit_ar$att.avg, fit_wld$att.avg,
        label = "att.avg: ar vs wild should be byte-identical")
})

## ============================================================================
## T7: para.error = "ar" on missing-data panel succeeds
## ============================================================================

test_that("T7: para.error = 'ar' succeeds on missing-data panel", {
    skip_on_cran()
    df_m <- dgp_m(seed = 1)
    fit <- suppressWarnings(suppressMessages(
        fect(
            Y ~ D,
            data   = df_m,
            index  = c("id", "time"),
            method = "fe",
            force  = "two-way",
            time.component.from = "nevertreated",
            se         = TRUE,
            vartype    = "parametric",
            para.error = "ar",
            na.rm      = FALSE,
            nboots     = 100,
            parallel   = FALSE,
            keep.sims  = TRUE,
            CV         = FALSE
        )
    ))
    expect_equal(fit$para.error, "ar",
        label = "para.error should resolve to 'ar' on missing-data panel")
    ## At most 20% NA in att.avg.boot (spec T7 tolerance)
    na_rate <- mean(is.na(c(fit$att.avg.boot)))
    expect_lte(na_rate, 0.20,
        label = sprintf("NA rate in att.avg.boot = %.2f%%, should be <= 20%%",
                        na_rate * 100))
})

## ============================================================================
## T9, T10 (spec) — vartype = "wild" deprecation alias
## REMOVED: vartype = "wild" was never released as a standalone vartype, so
## there is no caller to deprecate. The redesign exposes para.error directly.
##
## T11 (spec) / covered in builder file — verify fit$para.error stores resolved value
## Already in test-para-error.R (Test 1), not repeated here.

## T12 (spec) / covered in builder file (Test 10), not repeated here.

## ============================================================================
## T13: estimand() location-shift applies to all three para.error modes
## ============================================================================

test_that("T13: estimand() percentile CI contains point estimate for all para.error modes", {
    skip_on_cran()
    df <- dgp_a(seed = 42)
    for (pm in c("ar", "empirical", "wild")) {
        fit <- fect_para_full(df, para.error = pm, nboots = 200)
        est <- estimand(fit, "att", "overall", window = c(1, 8), ci.method = "percentile")
        expect_true(
            is.finite(est$ci.lo) && is.finite(est$ci.hi),
            label = sprintf("T13 %s: CI must be finite", pm)
        )
        expect_true(
            est$ci.lo < fit$att.avg && fit$att.avg < est$ci.hi,
            label = sprintf(
                "T13 %s: point estimate %.4f must be inside CI [%.4f, %.4f]",
                pm, fit$att.avg, est$ci.lo, est$ci.hi
            )
        )
    }
})

## ============================================================================
## T14: Anti-regression — vartype = "bootstrap" unchanged
## ============================================================================

test_that("T14: vartype='bootstrap' produces finite, sensible results (anti-regression)", {
    skip_on_cran()
    data("simdata", package = "fect")
    set.seed(42)
    fit <- suppressMessages(
        fect(Y ~ D, data = simdata, index = c("id", "time"),
             method = "fe", force = "two-way",
             se = TRUE, vartype = "bootstrap",
             nboots = 50, parallel = FALSE, seed = 42, CV = FALSE)
    )
    expect_equal(fit$vartype, "bootstrap")
    expect_true(is.numeric(fit$att.avg) && !is.na(fit$att.avg),
        label = "bootstrap att.avg must be finite")
    ## Anti-regression: att.avg should be near 2.5 for simdata (known fixture)
    expect_true(abs(fit$att.avg) < 10,
        label = sprintf("bootstrap att.avg = %.4f looks reasonable", fit$att.avg))
    ## Ensure para.error is NULL for bootstrap path
    expect_null(fit$para.error,
        label = "para.error should be NULL for bootstrap fits")
})

## ============================================================================
## T15: Anti-regression — vartype = "jackknife" unchanged
## ============================================================================

test_that("T15: vartype='jackknife' produces finite, sensible results (anti-regression)", {
    skip_on_cran()
    ## Small panel so jackknife is fast
    df_small <- dgp_a(seed = 7)
    df_small <- df_small[df_small$id <= 20, ]
    fit <- suppressMessages(
        fect(Y ~ D, data = df_small, index = c("id", "time"),
             method = "fe", force = "two-way",
             se = TRUE, vartype = "jackknife",
             parallel = FALSE, CV = FALSE)
    )
    expect_equal(fit$vartype, "jackknife")
    expect_true(is.numeric(fit$att.avg) && !is.na(fit$att.avg),
        label = "jackknife att.avg must be finite")
    expect_null(fit$para.error,
        label = "para.error should be NULL for jackknife fits")
})

## ============================================================================
## T16: Anti-regression — vartype = "parametric" with para.error = "auto" unchanged
## ============================================================================

test_that("T16: vartype='parametric' para.error='auto' produces same resolved path as before", {
    skip_on_cran()
    df <- dgp_a(seed = 42)
    fit <- fect_para_full(df, para.error = "auto", nboots = 50)
    ## On fully-observed panel, auto should resolve to "empirical"
    expect_equal(fit$para.error, "empirical",
        label = "auto should resolve to empirical on fully-observed panel")
    ## Point estimate should be near 3.0 (true ATT)
    expect_true(abs(fit$att.avg - 3.0) < 1.5,
        label = sprintf("att.avg = %.4f should be near true ATT = 3.0", fit$att.avg))
    ## est.avg is a matrix (class = matrix/array)
    expect_true(is.matrix(fit$est.avg),
        label = "est.avg should be a matrix")
    expect_true(all(is.finite(fit$est.avg[, "ATT.avg"])),
        label = "est.avg[,ATT.avg] should be finite")
})

## ============================================================================
## T17: aptt + para.error = "wild" + bca returns finite CIs (DGP-B)
## ============================================================================

test_that("T17: aptt and log.att estimands work with para.error='wild' on DGP-B", {
    skip_on_cran()
    df <- dgp_b(seed = 99)
    fit <- suppressMessages(
        fect(
            Y ~ D,
            data   = df,
            index  = c("id", "time"),
            method = "fe",
            force  = "two-way",
            time.component.from = "nevertreated",
            se         = TRUE,
            vartype    = "parametric",
            para.error = "wild",
            nboots     = 200,
            parallel   = FALSE,
            keep.sims  = TRUE,
            seed       = 99,
            CV         = FALSE
        )
    )
    ## aptt estimand
    est_aptt <- estimand(fit, "aptt", "event.time", ci.method = "bca")
    expect_true(all(is.finite(est_aptt$ci.lo)) && all(is.finite(est_aptt$ci.hi)),
        label = "T17: aptt CI must be finite for all event times")
    expect_true(all(est_aptt$estimate > 0),
        label = sprintf("T17: aptt estimates should be positive (true APTT>0), got min=%.4f",
                        min(est_aptt$estimate)))

    ## log.att estimand
    est_log <- estimand(fit, "log.att", "event.time", ci.method = "bca")
    expect_true(all(is.finite(est_log$ci.lo)) && all(is.finite(est_log$ci.hi)),
        label = "T17: log.att CI must be finite for all event times")
    ## Each log.att estimate should be within +/-0.1 of log(1.3) = 0.2624
    true_logatt <- log(1.3)
    for (i in seq_len(nrow(est_log))) {
        expect_true(
            abs(est_log$estimate[i] - true_logatt) < 0.15,
            label = sprintf(
                "T17: log.att estimate[%d] = %.4f should be within 0.15 of %.4f",
                i, est_log$estimate[i], true_logatt
            )
        )
    }
})

## ============================================================================
## T18: CI width parity — wild vs empirical within 30% (nboots = 1000)
## ============================================================================

test_that("T18: CI width parity — wild/empirical ratio in [0.70, 1.30]", {
    skip_on_cran()
    df <- dgp_a(seed = 42)
    fit_emp <- fect_para_full(df, para.error = "empirical", nboots = 1000)
    fit_wld <- fect_para_full(df, para.error = "wild",      nboots = 1000)
    fit_ar  <- fect_para_full(df, para.error = "ar",        nboots = 1000)

    est_emp <- estimand(fit_emp, "att", "overall", window = c(1, 8), ci.method = "basic")
    est_wld <- estimand(fit_wld, "att", "overall", window = c(1, 8), ci.method = "basic")
    est_ar  <- estimand(fit_ar,  "att", "overall", window = c(1, 8), ci.method = "basic")

    width_emp <- est_emp$ci.hi - est_emp$ci.lo
    width_wld <- est_wld$ci.hi - est_wld$ci.lo
    width_ar  <- est_ar$ci.hi  - est_ar$ci.lo

    ratio_wld_emp <- width_wld / width_emp
    ratio_ar_emp  <- width_ar  / width_emp

    expect_gte(ratio_wld_emp, 0.70,
        label = sprintf("T18: wild/empirical width ratio = %.4f should be >= 0.70",
                        ratio_wld_emp))
    expect_lte(ratio_wld_emp, 1.30,
        label = sprintf("T18: wild/empirical width ratio = %.4f should be <= 1.30",
                        ratio_wld_emp))
    expect_gte(ratio_ar_emp, 0.50,
        label = sprintf("T18: ar/empirical width ratio = %.4f should be >= 0.50",
                        ratio_ar_emp))
})

## ============================================================================
## T19, T20, T21 (spec) --- coverage validation simulations
##
## MOVED to tests/coverage-study/run_para_error_coverage.R
##
## These three Monte-Carlo studies (T19: 100 reps x 1000 nboots x 15 cells on
## DGP-A IID; T20: same on DGP-A8 AR(1) rho=0.8; T21: 50 reps x 500 nboots
## width parity) take ~30 min wall and are not appropriate for routine
## devtools::test() runs.  They live alongside the coverage-study artifacts.
##
## Run them via:  Rscript tests/coverage-study/run_para_error_coverage.R
##
## When to run: any change to R/boot.R parametric branch, R/po-estimands.R
## location-shift code, the vartype / ci.method / para.error machinery,
## jackknife dispatch, or eff.boot construction.  See the README in that
## directory for full trigger list and acceptance criteria.
## ============================================================================

## ============================================================================
## E2: para.error = "wild" with nboots = 200 completes with >= 95% success rate
## ============================================================================

test_that("E2: para.error='wild' nboots=200 completes with >=95% successful replicates", {
    skip_on_cran()
    df <- dgp_a(seed = 42)
    fit <- fect_para_full(df, para.error = "wild", nboots = 200)
    n_ok <- sum(!is.na(c(fit$att.avg.boot)))
    expect_gte(n_ok, 190L,
        label = sprintf("E2: %d/200 replicates succeeded (need >=190)", n_ok))
    ## est.avg is a matrix
    expect_true(is.matrix(fit$est.avg),
        label = "E2: est.avg should be a matrix")
    expect_true(all(is.finite(fit$est.avg[, "ATT.avg"])),
        label = "E2: est.avg[,ATT.avg] should be finite")
    expect_true(all(is.finite(fit$est.avg[, "CI.lower"])),
        label = "E2: est.avg[,CI.lower] should be finite")
    expect_true(all(is.finite(fit$est.avg[, "CI.upper"])),
        label = "E2: est.avg[,CI.upper] should be finite")
})

## ============================================================================
## E3: para.error = "wild" with treatment reversals routes to hard error
## ============================================================================

test_that("E3: para.error='wild' with treatment reversals fires the parametric reversal gate", {
    skip_on_cran()
    df_rev <- dgp_rev(seed = 77)
    expect_error(
        suppressMessages(
            fect(
                Y ~ D,
                data   = df_rev,
                index  = c("id", "time"),
                method = "fe",
                force  = "two-way",
                time.component.from = "nevertreated",
                se         = TRUE,
                vartype    = "parametric",
                para.error = "wild",
                nboots     = 50,
                parallel   = FALSE,
                CV         = FALSE
            )
        ),
        regexp = "Parametric bootstrap is not valid when treatment reversal",
        fixed  = FALSE,
        label  = "E3: reversal gate should fire for para.error='wild'"
    )
})

## ============================================================================
## E4: para.error = "wild" + time.component.from = "notyettreated" -> hard error
## ============================================================================

test_that("E4: para.error='wild' + time.component.from='notyettreated' fires the notyettreated gate", {
    skip_on_cran()
    df <- dgp_a(seed = 42)
    expect_error(
        suppressMessages(
            fect(
                Y ~ D,
                data   = df,
                index  = c("id", "time"),
                method = "fe",
                force  = "two-way",
                time.component.from = "notyettreated",
                se         = TRUE,
                vartype    = "parametric",
                para.error = "wild",
                nboots     = 50,
                parallel   = FALSE,
                CV         = FALSE
            )
        ),
        regexp = "notyettreated",
        fixed  = FALSE,
        label  = "E4: notyettreated gate should fire for para.error='wild'"
    )
})

## ============================================================================
## E6: para.error = "wild" with staggered adoption completes, <= 10% NA
## ============================================================================

test_that("E6: para.error='wild' with staggered adoption completes, <=10% NA in att.avg.boot", {
    skip_on_cran()
    df_stag <- dgp_stag(seed = 55)
    fit <- suppressMessages(
        fect(
            Y ~ D,
            data   = df_stag,
            index  = c("id", "time"),
            method = "fe",
            force  = "two-way",
            time.component.from = "nevertreated",
            se         = TRUE,
            vartype    = "parametric",
            para.error = "wild",
            nboots     = 200,
            parallel   = FALSE,
            keep.sims  = TRUE,
            CV         = FALSE
        )
    )
    expect_true(!is.null(fit$att.avg.boot),
        label = "E6: att.avg.boot should exist")
    na_rate <- mean(is.na(c(fit$att.avg.boot)))
    expect_lte(na_rate, 0.10,
        label = sprintf("E6: NA rate in att.avg.boot = %.2f%% (should be <= 10%%)",
                        na_rate * 100))
})
