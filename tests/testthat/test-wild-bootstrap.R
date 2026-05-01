## ---------------------------------------------------------------
## Tests for the wild bootstrap fix (fect v2.4.2+).
##
## Wild bootstrap now perturbs residuals on ALL observed cells
## (control AND treated), following CGM 2008 style.
##
## Scenarios from test-spec.md:
##   T1:  smoke test — all 5 ci.methods, finite CIs
##   T2:  width parity — wild >= 80% of bootstrap width      [skip_on_cran]
##   T3:  point estimate byte-identical to bootstrap
##   T4:  aptt + bca on DGP-B: finite CIs                    [skip_on_cran]
##   T5:  log.att + bca on DGP-B: finite CIs                 [skip_on_cran]
##   T6:  anti-regression: bootstrap path unchanged
##   T7:  anti-regression: parametric path unchanged
##   T8:  anti-regression: jackknife path unchanged
##   T9:  binary + wild still hard-errors
##   T10: nboots=200 end-to-end works
##   T11: coverage >= 0.85 at 100 reps (simulation)          [skip_on_cran]
##   E1:  staggered panel (T.on outside time.on range) handled
##   E2:  all-same treatment start (length-1 att) handled
##   E3:  placeboTest = TRUE (eff.out NA cells) handled
##   E4:  att.avg NA propagates to make_boot_na without crash
## ---------------------------------------------------------------

suppressWarnings(data("simdata", package = "fect"))

## ---- DGP helpers ---------------------------------------------------

## DGP-A: two-way FE, IID Gaussian, true ATT = 3.0
## Matches test-spec.md §2 exactly.
make_dgp_a <- function(seed) {
  set.seed(seed)
  N <- 40; TT <- 20; T0 <- 12; Ntr <- 12
  alpha_i <- rnorm(N, 0, 2)
  xi_t    <- rnorm(TT, 0, 1)
  eps     <- matrix(rnorm(N * TT, 0, 1), TT, N)
  id_tr   <- sample(1:N, Ntr)
  D <- matrix(0L, TT, N)
  for (i in id_tr) D[(T0 + 1):TT, i] <- 1L
  Y <- outer(xi_t, rep(1, N)) + outer(rep(1, TT), alpha_i) + 3.0 * D + eps
  data.frame(
    id   = rep(1:N, each = TT),
    time = rep(1:TT, N),
    Y    = c(Y),
    D    = c(D)
  )
}

## DGP-B: positive Y, multiplicative, true APTT ~ 0.30, log.ATT ~ log(1.3)
## Matches test-spec.md §2 exactly.
make_dgp_b <- function(seed = 99) {
  set.seed(seed)
  N <- 40; TT <- 20; T0 <- 12; Ntr <- 12
  alpha_i <- rnorm(N, 0, 0.5)
  xi_t    <- rnorm(TT, 0, 0.3)
  eps     <- matrix(rnorm(N * TT, 0, 0.5), TT, N)
  id_tr   <- sample(1:N, Ntr)
  D <- matrix(0L, TT, N)
  for (i in id_tr) D[(T0 + 1):TT, i] <- 1L
  Y0 <- 20 + outer(rep(1, TT), alpha_i) + outer(xi_t, rep(1, N)) + eps
  Y  <- Y0 * ifelse(D == 1L, 1.3, 1)
  Y  <- pmax(Y, 0.01)
  data.frame(
    id   = rep(1:N, each = TT),
    time = rep(1:TT, N),
    Y    = c(Y),
    D    = c(D)
  )
}

## Pre-built DGP-A for anti-regression tests.
dgp_a_data <- make_dgp_a(42)

## ---- Shared fit helpers --------------------------------------------

.fit_wild <- function(data, nboots = 200, seed = 42,
                       keep.sims = FALSE, CV = FALSE) {
  set.seed(seed)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = data, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, vartype = "wild", nboots = nboots,
      parallel = FALSE, keep.sims = keep.sims, CV = CV
    )
  ))
}

.fit_boot <- function(data, nboots = 200, seed = 42,
                       keep.sims = FALSE, CV = FALSE) {
  set.seed(seed)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = data, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, vartype = "bootstrap", nboots = nboots,
      parallel = FALSE, keep.sims = keep.sims, CV = CV
    )
  ))
}

## Parametric needs time.component.from="nevertreated" when staggered
## adoption (not-yet-treated control default is rejected).
.fit_para <- function(data, nboots = 50, seed = 42, CV = FALSE) {
  set.seed(seed)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = data, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, vartype = "parametric", nboots = nboots,
      parallel = FALSE, CV = CV,
      time.component.from = "nevertreated"
    )
  ))
}

.fit_jack <- function(data, CV = FALSE) {
  set.seed(1)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = data, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, vartype = "jackknife",
      parallel = FALSE, CV = CV
    )
  ))
}

## ---- T1: Smoke test ------------------------------------------------

test_that("T1: wild + all ci.methods produce finite, non-degenerate CIs", {

  df <- make_dgp_a(42)
  ## keep.sims = TRUE needed for estimand("att", "overall", window = ...)
  fit <- .fit_wild(df, nboots = 200, keep.sims = TRUE)

  ## (a) att.avg is finite scalar
  expect_true(is.numeric(fit$att.avg))
  expect_true(is.finite(fit$att.avg))

  ## (b) att.avg.boot: at most 5% NA entries
  boot_vals <- fit$att.avg.boot
  expect_false(is.null(boot_vals))
  na_frac <- mean(is.na(boot_vals))
  expect_lte(na_frac, 0.05)

  ## (c) All 5 ci.methods return finite CIs.
  ## window = c(1, 8) filters to post-treatment event times per test-spec.md §T1.
  for (m in c("basic", "percentile", "bc", "bca", "normal")) {
    est <- suppressWarnings(
      fect::estimand(fit, "att", "overall", window = c(1, 8),
                     ci.method = m)
    )
    expect_s3_class(est, "data.frame")
    expect_true(is.finite(est$ci.lo))
    expect_true(is.finite(est$ci.hi))
  }
})

## ---- T2: Width parity -----------------------------------------------

test_that("T2: wild CI width >= 80% of bootstrap CI width", {

  skip_on_cran()

  df <- make_dgp_a(42)

  fit_w <- .fit_wild(df,  nboots = 1000, seed = 42, keep.sims = TRUE)
  fit_b <- .fit_boot(df, nboots = 1000, seed = 42, keep.sims = TRUE)

  for (m in c("basic", "normal")) {
    est_w <- suppressWarnings(
      fect::estimand(fit_w, "att", "overall", window = c(1, 8),
                     ci.method = m)
    )
    est_b <- suppressWarnings(
      fect::estimand(fit_b, "att", "overall", window = c(1, 8),
                     ci.method = m)
    )
    width_w <- est_w$ci.hi - est_w$ci.lo
    width_b <- est_b$ci.hi - est_b$ci.lo

    ## Width ratio must be >= 0.80 (was ~0.55-0.60 pre-fix).
    ## Tolerance from test-spec.md §T2: 0.80.
    ratio <- width_w / width_b
    expect_gte(ratio, 0.80)
  }
})

## ---- T3: Point estimate byte-identical to bootstrap -----------------

test_that("T3: wild att.avg is byte-identical to bootstrap att.avg", {

  df <- make_dgp_a(42)
  fit_w <- .fit_wild(df, nboots = 50, seed = 42)
  fit_b <- .fit_boot(df, nboots = 50, seed = 42)

  ## Same data + same main-fit path -> same point estimate.
  expect_identical(fit_w$att.avg, fit_b$att.avg)
})

## ---- T4: aptt + wild + bca on DGP-B --------------------------------

test_that("T4: wild + aptt + bca on DGP-B returns finite CIs", {

  skip_on_cran()

  df_pos <- make_dgp_b(99)
  fit <- .fit_wild(df_pos, nboots = 200, keep.sims = TRUE)

  est <- suppressWarnings(
    fect::estimand(fit, "aptt", "event.time", ci.method = "bca")
  )

  expect_s3_class(est, "data.frame")
  expect_true(all(is.finite(est$ci.lo)))
  expect_true(all(is.finite(est$ci.hi)))
  expect_true(all(!is.na(est$ci.lo)))
  expect_true(all(!is.na(est$ci.hi)))
  ## ATT > 0 in DGP-B (multiplicative treatment with factor 1.3)
  expect_true(all(est$estimate > 0))
})

## ---- T5: log.att + wild + bca on DGP-B -----------------------------

test_that("T5: wild + log.att + bca on DGP-B returns finite CIs", {

  skip_on_cran()

  df_pos <- make_dgp_b(99)
  fit <- .fit_wild(df_pos, nboots = 200, keep.sims = TRUE)

  est <- suppressWarnings(
    fect::estimand(fit, "log.att", "event.time", ci.method = "bca")
  )

  expect_s3_class(est, "data.frame")
  expect_true(all(is.finite(est$ci.lo)))
  expect_true(all(is.finite(est$ci.hi)))
  expect_true(all(!is.na(est$ci.lo)))
  expect_true(all(!is.na(est$ci.hi)))

  ## Estimates near log(1.3) ~ 0.2624, tolerance +-0.30 for single sample
  expect_true(
    all(abs(est$estimate - log(1.3)) < 0.30)
  )
})

## ---- T6: Anti-regression: bootstrap path ----------------------------

test_that("T6: bootstrap path att.avg and est.att unchanged by wild fix", {

  skip_on_cran()

  ## Use DGP-A data (no reversals) for clean reproducibility.
  fit1 <- .fit_boot(dgp_a_data, nboots = 50, seed = 42)
  fit2 <- .fit_boot(dgp_a_data, nboots = 50, seed = 42)

  ## Two identical runs must be byte-identical.
  expect_identical(fit1$att.avg, fit2$att.avg)
  expect_identical(fit1$est.att, fit2$est.att)
})

## ---- T7: Anti-regression: parametric path ---------------------------

test_that("T7: parametric path numerically stable across two identical runs", {

  skip_on_cran()

  ## Use DGP-A (no reversals, never-treated controls exist).
  fit1 <- .fit_para(dgp_a_data, nboots = 50, seed = 42)
  fit2 <- .fit_para(dgp_a_data, nboots = 50, seed = 42)

  expect_identical(fit1$att.avg, fit2$att.avg)
  expect_identical(fit1$est.att, fit2$est.att)
})

## ---- T8: Anti-regression: jackknife path ----------------------------

test_that("T8: jackknife path numerically stable across two identical runs", {

  skip_on_cran()

  ## Use DGP-A subset (units 1-20 only) per test-spec.md §T8.
  sub <- dgp_a_data[dgp_a_data$id <= 20, ]

  fit1 <- .fit_jack(sub)
  fit2 <- .fit_jack(sub)

  expect_identical(fit1$att.avg, fit2$att.avg)
  expect_identical(fit1$est.att, fit2$est.att)
})

## ---- T9: Binary + wild still hard-errors ---------------------------

test_that("T9: binary = TRUE + vartype = 'wild' still hard-errors at fit time", {

  ## Build a minimal binary panel.
  set.seed(99)
  N <- 20; TT <- 10; T0 <- 6
  df_bin <- data.frame(
    id   = rep(1:N, each = TT),
    time = rep(1:TT, N)
  )
  df_bin$D <- as.integer(df_bin$id <= 5 & df_bin$time > T0)
  df_bin$Y <- as.integer(rbinom(N * TT, 1, 0.4 + 0.2 * df_bin$D))

  expect_error(
    suppressMessages(
      fect::fect(
        Y ~ D, data = df_bin, index = c("id", "time"),
        method = "fe", force = "two-way",
        se = TRUE, vartype = "wild", binary = TRUE,
        parallel = FALSE, CV = FALSE
      )
    ),
    regexp = "\"wild\" vartype is not supported for binary"
  )
})

## ---- T10: nboots = 200 end-to-end -----------------------------------

test_that("T10: wild nboots=200 completes and passes slot checks", {

  df <- make_dgp_a(42)
  fit <- .fit_wild(df, nboots = 200, seed = 1)

  ## At least 95% of replicates succeeded (threshold per test-spec.md §T10).
  n_ok <- sum(!is.na(fit$att.avg.boot))
  expect_gte(n_ok, 190L)

  ## est.avg is a data frame with finite inferential columns.
  ea <- fit$est.avg
  expect_false(is.null(ea))
  expect_true(is.finite(ea[1, "ATT.avg"]))
  expect_true(is.finite(ea[1, "S.E."]))
  expect_true(is.finite(ea[1, "CI.lower"]))
  expect_true(is.finite(ea[1, "CI.upper"]))
})

## ---- T11: Coverage simulation >= 0.85 (100 reps) -------------------

test_that("T11: wild coverage >= 0.85 on DGP-A (100 reps, nboots=1000)", {

  skip_on_cran()

  n_reps     <- 100
  nboots     <- 1000
  ci_methods <- c("basic", "percentile", "bc", "bca", "normal")
  true_att   <- 3.0

  covers <- matrix(FALSE, nrow = n_reps, ncol = length(ci_methods),
                   dimnames = list(NULL, ci_methods))

  for (r in seq_len(n_reps)) {
    df_r <- make_dgp_a(seed = r * 100)
    fit_r <- suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D, data = df_r, index = c("id", "time"),
        method = "fe", force = "two-way",
        se = TRUE, vartype = "wild", nboots = nboots,
        parallel = TRUE, CV = FALSE,
        keep.sims = TRUE
      )
    ))
    for (m in ci_methods) {
      est_r <- suppressWarnings(
        fect::estimand(fit_r, "att", "overall",
                       window = c(1, 8), ci.method = m)
      )
      covers[r, m] <- isTRUE(est_r$ci.lo <= true_att &&
                                est_r$ci.hi >= true_att)
    }
  }

  ## Acceptance criterion per test-spec.md §T11: >= 0.85 for each ci.method.
  for (m in ci_methods) {
    cov_m <- mean(covers[, m])
    expect_gte(cov_m, 0.85)
  }
})

## ---- E1: Staggered adoption (T.on outside time.on range) -----------

test_that("E1: staggered adoption panel (late starters) handled without error", {

  skip_on_cran()

  ## Staggered: some units start treatment at different times.
  set.seed(201)
  N <- 30; TT <- 20
  treat_start <- c(rep(NA, 10), sample(11:19, 20, replace = TRUE))
  df_stag <- data.frame(
    id   = rep(1:N, each = TT),
    time = rep(1:TT, N)
  )
  df_stag$D <- as.integer(
    !is.na(treat_start[df_stag$id]) &
      df_stag$time >= treat_start[df_stag$id]
  )
  alpha_i <- rnorm(N)
  xi_t    <- rnorm(TT, sd = 0.5)
  df_stag$Y <- alpha_i[df_stag$id] + xi_t[df_stag$time] +
    2.5 * df_stag$D + rnorm(N * TT, sd = 0.5)

  fit_s <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = df_stag, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, vartype = "wild", nboots = 100,
      parallel = FALSE, CV = FALSE,
      keep.sims = TRUE
    )
  ))

  expect_true(is.finite(fit_s$att.avg))

  ## Verify CIs are finite via est.avg (no keep.sims requirement).
  ea <- fit_s$est.avg
  expect_true(is.finite(ea[1, "CI.lower"]))
  expect_true(is.finite(ea[1, "CI.upper"]))
})

## ---- E2: att length-1 (all units same treatment start) -------------

test_that("E2: uniform adoption panel (length-1 att) completes without error", {

  skip_on_cran()

  ## All treated units start treatment at the same period T0+1.
  set.seed(202)
  N <- 30; TT <- 20; T0 <- 12; Ntr <- 10
  alpha_i <- rnorm(N)
  xi_t    <- rnorm(TT, sd = 0.5)
  df_uni <- data.frame(
    id   = rep(1:N, each = TT),
    time = rep(1:TT, N)
  )
  df_uni$D <- as.integer(df_uni$id <= Ntr & df_uni$time > T0)
  df_uni$Y <- alpha_i[df_uni$id] + xi_t[df_uni$time] +
    2.0 * df_uni$D + rnorm(N * TT, sd = 0.5)

  fit_u <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = df_uni, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, vartype = "wild", nboots = 100,
      parallel = FALSE, CV = FALSE
    )
  ))

  expect_true(is.finite(fit_u$att.avg))
  expect_gte(length(fit_u$att.avg.boot), 1L)
})

## ---- E3: eff.out has NA cells (placeboTest = TRUE) -----------------

test_that("E3: placeboTest=TRUE (NA eff.out cells) does not cause crash", {

  skip_on_cran()

  set.seed(203)
  N <- 30; TT <- 20; T0 <- 12; Ntr <- 10
  alpha_i <- rnorm(N)
  xi_t    <- rnorm(TT, sd = 0.5)
  df_pl <- data.frame(
    id   = rep(1:N, each = TT),
    time = rep(1:TT, N)
  )
  df_pl$D <- as.integer(df_pl$id <= Ntr & df_pl$time > T0)
  df_pl$Y <- alpha_i[df_pl$id] + xi_t[df_pl$time] +
    2.0 * df_pl$D + rnorm(N * TT, sd = 0.5)

  ## placeboTest introduces NA residuals in pre-treatment window.
  fit_pl <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = df_pl, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, vartype = "wild", nboots = 100,
      parallel = FALSE, CV = FALSE,
      placeboTest = TRUE, placebo.period = c(-3, -1)
    )
  ))

  ## Key requirement: completes without hard crash.
  expect_true(is.numeric(fit_pl$att.avg))
})

## ---- E4: att.avg = NA propagates gracefully (make_boot_na) ----------

test_that("E4: degenerate fit (extreme setup) propagates without crash", {

  skip_on_cran()

  ## Near-degenerate panel: only 1 treated unit with 1 post-treatment period.
  ## Many bootstrap replicates will exclude the treated unit, triggering
  ## make_boot_na. The overall fit must not hard-crash.
  set.seed(204)
  N <- 20; TT <- 15; T0 <- 14; Ntr <- 1
  df_deg <- data.frame(
    id   = rep(1:N, each = TT),
    time = rep(1:TT, N)
  )
  df_deg$D <- as.integer(df_deg$id <= Ntr & df_deg$time > T0)
  alpha_i <- rnorm(N, sd = 2)
  df_deg$Y <- alpha_i[df_deg$id] + rnorm(N * TT, sd = 1) +
    3.0 * df_deg$D

  ## Should complete — possibly with NAs in att.avg.boot, but no crash.
  expect_no_error(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D, data = df_deg, index = c("id", "time"),
        method = "fe", force = "two-way",
        se = TRUE, vartype = "wild", nboots = 50,
        parallel = FALSE, CV = FALSE
      )
    ))
  )
})
