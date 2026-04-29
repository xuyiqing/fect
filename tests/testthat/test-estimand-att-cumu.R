## ---------------------------------------------------------------
## Tests for estimand(fit, "att.cumu", ...).
##
## Locked invariants per ref/po-estimands-contract.md §2:
##   estimand(fit, "att.cumu", "event.time") matches effect(fit) on
##     catt / S.E. / CI.lower / CI.upper.
##   estimand(fit, "att.cumu", "overall", window = period) matches the
##     final row of att.cumu(fit, period = period) on catt / S.E. / CIs.
##
## simdata has treatment reversals so effect() / att.cumu() warn or
## error on it. We construct a small canonical no-reversal panel
## inline for the byte-equality tests.
## ---------------------------------------------------------------

.canonical_no_reversal_panel <- function(N = 30, TT = 20, seed = 1) {
  set.seed(seed)
  df <- expand.grid(id = 1:N, time = 1:TT)
  treat_start <- sample(c(NA, 8:15), N, replace = TRUE)
  df$D <- ifelse(is.na(treat_start[df$id]) | df$time < treat_start[df$id],
                 0, 1)
  df$Y <- rnorm(nrow(df)) + 0.5 * df$D + 0.1 * df$time
  df
}

.fit_no_reversal <- function(nboots = 50, keep.sims = TRUE) {
  df <- .canonical_no_reversal_panel()
  set.seed(42)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, nboots = nboots, parallel = FALSE,
      keep.sims = keep.sims
    )
  ))
}


## -- AC.1  Schema: returns documented columns -----------------------

test_that("AC.1: estimand(fit, 'att.cumu', 'event.time') has the documented schema", {

  skip_on_cran()

  fit <- .fit_no_reversal()
  est <- fect::estimand(fit, "att.cumu", "event.time")

  expect_s3_class(est, "data.frame")
  expected_cols <- c("event.time", "estimate", "se",
                     "ci.lo", "ci.hi", "n_cells", "vartype")
  expect_setequal(names(est), expected_cols)

  expect_true(all(est$event.time >= 1))
  expect_true(is.numeric(est$estimate))
})


## -- AC.2  Backcompat: matches effect(fit) numerically (LOCKED INVARIANT) --

test_that("AC.2: estimand(fit, 'att.cumu', 'event.time') matches effect(fit)", {

  skip_on_cran()

  fit <- .fit_no_reversal()
  est <- fect::estimand(fit, "att.cumu", "event.time")
  out_eff <- suppressWarnings(suppressMessages(
    fect::effect(fit, cumu = TRUE, plot = FALSE)
  ))

  ## effect() returns the augmented fit with the cumulative matrix at
  ## $effect.est.att (columns ATT / S.E. / CI.lower / CI.upper / p.value,
  ## rownames = event time).
  expect_false(is.null(out_eff))
  expect_false(is.null(out_eff$effect.est.att))
  M <- out_eff$effect.est.att
  eff_df <- as.data.frame(M)
  eff_df$event.time <- as.numeric(rownames(M))
  eff_df <- eff_df[order(eff_df$event.time), ]
  est_sorted <- est[order(est$event.time), ]

  ## Match on the post-treatment event times that effect() reports.
  shared <- intersect(eff_df$event.time, est_sorted$event.time)
  expect_true(length(shared) > 0)

  eff_sub <- eff_df[eff_df$event.time %in% shared, ]
  est_sub <- est_sorted[est_sorted$event.time %in% shared, ]

  expect_equal(est_sub$estimate, unname(eff_sub$ATT),      tolerance = 1e-8)
  expect_equal(est_sub$se,       unname(eff_sub$S.E.),     tolerance = 1e-8)
  expect_equal(est_sub$ci.lo,    unname(eff_sub$CI.lower), tolerance = 1e-8)
  expect_equal(est_sub$ci.hi,    unname(eff_sub$CI.upper), tolerance = 1e-8)
})


## -- AC.3  Backcompat: matches att.cumu(fit, period) on overall window --

test_that("AC.3: estimand(fit, 'att.cumu', 'overall', window) matches att.cumu()", {

  skip_on_cran()

  fit <- .fit_no_reversal()
  period <- c(1, 5)

  est <- fect::estimand(fit, "att.cumu", "overall", window = period)
  out_acc <- suppressWarnings(suppressMessages(
    fect::att.cumu(fit, period = period, plot = FALSE)
  ))

  ## att.cumu() returns a matrix; final row is the cumulative through
  ## period[2] starting at period[1].
  final_row <- out_acc[nrow(out_acc), ]

  expect_equal(est$estimate, unname(final_row["catt"]),     tolerance = 1e-8)
  expect_equal(est$se,       unname(final_row["S.E."]),     tolerance = 1e-8)
  expect_equal(est$ci.lo,    unname(final_row["CI.lower"]), tolerance = 1e-8)
  expect_equal(est$ci.hi,    unname(final_row["CI.upper"]), tolerance = 1e-8)
})


## -- AC.4  Reversal panels: explicit error -------------------------

test_that("AC.4: cumulative on reversal panels errors with the contract message", {

  skip_on_cran()

  ## simdata has reversal.
  data("simdata", package = "fect")
  set.seed(42)
  fit_rev <- suppressWarnings(suppressMessages(
    fect::fect(Y ~ D, data = simdata, index = c("id", "time"),
               method = "fe", force = "two-way",
               se = TRUE, nboots = 20, parallel = FALSE,
               keep.sims = TRUE)
  ))
  expect_true(isTRUE(fit_rev$hasRevs == 1))
  expect_error(
    fect::estimand(fit_rev, "att.cumu", "event.time"),
    "treatment reversals",
    fixed = FALSE
  )
})


## -- AC.5  keep.sims = FALSE: clear error path ---------------------

test_that("AC.5: cumulative without keep.sims errors with the locked wording", {

  skip_on_cran()

  fit_no_keep <- .fit_no_reversal(keep.sims = FALSE)
  ## Event-time cumulative delegates to effect(), which requires
  ## eff.boot. Without keep.sims the fit has no eff.boot; we surface
  ## the canonical fect error wording.
  expect_error(
    fect::estimand(fit_no_keep, "att.cumu", "event.time"),
    "keep\\.sims = TRUE",
    fixed = FALSE
  )
})
