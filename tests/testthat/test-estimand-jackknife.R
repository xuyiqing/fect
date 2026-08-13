## ---------------------------------------------------------------
## Tests for estimand() and imputed_outcomes() behavior on jackknife fits.
##
## Covers:
##   J.1  Slot contract: jackknife eff.boot passes .validate_po_contract
##   J.2  Fast path: event.time + normal is byte-equal to fit$est.att
##   J.3  Overall path: returns finite Wald CI
##   J.4  ci.method guard: bca hard-errors with explanation
##   J.5  ci.method guard: percentile hard-errors with explanation
##   J.6  ci.method guard: basic hard-errors
##   J.7  ci.method guard: bc hard-errors
##   J.8  imputed_outcomes(replicates = TRUE) hard-errors
##   J.9  Overall with cells window: returns finite Wald CI
##   J.10 aptt event.time under jackknife
##   J.11 log.att event.time under jackknife
##   J.12 Anti-regression: bootstrap estimand still works after fix
## ---------------------------------------------------------------

## DGP-A: balanced two-way FE, true ATT = 3.0
.make_jackknife_data <- function(seed = 42, N = 40, TT = 20, tr_start = 11,
                                  tau = 3.0) {
  set.seed(seed)
  time_idx  <- rep(1:TT, each = N)
  unit_idx  <- rep(1:N, times = TT)
  is_treated <- unit_idx <= (N / 2)
  D  <- as.integer(is_treated & time_idx >= tr_start)
  alpha_i <- rnorm(N, sd = 1)
  gamma_t <- rnorm(TT, sd = 0.5)
  Y0 <- outer(gamma_t, alpha_i, "+") + rnorm(TT * N, sd = 0.5)
  Y  <- as.vector(Y0) + tau * D
  data.frame(Y = Y, D = D, id = unit_idx, time = time_idx)
}

.fit_jack <- function(simdata, keep.sims = TRUE) {
  set.seed(1)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = simdata, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, vartype = "jackknife",
      parallel = FALSE, keep.sims = keep.sims, CV = FALSE
    )
  ))
}

.fit_boot <- function(simdata, nboots = 50, keep.sims = TRUE) {
  set.seed(1)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = simdata, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, vartype = "bootstrap", nboots = nboots,
      parallel = FALSE, keep.sims = keep.sims, CV = FALSE
    )
  ))
}

simdata_A <- .make_jackknife_data()


## -- J.1  Slot contract passes for jackknife eff.boot ---------------

test_that("J.1: .validate_po_contract passes for jackknife fit", {
  skip_on_cran()
  fit <- .fit_jack(simdata_A)
  ## Slot contract is called inside estimand(); if it fails we'd get an error.
  ## Calling the fast path (no eff.boot needed) is sufficient to confirm no error.
  expect_no_error(fect::estimand(fit, "att", "event.time", ci.method = "normal"))
  ## Also confirm eff.boot dimensions match jackknife expectation.
  expect_equal(dim(fit$eff.boot)[3], ncol(fit$Y.dat))  ## third dim = N
  expect_equal(dim(fit$eff.boot)[2], ncol(fit$Y.dat) - 1L)  ## second dim = N-1
})


## -- J.2  Fast path: event.time + normal byte-equals fit$est.att ----

test_that("J.2: event.time + normal is byte-equal to fit$est.att", {
  skip_on_cran()
  fit <- .fit_jack(simdata_A)
  est <- fect::estimand(fit, "att", "event.time", ci.method = "normal")

  tol <- 1e-10
  expect_equal(nrow(est), nrow(fit$est.att))
  expect_equal(est$event.time, as.numeric(rownames(fit$est.att)))
  expect_lt(max(abs(est$estimate - unname(fit$est.att[, "ATT"]))), tol)
  expect_lt(max(abs(est$se      - unname(fit$est.att[, "S.E."]))), tol)
  expect_lt(max(abs(est$ci.lo   - unname(fit$est.att[, "CI.lower"]))), tol)
  expect_lt(max(abs(est$ci.hi   - unname(fit$est.att[, "CI.upper"]))), tol)
  expect_true(all(est$vartype == "jackknife"))
})


## -- J.3  Overall path: returns finite Wald CI ----------------------

test_that("J.3: overall + normal returns finite Wald CI", {
  skip_on_cran()
  fit <- .fit_jack(simdata_A)
  est <- fect::estimand(fit, "att", "overall", ci.method = "normal")

  expect_equal(nrow(est), 1L)
  expect_true(is.finite(est$estimate))
  expect_true(is.finite(est$se) && est$se > 0)
  expect_true(est$ci.lo < est$estimate)
  expect_true(est$estimate < est$ci.hi)
  expect_equal(est$vartype, "jackknife")

  ## Wald consistency: ci width == 2 * z * se
  expected_width <- 2 * stats::qnorm(0.975) * est$se
  expect_lt(abs((est$ci.hi - est$ci.lo) - expected_width), 1e-10)

  ## Point estimate matches mean(eff[treated])
  D_mat   <- fit$D.dat
  treated_mean <- mean(fit$eff[!is.na(D_mat) & D_mat == 1], na.rm = TRUE)
  expect_lt(abs(est$estimate - treated_mean), 1e-10)
})


## -- J.4  ci.method guard: bca hard-errors --------------------------

test_that("J.4: ci.method = 'bca' hard-errors for jackknife fit", {
  skip_on_cran()
  fit <- .fit_jack(simdata_A)
  err <- tryCatch(
    fect::estimand(fit, "att", "event.time", ci.method = "bca"),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "jackknife", ignore.case = TRUE)
  expect_match(err, "standard error|SE|sampling distribution", ignore.case = TRUE)
  expect_match(err, "normal", ignore.case = TRUE)
})


## -- J.5  ci.method guard: percentile hard-errors -------------------

test_that("J.5: ci.method = 'percentile' hard-errors for jackknife fit", {
  skip_on_cran()
  fit <- .fit_jack(simdata_A)
  err <- tryCatch(
    fect::estimand(fit, "att", "event.time", ci.method = "percentile"),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "jackknife", ignore.case = TRUE)
  expect_match(err, "normal", ignore.case = TRUE)
})


## -- J.6  ci.method guard: basic hard-errors ------------------------

test_that("J.6: ci.method = 'basic' hard-errors for jackknife fit", {
  skip_on_cran()
  fit <- .fit_jack(simdata_A)
  err <- tryCatch(
    fect::estimand(fit, "att", "event.time", ci.method = "basic"),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "jackknife", ignore.case = TRUE)
})


## -- J.7  ci.method guard: bc hard-errors ---------------------------

test_that("J.7: ci.method = 'bc' hard-errors for jackknife fit", {
  skip_on_cran()
  fit <- .fit_jack(simdata_A)
  err <- tryCatch(
    fect::estimand(fit, "att", "event.time", ci.method = "bc"),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "jackknife", ignore.case = TRUE)
})


## -- J.8  imputed_outcomes(replicates = TRUE) hard-errors -----------

test_that("J.8: imputed_outcomes(replicates = TRUE) hard-errors for jackknife", {
  skip_on_cran()
  fit <- .fit_jack(simdata_A)
  err <- tryCatch(
    fect::imputed_outcomes(fit, replicates = TRUE),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "jackknife", ignore.case = TRUE)
  expect_match(err, "column|N-1|dimension", ignore.case = TRUE)
})


## -- J.9  Overall + cells window: returns finite Wald CI ------------

test_that("J.9: overall + window returns finite Wald CI from cells-filter branch", {
  skip_on_cran()
  fit <- .fit_jack(simdata_A)
  est_full   <- fect::estimand(fit, "att", "overall", ci.method = "normal")
  est_window <- fect::estimand(fit, "att", "overall", ci.method = "normal",
                               window = c(1, 5))

  expect_equal(nrow(est_window), 1L)
  expect_true(is.finite(est_window$estimate))
  expect_true(is.finite(est_window$se) && est_window$se > 0)
  expect_true(est_window$ci.lo < est_window$estimate)
  expect_true(est_window$estimate < est_window$ci.hi)

  ## Wald consistency
  expected_width <- 2 * stats::qnorm(0.975) * est_window$se
  expect_lt(abs((est_window$ci.hi - est_window$ci.lo) - expected_width), 1e-10)

  ## Window estimate differs from full estimate (different cell set)
  expect_false(isTRUE(all.equal(est_window$estimate, est_full$estimate)))
})


## -- J.10 aptt event.time under jackknife ---------------------------

test_that("J.10: aptt event.time returns finite SE/CI under jackknife", {
  skip_on_cran()
  ## DGP-B: Y > 0 everywhere for aptt
  set.seed(42)
  simdata_B <- simdata_A
  simdata_B$Y <- abs(simdata_A$Y) + 5
  fit <- .fit_jack(simdata_B)

  est <- fect::estimand(fit, "aptt", "event.time", ci.method = "normal")
  expect_s3_class(est, "data.frame")
  expect_true(all(est$vartype == "jackknife"))

  post_rows <- est[!is.na(est$se) & is.finite(est$se), ]
  expect_gt(nrow(post_rows), 0L)
  ## All finite-SE rows have valid CIs
  expect_true(all(post_rows$ci.lo < post_rows$estimate))
  expect_true(all(post_rows$estimate < post_rows$ci.hi))
  ## Wald consistency for all finite-SE rows
  widths <- post_rows$ci.hi - post_rows$ci.lo
  expected <- 2 * stats::qnorm(0.975) * post_rows$se
  expect_lt(max(abs(widths - expected)), 1e-10)
})


## -- J.11 log.att event.time under jackknife ------------------------

test_that("J.11: log.att event.time returns finite SE/CI under jackknife", {
  skip_on_cran()
  set.seed(42)
  simdata_B <- simdata_A
  simdata_B$Y <- abs(simdata_A$Y) + 5
  fit <- .fit_jack(simdata_B)

  est <- fect::estimand(fit, "log.att", "event.time", ci.method = "normal")
  expect_s3_class(est, "data.frame")
  expect_true(all(est$vartype == "jackknife"))

  post_rows <- est[!is.na(est$se) & is.finite(est$se), ]
  expect_gt(nrow(post_rows), 0L)
  expect_true(all(post_rows$ci.lo < post_rows$estimate))
  expect_true(all(post_rows$estimate < post_rows$ci.hi))
  ## Wald consistency
  widths <- post_rows$ci.hi - post_rows$ci.lo
  expected <- 2 * stats::qnorm(0.975) * post_rows$se
  expect_lt(max(abs(widths - expected)), 1e-10)
})


## -- J.12 Anti-regression: bootstrap estimand still works ----------

test_that("J.12: bootstrap estimand (att, event.time) unaffected by jackknife fix", {
  skip_on_cran()
  fit <- .fit_boot(simdata_A)
  est <- fect::estimand(fit, "att", "event.time")

  tol <- 1e-10
  expect_equal(nrow(est), nrow(fit$est.att))
  expect_lt(max(abs(est$estimate - unname(fit$est.att[, "ATT"]))), tol)
  expect_true(all(est$vartype == "bootstrap"))

  ## overall also works
  est_ov <- fect::estimand(fit, "att", "overall")
  expect_true(is.finite(est_ov$estimate))
  expect_true(is.finite(est_ov$se))
})




## -- S-11 extended  Anti-regression: bootstrap path unchanged -------------------
## Tests the three calls from test-spec.md §13:
##   est_boot_et  (att, event.time, default ci.method)
##   est_boot_ov  (att, overall, default ci.method)
##   est_boot_pct (att, overall, ci.method = "percentile")
##
## Note: estimand(boot, "att", "event.time", ci.method="percentile") currently
## raises "not yet implemented" for by="event.time" + non-normal ci.method on
## a non-placebo/carryover test (pre-existing limitation, not a regression
## introduced by the jackknife fix). We use by="overall" + "percentile" instead,
## which IS supported and exercises the same downstream .compute_ci() code path.

test_that("S-11: bootstrap att/event.time byte-equals fit$est.att (anti-regression)", {
  skip_on_cran()
  fit <- .fit_boot(simdata_A, nboots = 200)
  est_boot_et <- fect::estimand(fit, "att", "event.time")
  expect_s3_class(est_boot_et, "data.frame")
  tol <- 1e-10
  expect_lt(max(abs(est_boot_et$estimate - unname(fit$est.att[, "ATT"]))), tol)
  expect_true(all(est_boot_et$vartype == "bootstrap"))
})

test_that("S-11: bootstrap att/overall is finite (anti-regression)", {
  skip_on_cran()
  fit <- .fit_boot(simdata_A, nboots = 200)
  est_boot_ov <- fect::estimand(fit, "att", "overall")
  expect_true(is.finite(est_boot_ov$estimate))
  expect_true(is.finite(est_boot_ov$se))
  expect_equal(est_boot_ov$vartype, "bootstrap")
})

test_that("S-11: bootstrap att/overall percentile returns valid CI (anti-regression)", {
  skip_on_cran()
  fit <- .fit_boot(simdata_A, nboots = 200)
  est_boot_pct <- fect::estimand(fit, "att", "overall",
                                  ci.method = "percentile")
  expect_s3_class(est_boot_pct, "data.frame")
  expect_equal(est_boot_pct$vartype, "bootstrap")
  expect_true(is.finite(est_boot_pct$ci.lo))
  expect_true(is.finite(est_boot_pct$ci.hi))
})


## -- S-12  Anti-regression: parametric fit unchanged -------------------------
## Parametric bootstrap requires never-treated control units, so we use a DGP
## where only 12 of 40 units are treated (units 1..12, periods 13..20), leaving
## 28 never-treated controls. This matches the DGP convention from
## test-estimand-parametric-cifix.R (make_panel_A).

.make_param_data <- function(seed = 42) {
  set.seed(seed)
  N <- 40L; TT <- 20L; T0 <- 12L; Ntr <- 12L
  alpha_i <- rnorm(N, 0, 2)
  xi_t    <- rnorm(TT, 0, 1)
  D       <- matrix(0L, TT, N)
  D[(T0 + 1L):TT, 1L:Ntr] <- 1L
  eps     <- matrix(rnorm(N * TT, 0, 1), TT, N)
  Y       <- outer(xi_t, rep(1, N)) + outer(rep(1, TT), alpha_i) + 3.0 * D + eps
  data.frame(id   = rep(1:N, each = TT),
             time = rep(1:TT, N),
             Y    = as.vector(Y),
             D    = as.vector(D))
}

.fit_param_s12 <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    skip_on_cran()
    d <- .make_param_data()
    set.seed(42)
    cached <<- suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D, data = d, index = c("id", "time"),
        method = "fe", force = "two-way",
        se = TRUE, vartype = "parametric", nboots = 200,
        time.component.from = "nevertreated",
        parallel = FALSE, keep.sims = TRUE, CV = FALSE
      )
    ))
    cached
  }
})

test_that("S-12: parametric fit estimand(att, event.time) returns 'parametric' vartype", {
  skip_on_cran()
  fit <- .fit_param_s12()
  est <- fect::estimand(fit, "att", "event.time")
  expect_s3_class(est, "data.frame")
  expect_true(all(est$vartype == "parametric"))
})

test_that("S-12: parametric fit estimand(att, overall) is finite and 'parametric'", {
  skip_on_cran()
  fit <- .fit_param_s12()
  est <- fect::estimand(fit, "att", "overall")
  expect_s3_class(est, "data.frame")
  expect_true(is.finite(est$estimate))
  expect_true(is.finite(est$se))
  expect_equal(est$vartype, "parametric")
})


## -- S-SIM  100-rep coverage simulation: jackknife + normal, DGP-A -----------
##
## Acceptance criterion: coverage_jack >= 0.85 (true ATT = 3.0).
## Uses fit$est.avg to avoid any slot-contract overhead inside the loop;
## after the fix this equals estimand(fit, "att", "overall")$ci.lo/ci.hi.

.run_coverage_jack <- function(seed, N = 40, TT = 20, tau = 3.0,
                                tr_start = 11) {
  set.seed(seed)
  n_treated  <- N / 2
  unit_idx   <- rep(1:N, times = TT)
  time_idx   <- rep(1:TT, each = N)
  is_treated <- unit_idx <= n_treated
  D          <- as.integer(is_treated & time_idx >= tr_start)
  alpha_i    <- rnorm(N, sd = 1)
  gamma_t    <- rnorm(TT, sd = 0.5)
  Y0  <- outer(gamma_t, alpha_i, "+")[cbind(time_idx, unit_idx)] +
          rnorm(N * TT, sd = 0.5)
  Y   <- Y0 + tau * D
  dat <- data.frame(Y = Y, D = D, id = unit_idx, time = time_idx)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(Y ~ D, data = dat, index = c("id", "time"),
               method = "fe", force = "two-way",
               se = TRUE, vartype = "jackknife",
               parallel = FALSE, keep.sims = FALSE, CV = FALSE)
  ))
  att_avg <- fit$est.avg[1, "ATT.avg"]
  ci_lo   <- fit$est.avg[1, "CI.lower"]
  ci_hi   <- fit$est.avg[1, "CI.upper"]
  list(att = att_avg, covered = tau >= ci_lo & tau <= ci_hi)
}

test_that("S-SIM: jackknife + normal CI coverage >= 0.85 on DGP-A (100 reps)", {
  skip_on_cran()
  R             <- 100L
  results       <- lapply(seq_len(R), .run_coverage_jack)
  coverage_jack <- mean(sapply(results, `[[`, "covered"))
  expect_gte(coverage_jack, 0.85,
             label = paste0("jackknife coverage = ", round(coverage_jack, 3),
                            " (must be >= 0.85)"))
})
