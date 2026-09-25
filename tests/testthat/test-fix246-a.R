## ---------------------------------------------------------------
## fix246 section A: inference and bootstrap (items A1-A7).
##
## Regression tests for the 2.4.6 correctness run. Each block fails on
## fect dev @ 412d7ae and passes after the fix. Fits are small and serial
## (parallel = FALSE); helpers use pkg::fn for non-base functions.
## ---------------------------------------------------------------

.fx_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## turnout: 47 states x 24 years. The treated states are not the first
## columns and adopt at different times, so replicate columns (bootstrap
## order) differ from the panel's column order.
.fx_turnout <- function() {
  e <- new.env()
  utils::data("turnout", package = "fect", envir = e)
  e$turnout
}

## One cached fit per setting, so the heavier checks share fits.
.fx_turnout_fit <- local({
  cache <- list()
  function(vartype, placebo = FALSE) {
    key <- paste(vartype, placebo)
    if (!is.null(cache[[key]])) return(cache[[key]])
    args <- list(turnout ~ policy_edr + policy_mail_in + policy_motor,
                 data = .fx_turnout(), index = c("abb", "year"),
                 method = "gsynth", force = "two-way", r = 1, CV = FALSE,
                 se = TRUE, vartype = vartype, nboots = 40,
                 parallel = FALSE, keep.sims = TRUE, seed = 11)
    if (placebo) {
      args$placeboTest <- TRUE
      args$placebo.period <- c(-2, 0)
    }
    fit <- .fx_quiet(do.call(fect::fect, args))
    cache[[key]] <<- fit
    fit
  }
})

## Reference draws computed straight from the replicate slots: replicate b's
## columns are the original units fit$colnames.boot[[b]], so its event-time
## matrix is fit$T.on[, ids] and its outcomes fit$Y.dat[, ids].
.fx_rep_slices <- function(fit, b) {
  ids <- fit$colnames.boot[[b]]
  list(eff = matrix(fit$eff.boot[, seq_along(ids), b], nrow = nrow(fit$Y.dat)),
       Tev = fit$T.on[, ids, drop = FALSE],
       D   = fit$D.dat[, ids, drop = FALSE],
       Y   = fit$Y.dat[, ids, drop = FALSE])
}


## -- A1  estimand() / imputed_outcomes() read each replicate's own cells --

test_that("A1: estimand(att, overall) SE equals est.avg SE (bootstrap and parametric)", {
  skip_on_cran()
  for (vt in c("bootstrap", "parametric")) {
    fit <- .fx_turnout_fit(vt)
    est <- fect::estimand(fit, "att", "overall")
    ## the fit's own replicate ATTs: est.avg S.E. = sd(att.avg.boot)
    expect_equal(est$se, unname(fit$est.avg[1, "S.E."]), tolerance = 1e-8,
                 info = vt)
    expect_equal(est$estimate, fit$att.avg, tolerance = 1e-10, info = vt)
  }
})

test_that("A1: imputed_outcomes(replicates = TRUE) rows are each replicate's own cells", {
  skip_on_cran()
  for (vt in c("bootstrap", "parametric")) {
    fit  <- .fx_turnout_fit(vt)
    B    <- dim(fit$eff.boot)[3]
    po   <- fect::imputed_outcomes(fit)
    po_b <- fect::imputed_outcomes(fit, replicates = TRUE)
    expect_setequal(unique(po_b$replicate), seq_len(B))
    ## averaging eff within a replicate gives that replicate's att.avg
    rep_mean <- tapply(po_b$eff, po_b$replicate, mean)
    expect_equal(as.numeric(rep_mean), as.numeric(fit$att.avg.boot),
                 tolerance = 1e-10, info = vt)
    expect_equal(po_b$Y0_hat, po_b$Y_obs - po_b$eff)
    ## every replicate row points at a real treated cell of the panel
    key   <- paste(po$id, po$time)
    key_b <- paste(po_b$id, po_b$time)
    expect_true(all(key_b %in% key))
    expect_equal(po_b$Y_obs, po$Y_obs[match(key_b, key)])
  }
  ## parametric: every treated unit appears once per replicate
  fit <- .fx_turnout_fit("parametric")
  expect_equal(nrow(fect::imputed_outcomes(fit, replicates = TRUE)),
               nrow(fect::imputed_outcomes(fit)) * dim(fit$eff.boot)[3])
})

test_that("A1: placebo event-time SEs use each replicate's own cells", {
  skip_on_cran()
  fit <- .fx_turnout_fit("bootstrap", placebo = TRUE)
  pp  <- fit$placebo.period
  est <- fect::estimand(fit, "att", "event.time", test = "placebo")
  B   <- dim(fit$eff.boot)[3]
  for (k in seq_len(nrow(est))) {
    et <- est$event.time[k]
    att_b <- vapply(seq_len(B), function(b) {
      s <- .fx_rep_slices(fit, b)
      m <- !is.na(s$Tev) & s$Tev >= pp[1] & s$Tev <= pp[2] & s$Tev == et
      if (!any(m)) NA_real_ else mean(s$eff[m], na.rm = TRUE)
    }, numeric(1))
    expect_equal(est$se[k], stats::sd(att_b, na.rm = TRUE), tolerance = 1e-8,
                 info = paste("event time", et))
  }
})

test_that("A1: aptt and log.att event-time SEs use each replicate's own cells", {
  skip_on_cran()
  fit <- .fx_turnout_fit("bootstrap")
  B   <- dim(fit$eff.boot)[3]
  ap  <- fect::estimand(fit, "aptt", "event.time", ci.method = "normal")
  la  <- fect::estimand(fit, "log.att", "event.time", ci.method = "normal")
  for (k in seq_len(nrow(ap))) {
    et <- ap$event.time[k]
    draws <- vapply(seq_len(B), function(b) {
      s <- .fx_rep_slices(fit, b)
      m <- !is.na(s$D) & s$D == 1 & !is.na(s$Tev) & s$Tev == et
      if (!any(m)) return(c(NA_real_, NA_real_))
      eff <- s$eff[m]; Y <- s$Y[m]; Y0 <- Y - eff
      c(mean(eff) / mean(Y0), mean(log(Y) - log(Y0)))
    }, numeric(2))
    expect_equal(ap$se[k], stats::sd(draws[1, ], na.rm = TRUE),
                 tolerance = 1e-8, info = paste("aptt, event time", et))
    expect_equal(la$se[k], stats::sd(draws[2, ], na.rm = TRUE),
                 tolerance = 1e-8, info = paste("log.att, event time", et))
  }
})

test_that("A1: a cells formula can use variables of the calling function", {
  skip_on_cran()
  fit <- .fx_turnout_fit("bootstrap")
  by_local <- function(fit) {
    last_et_fx246 <- 3
    list(est = fect::estimand(fit, "att", "overall",
                              cells = ~ event.time <= last_et_fx246),
         po  = fect::imputed_outcomes(fit, cells = ~ event.time <= last_et_fx246))
  }
  out <- by_local(fit)
  ref <- fect::estimand(fit, "att", "overall", window = c(1, 3))
  expect_equal(out$est$estimate, ref$estimate)
  expect_equal(out$est$se, ref$se)
  expect_true(all(out$po$event.time <= 3))
})
