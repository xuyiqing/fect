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


## -- A2  effect() / att.cumu() use the fit's stored vartype ------------------

.fx_simgsynth <- function() {
  e <- new.env()
  utils::data("simgsynth", package = "fect", envir = e)
  e$simgsynth
}

## Parametric gsynth fit with vartype passed as a variable, so fit$call$vartype
## is a symbol (as with any wrapper that forwards its own argument).
.fx_param_fit <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    vt <- "parametric"
    cached <<- .fx_quiet(fect::fect(
      Y ~ D + X1 + X2, data = .fx_simgsynth(), index = c("id", "time"),
      method = "gsynth", force = "two-way", r = 2, CV = FALSE, se = TRUE,
      vartype = vt, nboots = 50, parallel = FALSE, keep.sims = TRUE, seed = 3))
    cached
  }
})

.fx_covers <- function(tab, est_col = "ATT") {
  ok <- !is.na(tab[, est_col])
  all(tab[ok, "CI.lower"] <= tab[ok, est_col] & tab[ok, est_col] <= tab[ok, "CI.upper"])
}

test_that("A2: effect() CIs are centred on the estimate for parametric fits", {
  skip_on_cran()
  fit <- .fx_param_fit()
  expect_false(is.character(fit$call$vartype))   # vartype passed as a variable
  expect_identical(fit$vartype, "parametric")
  withr::local_options(fect.suppress_estimand_deprecation = TRUE)
  ef <- fect::effect(fit)$effect.est.att
  expect_true(.fx_covers(ef))
  ## the same when the call carries no vartype at all (gsynth's default)
  fit2 <- fit
  fit2$call$vartype <- NULL
  ef2 <- fect::effect(fit2)$effect.est.att
  expect_equal(ef2, ef)
  ## the stored vartype wins over a different literal in the call
  fit3 <- fit
  fit3$call$vartype <- "bootstrap"
  expect_equal(fect::effect(fit3)$effect.est.att, ef)
})

test_that("A2: effect() works for one selected unit and for a one-treated-unit fit", {
  skip_on_cran()
  withr::local_options(fect.suppress_estimand_deprecation = TRUE)
  fit <- .fx_param_fit()
  tr_ids <- colnames(fit$eff)[colSums(fit$D.dat) > 0]
  one <- fect::effect(fit, id = tr_ids[1])$effect.est.att
  expect_true(is.matrix(one))
  expect_true(all(is.finite(one[, "S.E."])))
  expect_true(.fx_covers(one))

  d <- .fx_simgsynth()
  tr_all <- unique(d$id[d$D == 1])
  d1 <- d[!(d$id %in% tr_all[-1]), ]           # one treated unit + controls
  fit1 <- .fx_quiet(fect::fect(
    Y ~ D + X1 + X2, data = d1, index = c("id", "time"), method = "gsynth",
    force = "two-way", r = 2, CV = FALSE, se = TRUE, vartype = "parametric",
    nboots = 30, parallel = FALSE, keep.sims = TRUE, seed = 3))
  ef1 <- fect::effect(fit1)$effect.est.att
  expect_true(all(is.finite(ef1[, "S.E."])))
})

test_that("A2: att.cumu() and estimand('att.cumu') use normal CIs for parametric draws", {
  skip_on_cran()
  withr::local_options(fect.suppress_estimand_deprecation = TRUE)
  fit <- .fx_param_fit()
  ac <- fect::att.cumu(fit, period = c(1, 5), plot = FALSE)
  z <- stats::qnorm(0.975)
  for (i in 2:nrow(ac)) {
    expect_equal(unname(ac[i, "CI.lower"]), unname(ac[i, "catt"] - z * ac[i, "S.E."]))
    expect_equal(unname(ac[i, "CI.upper"]), unname(ac[i, "catt"] + z * ac[i, "S.E."]))
  }
  es <- fect::estimand(fit, "att.cumu", "overall", window = c(1, 5))
  expect_true(es$ci.lo <= es$estimate && es$estimate <= es$ci.hi)
  expect_equal(es$estimate, unname(ac[nrow(ac), "catt"]))
})

test_that("A2: every se = TRUE fit stores its vartype; se = FALSE fits do not", {
  skip_on_cran()
  d <- .fx_simgsynth()
  f0 <- .fx_quiet(fect::fect(Y ~ D, data = d, index = c("id", "time"),
                             method = "fe", se = FALSE, parallel = FALSE))
  expect_null(f0$vartype)
  expect_identical(.fx_param_fit()$vartype, "parametric")
})


## -- A3  normalize = TRUE leaves parametric SEs on the outcome's scale -------

test_that("A3: parametric SEs are the same with and without normalize", {
  skip_on_cran()
  fits <- lapply(c(FALSE, TRUE), function(nz) {
    .fx_quiet(fect::fect(
      turnout ~ policy_edr + policy_mail_in + policy_motor,
      data = .fx_turnout(), index = c("abb", "year"), method = "gsynth",
      force = "two-way", r = 1, CV = FALSE, se = TRUE,
      vartype = "parametric", nboots = 30, parallel = FALSE, seed = 11,
      normalize = nz))
  })
  expect_equal(fits[[2]]$att.avg, fits[[1]]$att.avg, tolerance = 1e-6)
  ## normalization only rescales the problem; the SEs must not move
  ## (the EM tolerance leaves differences far below 1e-4)
  expect_equal(fits[[2]]$est.avg[1, "S.E."], fits[[1]]$est.avg[1, "S.E."],
               tolerance = 1e-4)
  expect_equal(fits[[2]]$est.att[, "S.E."], fits[[1]]$est.att[, "S.E."],
               tolerance = 1e-4)
})


## -- A4  weights with the parametric bootstrap ------------------------------

.fx_param_w_fit <- function(d, ...) {
  .fx_quiet(fect::fect(
    Y ~ D + X1 + X2, data = d, index = c("id", "time"), method = "gsynth",
    force = "two-way", CV = FALSE, r = 2, se = TRUE, vartype = "parametric",
    nboots = 20, parallel = FALSE, seed = 5, ...))
}

test_that("A4: weighted parametric bootstrap runs and each replicate uses the weights", {
  skip_on_cran()
  d <- .fx_simgsynth()
  d$w2 <- 0.5 + (as.integer(factor(d$id)) %% 7) / 4
  fit <- .fx_param_w_fit(d, W = "w2", keep.sims = TRUE)
  expect_true(is.finite(fit$est.avg[1, "S.E."]))
  expect_equal(unname(fit$est.avg[1, "S.E."]), stats::sd(fit$att.avg.boot))
  ## replicate b's ATT is the W-weighted mean of its treated cells, with the
  ## weights of its own units (colnames.boot[[b]])
  Wm <- matrix(NA_real_, nrow(fit$Y.dat), ncol(fit$Y.dat))
  Wm[cbind(match(d$time, fit$rawtime), match(d$id, fit$id))] <- d$w2
  att_w <- vapply(seq_len(dim(fit$eff.boot)[3]), function(b) {
    ids <- fit$colnames.boot[[b]]
    eb  <- fit$eff.boot[, seq_along(ids), b]
    trt <- fit$D.dat[, ids] == 1
    sum(eb * Wm[, ids] * trt) / sum(Wm[, ids] * trt)
  }, numeric(1))
  expect_equal(att_w, as.numeric(fit$att.avg.boot), tolerance = 1e-10)
  ## W.agg alone (weights in the aggregation only) runs as well
  fa <- .fx_param_w_fit(d, W.agg = "w2")
  expect_true(is.finite(fa$est.avg[1, "S.E."]))
})

test_that("A4: weights equal to 1 reproduce the unweighted parametric SEs", {
  skip_on_cran()
  d <- .fx_simgsynth()
  d$w1 <- 1
  f0 <- .fx_param_w_fit(d)
  f1 <- .fx_param_w_fit(d, W = "w1")
  expect_equal(f1$est.avg, f0$est.avg)
  expect_equal(f1$est.att, f0$est.att)
})


## -- A5  resampling from length-1 vectors; dropped replicates are counted ----

## simgsynth with only its first treated unit, relabelled so that it is the
## LAST column of the panel (sample(k, ...) on a length-1 vector draws 1:k).
.fx_one_treated <- function() {
  d  <- .fx_simgsynth()
  tr <- sort(unique(d$id[d$D == 1]))
  d  <- d[!(d$id %in% tr[-1]), ]
  d$id[d$id == tr[1]] <- 999
  d
}

.fx_messages <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(
    suppressWarnings(expr),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  list(value = val, messages = msgs)
}

test_that("A5: case bootstrap with one treated unit keeps every draw", {
  skip_on_cran()
  run <- .fx_messages(fect::fect(
    Y ~ D, data = .fx_one_treated(), index = c("id", "time"), method = "fe",
    force = "two-way", se = TRUE, vartype = "bootstrap", nboots = 30,
    parallel = FALSE, seed = 7))
  fit <- run$value
  expect_equal(length(fit$att.avg.boot), 30L)
  expect_true(is.finite(fit$est.avg[1, "S.E."]))
  expect_false(any(grepl("replicates failed", run$messages)))
})

test_that("A5: dropped replicates are reported with their count", {
  skip_on_cran()
  ## jackknife with one treated unit: the replicate that leaves it out has
  ## no treated unit and is dropped
  d <- .fx_one_treated()
  run <- .fx_messages(fect::fect(
    Y ~ D, data = d, index = c("id", "time"), method = "fe",
    force = "two-way", se = TRUE, vartype = "jackknife", parallel = FALSE))
  N <- length(unique(d$id))
  expect_true(any(grepl(sprintf(
    "1 of %d jackknife replicates failed and were dropped; uncertainty estimates use the remaining %d",
    N, N - 1L), run$messages, fixed = TRUE)))
  expect_equal(length(run$value$att.avg.boot), N - 1L)
})

test_that("A5: parametric bootstrap stops clearly with fewer than two never-treated units", {
  skip_on_cran()
  d  <- .fx_simgsynth()
  tr <- unique(d$id[d$D == 1])
  d2 <- d[d$id %in% c(tr, max(d$id)), ]           # one never-treated unit
  expect_error(
    .fx_quiet(fect::fect(
      Y ~ D, data = d2, index = c("id", "time"), method = "gsynth",
      force = "two-way", r = 0, CV = FALSE, se = TRUE,
      vartype = "parametric", nboots = 20, parallel = FALSE, seed = 7)),
    "needs at least two never-treated units (found 1)", fixed = TRUE)
})


## -- A6  cluster bootstrap: varying replicate width, cl checks ---------------

## turnout with clusters of unequal size (first letter of the state code).
.fx_turnout_cl <- function() {
  d <- .fx_turnout()
  d$region <- substr(d$abb, 1, 1)
  d
}

.fx_cl_fit <- function(d, ...) {
  fect::fect(turnout ~ policy_edr + policy_mail_in + policy_motor,
             data = d, index = c("abb", "year"), method = "gsynth",
             force = "two-way", r = 1, CV = FALSE, se = TRUE, nboots = 30,
             parallel = FALSE, seed = 11, ...)
}

.fx_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(
    suppressMessages(expr),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  list(value = val, warnings = w)
}

test_that("A6: unequal clusters run with keep.sims = TRUE; effect()/estimand() work", {
  skip_on_cran()
  d <- .fx_turnout_cl()
  run <- .fx_warnings(.fx_cl_fit(d, cl = "region", keep.sims = TRUE))
  fit <- run$value
  N <- ncol(fit$Y.dat)
  widths <- lengths(fit$colnames.boot)
  expect_equal(length(widths), dim(fit$eff.boot)[3])
  expect_true(length(unique(widths)) > 1)          # replicate widths vary
  expect_equal(dim(fit$eff.boot)[2], max(N, widths))
  ## columns beyond a replicate's own width are NA padding
  b <- which.min(widths)
  expect_true(all(is.na(fit$eff.boot[, -seq_len(widths[b]), b])))
  ## same SEs as the fit without keep.sims
  fit0 <- .fx_quiet(.fx_cl_fit(d, cl = "region"))
  expect_equal(fit$est.avg, fit0$est.avg)
  ## post-hoc readers use each replicate's own width and units
  est <- fect::estimand(fit, "att", "overall")
  expect_equal(est$se, unname(fit$est.avg[1, "S.E."]), tolerance = 1e-8)
  withr::local_options(fect.suppress_estimand_deprecation = TRUE)
  ef <- fect::effect(fit)$effect.est.att
  expect_true(all(is.finite(ef[, "S.E."])))
  ## no "Skipping replicate" warnings from the average-outcome block
  expect_false(any(grepl("Skipping replicate", run$warnings)))
})

test_that("A6: counterfactual bands use every cluster-bootstrap replicate", {
  skip_on_cran()
  run <- .fx_warnings(.fx_cl_fit(.fx_turnout_cl(), cl = "region"))
  expect_false(any(grepl("Skipping replicate", run$warnings)))
  ya <- run$value$Y.avg
  expect_true(sum(is.finite(ya$lower.ct)) > 0)
})

test_that("A6: cl is validated (one name, a column, constant within units, two clusters)", {
  skip_on_cran()
  d <- .fx_turnout_cl()
  expect_error(.fx_cl_fit(d, cl = c("region", "abb")),
               "\"cl\" must be a single column name.", fixed = TRUE)
  expect_error(.fx_cl_fit(d, cl = "no_such_col"),
               "is not a column of data", fixed = TRUE)
  d$clv <- d$region
  d$clv[d$year == 2000] <- "Z"                      # changes within units
  expect_error(.fx_cl_fit(d, cl = "clv"), "must be constant within each unit")
  d$one <- "A"
  expect_error(.fx_cl_fit(d, cl = "one"), "needs at least two clusters")
})

test_that("A6: cl under the parametric bootstrap warns and is not printed as clustered", {
  skip_on_cran()
  run <- .fx_warnings(.fx_cl_fit(.fx_turnout_cl(), cl = "region",
                                 vartype = "parametric"))
  expect_true(any(grepl("vartype = \"parametric\" with cl = ...: the cl argument is ignored",
                        run$warnings, fixed = TRUE)))
  out <- utils::capture.output(print(run$value))
  expect_false(any(grepl("^Cluster SE", out)))
})


## -- A7  degenerate replicates are dropped and counted, never crash ----------

## The small-control-pool panel of the issue sweep (work/batchA/diag032c.R):
## 3 treated units, Nco never-treated units, 30 periods, 2 factors.
.fx_small_pool <- function(Nco, Ntr = 3, TT = 30, T0 = 20, seed = 32) {
  set.seed(seed)
  N <- Ntr + Nco
  d <- expand.grid(time = 1:TT, id = 1:N)[, c("id", "time")]
  F <- matrix(stats::rnorm(TT * 2), TT, 2)
  L <- matrix(stats::rnorm(N * 2), N, 2)
  d$X1 <- stats::rnorm(nrow(d))
  d$X2 <- stats::rnorm(nrow(d))
  d$D <- as.numeric(d$id <= Ntr & d$time > T0)
  d$Y <- 1 + stats::rnorm(N)[d$id] + stats::rnorm(TT)[d$time] +
    rowSums(F[d$time, ] * L[d$id, ]) + 0.5 * d$X1 + 0.3 * d$X2 + 2 * d$D +
    stats::rnorm(nrow(d), sd = 0.3)
  d
}

.fx_drop_count <- function(msgs, nboots) {
  pat <- sprintf("^([0-9]+) of %d (bootstrap|jackknife) replicates failed and were dropped", nboots)
  hit <- grep(pat, msgs, value = TRUE)
  if (length(hit) == 0) return(0L)
  as.integer(sub(paste0(pat, ".*$"), "\\1", hit[1]))
}

test_that("A7: runs with four never-treated units finish; failed replicates are counted", {
  skip_on_cran()
  d <- .fx_small_pool(Nco = 4)
  for (vt in c("bootstrap", "parametric")) {
    for (ks in c(FALSE, TRUE)) {
      run <- .fx_messages(fect::fect(
        Y ~ D, data = d, index = c("id", "time"), method = "gsynth",
        force = "two-way", r = 2, CV = FALSE, se = TRUE, vartype = vt,
        nboots = 50, parallel = FALSE, seed = 1, keep.sims = ks))
      fit <- run$value
      k <- .fx_drop_count(run$messages, 50)
      expect_gt(k, 0)
      expect_equal(length(fit$att.avg.boot), 50L - k)
      expect_true(is.finite(fit$est.avg[1, "S.E."]))
      if (ks) expect_equal(dim(fit$eff.boot)[3], 50L - k)
    }
  }
})

test_that("A7: collinear factors in the main fit stop with a plain message", {
  skip_on_cran()
  ## treated units with 2 pre-treatment periods cannot identify 2 loadings
  ## plus the unit intercept
  d <- .fx_small_pool(Nco = 10, T0 = 2)
  for (m in c("gsynth", "cfe")) {
    expect_error(
      .fx_quiet(fect::fect(
        Y ~ D, data = d, index = c("id", "time"), method = m,
        force = "two-way", r = 2, CV = FALSE, se = FALSE, min.T0 = 2,
        time.component.from = "nevertreated", parallel = FALSE)),
      "factor loadings cannot be estimated", fixed = TRUE, info = m)
  }
})

## Replicate refits that return (instead of failing) with the 3-field list the
## pre-fix early returns gave: every 5th bootstrap refit.
.fx_degenerate_refits <- function() {
  real_fn <- fect:::fect_nevertreated
  calls <- 0L
  function(..., boot = 0) {
    if (isTRUE(boot == 1)) {
      calls <<- calls + 1L
      if (calls %% 5L == 0L) {
        return(list(att = rep(NA, 30), att.avg = NA, beta = matrix(NA, 0, 1)))
      }
    }
    real_fn(..., boot = boot)
  }
}

.fx_boot_simgsynth <- function(keep.sims) {
  .fx_messages(fect::fect(
    Y ~ D, data = .fx_simgsynth(), index = c("id", "time"),
    method = "gsynth", force = "two-way", r = 2, CV = FALSE, se = TRUE,
    vartype = "bootstrap", nboots = 20, parallel = FALSE, seed = 2,
    keep.sims = keep.sims))
}

test_that("A7: a refit that returns a degenerate result counts as a failed replicate", {
  skip_on_cran()
  testthat::local_mocked_bindings(fect_nevertreated = .fx_degenerate_refits(),
                                  .package = "fect")
  for (ks in c(FALSE, TRUE)) {
    run <- .fx_boot_simgsynth(keep.sims = ks)
    expect_equal(.fx_drop_count(run$messages, 20), 4L)
    expect_equal(length(run$value$att.avg.boot), 16L)
  }
})

test_that("A7: the collectors drop a replicate whose fields do not fit, instead of aborting", {
  skip_on_cran()
  ## switch the shape check off, so the degenerate refits reach the collectors
  testthat::local_mocked_bindings(
    fect_nevertreated = .fx_degenerate_refits(),
    .fect_boot_result_ok = function(...) TRUE,
    .package = "fect")
  for (ks in c(FALSE, TRUE)) {
    run <- .fx_boot_simgsynth(keep.sims = ks)
    expect_equal(.fx_drop_count(run$messages, 20), 4L)
    expect_equal(length(run$value$att.avg.boot), 16L)
    expect_true(is.finite(run$value$est.avg[1, "S.E."]))
  }
})
