## ---------------------------------------------------------------
## fect 2.4.7: follow-up fixes on PR #151 (c4331a4), issues #152-#161.
## One test_that block per item (K1-K10). Every block fails on c4331a4
## and passes after the fix. Self-contained: helpers prefixed .fu_.
## ---------------------------------------------------------------

.fu_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## One of fect's datasets, loaded into a local environment.
.fu_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "fect", envir = e)
  e[[name]]
}

## The basic interval at level 1 - a from draws b (quantile() type 7, the
## default), and whether it excludes 0.
.fu_basic_ci <- function(est, b, a) {
  q <- stats::quantile(b[!is.na(b)], c(1 - a / 2, a / 2), names = FALSE)
  c(2 * est - q[1], 2 * est - q[2])
}
.fu_excludes0 <- function(ci) ci[1] > 0 || ci[2] < 0


## -- K1  fect_mspe(): data named like a function (#152) ----------------------

test_that("K1: fect_mspe() scores a fit whose data is named df or data", {
  skip_on_cran()
  simgsynth <- .fu_data("simgsynth")
  fit_in_function <- function() {
    df <- simgsynth
    fect::fect(Y ~ D + X1 + X2, data = df, index = c("id", "time"),
               method = "ife", r = 2, CV = FALSE, se = FALSE,
               parallel = FALSE)
  }
  fit_and_score <- function() {
    df <- simgsynth
    f <- fect::fect(Y ~ D + X1 + X2, data = df, index = c("id", "time"),
                    method = "ife", r = 2, CV = FALSE, se = FALSE,
                    parallel = FALSE)
    fect::fect_mspe(f, seed = 1, k = 3)$summary$MSPE
  }
  wrap <- function(data) {
    fect::fect(Y ~ D + X1 + X2, data = data, index = c("id", "time"),
               method = "ife", r = 2, CV = FALSE, se = FALSE,
               parallel = FALSE)
  }
  inside <- .fu_quiet(fit_and_score())
  fa <- .fu_quiet(fit_in_function())
  fb <- .fu_quiet(wrap(simgsynth))
  ## scored outside the function that made the fit: the same score as inside
  expect_equal(.fu_quiet(fect::fect_mspe(fa, seed = 1, k = 3))$summary$MSPE,
               inside, tolerance = 1e-10)
  expect_equal(.fu_quiet(fect::fect_mspe(fb, seed = 1, k = 3))$summary$MSPE,
               inside, tolerance = 1e-10)
  ## a fit without a formula keeps no environment of its own: the name finds
  ## only the function stats::df, and the stop says so
  fit_no_formula <- function() {
    df <- simgsynth
    fect::fect(Y = "Y", D = "D", X = c("X1", "X2"), data = df,
               index = c("id", "time"), method = "ife", r = 2, CV = FALSE,
               se = FALSE, parallel = FALSE)
  }
  fc <- .fu_quiet(fit_no_formula())
  expect_error(.fu_quiet(fect::fect_mspe(fc, seed = 1, k = 3)),
               "`df` is not a data frame")
})


## -- K2  att.cumu(): row 1 and one-period windows (#153) --------------------

test_that("K2: att.cumu() row 1 uses the rule of its other rows and of effect(); period = c(1, 1) works", {
  skip_on_cran()
  simgsynth <- .fu_data("simgsynth")
  fit <- .fu_quiet(fect::fect(
    Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = TRUE, vartype = "bootstrap",
    nboots = 50, keep.sims = TRUE, parallel = FALSE, seed = 1))
  ac <- .fu_quiet(fect::att.cumu(fit, period = c(1, 3)))
  ef <- .fu_quiet(fect::effect(fit, period = c(1, 3)))$effect.est.att
  expect_equal(unname(ac[, c("catt", "S.E.", "CI.lower", "CI.upper", "p.value")]),
               unname(as.matrix(ef[, c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value")])),
               tolerance = 1e-8)
  ## row 1: the percentiles of the draws at event time 1
  b1 <- fit$att.boot[which(fit$time == 1), ]
  expect_equal(unname(ac[1, c("CI.lower", "CI.upper")]),
               unname(stats::quantile(b1, c(0.025, 0.975), na.rm = TRUE)),
               tolerance = 1e-10)
  ## a one-period window: that single row
  one <- .fu_quiet(fect::att.cumu(fit, period = c(1, 1)))
  expect_equal(nrow(one), 1L)
  expect_equal(unname(one[1, ]), unname(ac[1, ]), tolerance = 1e-12)
  est1 <- .fu_quiet(fect::estimand(fit, "att.cumu", "overall", window = c(1, 1)))
  expect_equal(c(est1$estimate, est1$se, est1$ci.lo, est1$ci.hi),
               unname(ac[1, c("catt", "S.E.", "CI.lower", "CI.upper")]),
               tolerance = 1e-12)
  ## no SEs: one row too
  fit0 <- .fu_quiet(fect::fect(
    Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
    method = "ife", r = 2, CV = FALSE, se = FALSE, parallel = FALSE))
  one0 <- .fu_quiet(fect::att.cumu(fit0, period = c(1, 1)))
  expect_equal(nrow(one0), 1L)
  expect_equal(unname(one0[1, 3]), unname(fit0$att[fit0$time == 1]))
})


## -- K3  jackknife + weights + placebo / carryover test (#154) --------------

test_that("K3: weighted jackknife fits run placebo and carryover tests with 7 columns", {
  skip_on_cran()
  simgsynth <- .fu_data("simgsynth")
  d <- simgsynth
  d$w <- 1 + (d$id %% 3)
  ## the columns of the unweighted jackknife fit (with covariates fect names
  ## the first one "Coef")
  ref <- .fu_quiet(fect::fect(Y ~ D + X1 + X2, data = d, index = c("id", "time"),
                              method = "fe", se = TRUE, vartype = "jackknife",
                              placeboTest = TRUE, placebo.period = c(-2, 0),
                              parallel = FALSE))
  cols7 <- colnames(ref$est.placebo)
  expect_equal(length(cols7), 7L)
  z <- stats::qnorm(0.95)
  for (warg in c("W", "W.agg")) {
    args <- list(Y ~ D + X1 + X2, data = d, index = c("id", "time"),
                 method = "fe", se = TRUE, vartype = "jackknife",
                 placeboTest = TRUE, placebo.period = c(-2, 0),
                 parallel = FALSE)
    args[[warg]] <- "w"
    f <- .fu_quiet(do.call(fect::fect, args))
    ep <- f$est.placebo
    expect_equal(colnames(ep), cols7)
    expect_equal(unname(ep[1, c("CI.lower(90%)", "CI.upper(90%)")]),
                 unname(ep[1, 1] + c(-1, 1) * z * ep[1, "S.E."]),
                 tolerance = 1e-10)
  }
  ## without covariates the fit ran before, with 5 columns
  f0 <- .fu_quiet(fect::fect(Y ~ D, data = d, index = c("id", "time"),
                             method = "fe", se = TRUE, vartype = "jackknife",
                             W = "w", placeboTest = TRUE,
                             placebo.period = c(-2, 0), parallel = FALSE))
  expect_equal(colnames(f0$est.placebo),
               c("ATT.placebo", "S.E.", "CI.lower", "CI.upper", "p.value",
                 "CI.lower(90%)", "CI.upper(90%)"))
  ## carryover test on a panel with reversals
  simdata <- .fu_data("simdata")
  s <- simdata
  s$w <- 1 + (s$id %% 3)
  fc <- .fu_quiet(fect::fect(Y ~ D + X1 + X2, data = s, index = c("id", "time"),
                             method = "fe", se = TRUE, vartype = "jackknife",
                             W = "w", carryoverTest = TRUE,
                             carryover.period = c(1, 2), parallel = FALSE))
  ec <- fc$est.carryover
  expect_equal(colnames(ec), cols7)
  expect_equal(unname(ec[1, c("CI.lower(90%)", "CI.upper(90%)")]),
               unname(ec[1, 1] + c(-1, 1) * z * ec[1, "S.E."]),
               tolerance = 1e-10)
})


## -- K4  est.group.att follows ci.method (#155) -----------------------------

test_that("K4: est.group.att shows the normal interval under 'normal' and the basic one under 'basic'", {
  skip_on_cran()
  simdata <- .fu_data("simdata")
  simdata$grp <- 1 + simdata$id %% 2
  fit_g <- function(cim) {
    .fu_quiet(fect::fect(Y ~ D + X1 + X2, data = simdata,
                         index = c("id", "time"), method = "fe",
                         group = "grp", se = TRUE, nboots = 50,
                         parallel = FALSE, seed = 1, keep.sims = TRUE,
                         ci.method = cim))
  }
  fn <- fit_g("normal")
  fb <- fit_g("basic")
  est <- fn$est.group.att[, "ATT"]
  se <- fn$est.group.att[, "S.E."]
  z <- stats::qnorm(0.975)
  expect_equal(unname(fn$est.group.att[, "CI.lower"]), unname(est - z * se),
               tolerance = 1e-10)
  expect_equal(unname(fn$est.group.att[, "CI.upper"]), unname(est + z * se),
               tolerance = 1e-10)
  expect_equal(unname(fn$est.group.att[, "p.value"]),
               unname(2 * (1 - stats::pnorm(abs(est / se)))), tolerance = 1e-10)
  ## the basic interval of each replicate's group ATTs, rebuilt from the
  ## saved replicates: the mean effect over the replicate's treated cells of
  ## the units in each group
  grp_of_unit <- tapply(simdata$grp, simdata$id, function(v) v[1])[as.character(fb$id)]
  nb <- dim(fb$eff.boot)[3]
  B <- sapply(seq_len(nb), function(b) {
    u <- as.integer(fb$colnames.boot[[b]])
    w <- length(u)
    e <- fb$eff.boot[, seq_len(w), b]
    dd <- fb$D.boot[, seq_len(w), b]
    g <- matrix(grp_of_unit[u], nrow(e), w, byrow = TRUE)
    ok <- !is.na(e) & dd == 1
    tapply(e[ok], g[ok], mean)[c("1", "2")]
  })
  expect_equal(unname(apply(B, 1, stats::sd)), unname(fb$est.group.att[, "S.E."]),
               tolerance = 1e-10)
  basic <- t(sapply(1:2, function(k) .fu_basic_ci(est[[k]], B[k, ], 0.05)))
  expect_equal(unname(fb$est.group.att[, c("CI.lower", "CI.upper")]),
               unname(basic), tolerance = 1e-10)
})


## -- K5  att.avg.unit with missing cells (#156) -----------------------------

test_that("K5: att.avg.unit averages each treated unit over its observed treated cells", {
  skip_on_cran()
  simgsynth <- .fu_data("simgsynth")
  d_a <- simgsynth[!(simgsynth$id == 101 & simgsynth$time == 5), ]
  d_b <- simgsynth[!(simgsynth$id %in% 101:105 & simgsynth$time == 5), ]
  by_hand <- function(fit) {
    tr <- which(colSums(fit$D.dat, na.rm = TRUE) > 0)
    mean(sapply(tr, function(i) {
      ok <- fit$D.dat[, i] == 1 & !is.na(fit$eff[, i])
      mean(fit$eff[ok, i])
    }))
  }
  specs <- list(
    fe     = list(method = "fe"),
    gsynth = list(method = "gsynth", r = 2, CV = FALSE),
    ife_cv = list(method = "ife", CV = TRUE, r = c(0, 2)),
    mc     = list(method = "mc", lambda = 0.1, CV = FALSE),
    cfe    = list(method = "cfe")
  )
  for (nm in names(specs)) {
    for (d in list(d_a, d_b)) {
      args <- c(list(Y ~ D + X1 + X2, data = d, index = c("id", "time"),
                     se = FALSE, parallel = FALSE), specs[[nm]])
      f <- .fu_quiet(do.call(fect::fect, args))
      expect_false(is.nan(f$att.avg.unit), label = nm)
      expect_equal(f$att.avg.unit, by_hand(f), tolerance = 1e-10, label = nm)
    }
  }
  ## the replicates too: est.avg.unit has an S.E.
  fb <- .fu_quiet(fect::fect(Y ~ D + X1 + X2, data = d_b, index = c("id", "time"),
                             method = "fe", se = TRUE, nboots = 20, seed = 1,
                             parallel = FALSE))
  expect_true(all(is.finite(fb$est.avg.unit[1, 1:4])))
})


## -- K6  normalize = TRUE: est.cm and Y.dat on the outcome's scale (#157) ----

test_that("K6: normalize = TRUE gives est.cm and Y.dat on the outcome's scale", {
  skip_on_cran()
  sim_base <- .fu_data("sim_base")
  sb <- sim_base
  sb$Ypos <- sb$Y + 30
  fit_n <- function(nrm) {
    .fu_quiet(fect::fect(Ypos ~ D + X1 + X2, data = sb, index = c("id", "time"),
                         method = "fe", force = "two-way", cm = TRUE,
                         se = FALSE, parallel = FALSE, normalize = nrm))
  }
  f0 <- fit_n(FALSE)
  f1 <- fit_n(TRUE)
  expect_equal(f1$Y.dat, f0$Y.dat, tolerance = 1e-8)
  expect_equal(f1$est.cm$fit, f0$est.cm$fit, tolerance = 1e-6)
  io0 <- fect::imputed_outcomes(f0)
  io1 <- fect::imputed_outcomes(f1)
  expect_equal(io1$Y_obs, io0$Y_obs, tolerance = 1e-8)
  expect_equal(io1$Y0_hat, io0$Y0_hat, tolerance = 1e-6)
  a0 <- .fu_quiet(fect::estimand(f0, "aptt", "event.time", vartype = "none"))
  a1 <- .fu_quiet(fect::estimand(f1, "aptt", "event.time", vartype = "none"))
  expect_equal(a1$estimate, a0$estimate, tolerance = 1e-6)
  l0 <- .fu_quiet(fect::estimand(f0, "log.att", "event.time", vartype = "none"))
  l1 <- .fu_quiet(fect::estimand(f1, "log.att", "event.time", vartype = "none"))
  expect_equal(l1$estimate, l0$estimate, tolerance = 1e-6)
  i0 <- .fu_quiet(fect::fect_iden(f0, moderator = "X1"))
  i1 <- .fu_quiet(fect::fect_iden(f1, moderator = "X1"))
  expect_equal(i1$e0$stat, i0$e0$stat, tolerance = 1e-6)
  expect_equal(i1$e1$stat, i0$e1$stat, tolerance = 1e-6)
  ## the covariate array (the fit's second `X`) and the cm plot
  expect_equal(f1[names(f1) == "X"][[2]], f0[names(f0) == "X"][[2]], tolerance = 1e-8)
  cm_y <- function(fit) {
    p <- .fu_quiet(plot(fit, type = "hte", covariate = "X1", cm = TRUE,
                        loess.fit = FALSE))
    k <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"),
                      logical(1)))[1]
    ggplot2::layer_data(p, k)$y
  }
  expect_equal(cm_y(f1), cm_y(f0), tolerance = 1e-6)
})


## -- K7  basic p-values of bootstrap fits match the interval (#158) ----------

test_that("K7: with ci.method = 'basic', p < alpha exactly when 0 is outside the basic interval", {
  skip_on_cran()
  simgsynth <- .fu_data("simgsynth")
  fit <- .fu_quiet(fect::fect(
    Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
    method = "gsynth", force = "two-way", r = 2, CV = FALSE, se = TRUE,
    vartype = "bootstrap", nboots = 200, ci.method = "basic",
    placeboTest = TRUE, placebo.period = c(-2, 0),
    parallel = FALSE, seed = 1))
  ## rows: event-time ATTs, overall ATT, coefficients, placebo ATT
  rows <- list(
    list(est = fit$est.att[, "ATT"], p = fit$est.att[, "p.value"],
         ci = fit$est.att[, c("CI.lower", "CI.upper"), drop = FALSE],
         b = fit$att.boot),
    list(est = fit$est.avg[1, "ATT.avg"], p = fit$est.avg[1, "p.value"],
         ci = fit$est.avg[, c("CI.lower", "CI.upper"), drop = FALSE],
         b = matrix(fit$att.avg.boot, 1)),
    list(est = fit$est.beta[, "Coef"], p = fit$est.beta[, "p.value"],
         ci = fit$est.beta[, c("CI.lower", "CI.upper"), drop = FALSE],
         b = fit$beta.boot),
    list(est = fit$est.placebo[1, 1], p = fit$est.placebo[1, "p.value"],
         ci = fit$est.placebo[, c("CI.lower", "CI.upper"), drop = FALSE],
         b = matrix(fit$att.placebo.boot, 1))
  )
  levels <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8)
  for (r in rows) {
    for (k in seq_along(r$est)) {
      ## the fit's interval is the basic interval at the 95% level
      expect_equal(unname(r$ci[k, ]), .fu_basic_ci(r$est[[k]], r$b[k, ], 0.05),
                   tolerance = 1e-10)
      for (a in levels) {
        ci <- .fu_basic_ci(r$est[[k]], r$b[k, ], a)
        expect_identical(unname(r$p[[k]] < a), .fu_excludes0(ci))
      }
    }
  }
  ## the dloo pre-treatment placebo rows (dloo = TRUE) follow the same rule
  fd <- .fu_quiet(fect::fect(
    Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"), method = "fe",
    se = TRUE, vartype = "bootstrap", nboots = 200, ci.method = "basic",
    dloo = TRUE, parallel = FALSE, seed = 1))
  pe <- fd$pre.est.att
  pb <- fd$pre.att.boot
  expect_equal(nrow(pb), nrow(pe))
  for (k in seq_len(nrow(pe))) {
    expect_equal(unname(pe[k, c("CI.lower", "CI.upper")]),
                 .fu_basic_ci(pe[k, "ATT"], pb[k, ], 0.05), tolerance = 1e-10)
    for (a in levels) {
      expect_identical(unname(pe[k, "p.value"] < a),
                       .fu_excludes0(.fu_basic_ci(pe[k, "ATT"], pb[k, ], a)))
    }
  }
  ## subgroup switch-off effects under ci.method = "normal": the normal p-value
  simdata <- .fu_data("simdata")
  simdata$grp <- 1 + simdata$id %% 2
  fg <- .fu_quiet(fect::fect(Y ~ D + X1 + X2, data = simdata,
                             index = c("id", "time"), method = "fe",
                             group = "grp", se = TRUE, nboots = 30,
                             parallel = FALSE, seed = 1))
  for (g in names(fg$est.group.output)) {
    off <- fg$est.group.output[[g]]$att.off
    expect_false(is.null(off))
    expect_equal(unname(off[, "p.value"]),
                 unname(2 * (1 - stats::pnorm(abs(off[, "ATT.OFF"] / off[, "S.E."])))),
                 tolerance = 1e-10)
  }
})


## -- K8  estimand() default ci.method on jackknife fits (#159) --------------

test_that("K8: estimand() on a jackknife fit defaults to ci.method = 'normal' for every type", {
  skip_on_cran()
  turnout <- .fu_data("turnout")
  fit <- .fu_quiet(fect::fect(
    turnout ~ policy_edr + policy_mail_in + policy_motor, data = turnout,
    index = c("abb", "year"), method = "fe", se = TRUE, vartype = "jackknife",
    keep.sims = TRUE, parallel = FALSE))
  calls <- list(c("att.cumu", "event.time"), c("att.cumu", "overall"),
                c("aptt", "event.time"), c("log.att", "event.time"))
  for (cl in calls) {
    d <- .fu_quiet(fect::estimand(fit, cl[1], cl[2]))
    n <- .fu_quiet(fect::estimand(fit, cl[1], cl[2], ci.method = "normal"))
    expect_identical(d, n)
  }
  ## an explicit non-normal method still stops
  expect_error(fect::estimand(fit, "att.cumu", "event.time", ci.method = "basic"),
               "not supported")
})


## -- K9  binary = TRUE stops at once (#160) --------------------------------

test_that("K9: binary = TRUE stops at once in fect() and interFE() with a plain message", {
  skip_on_cran()
  simgsynth <- .fu_data("simgsynth")
  d <- simgsynth
  d$Yb <- as.numeric(d$Y > stats::median(d$Y))
  for (frc in c("two-way", "none")) {
    expect_error(fect::fect(Yb ~ D + X1 + X2, data = d, index = c("id", "time"),
                            binary = TRUE, force = frc, r = 0, CV = FALSE,
                            se = FALSE, parallel = FALSE),
                 "binary = TRUE is not supported in this version of fect")
  }
  expect_error(fect::interFE(Yb ~ D + X1 + X2, data = d, index = c("id", "time"),
                             r = 1, force = "two-way", binary = TRUE, se = FALSE),
               "binary = TRUE is not supported in this version of fect")
})


## -- K10  estimand(): cells for att.cumu, and the stop messages (#161) ------

test_that("K10: estimand() stops on cells for att.cumu and names the real cause", {
  skip_on_cran()
  simgsynth <- .fu_data("simgsynth")
  d <- simgsynth
  d$Y <- d$Y + 20
  fit <- .fu_quiet(fect::fect(Y ~ D + X1 + X2, data = d, index = c("id", "time"),
                              method = "fe", se = TRUE, nboots = 30,
                              keep.sims = TRUE, parallel = FALSE, seed = 1))
  fit_nosims <- .fu_quiet(fect::fect(Y ~ D + X1 + X2, data = d,
                                     index = c("id", "time"), method = "fe",
                                     se = TRUE, nboots = 30, parallel = FALSE,
                                     seed = 1))
  ## att.cumu: `cells` is not silently ignored any more
  expect_error(.fu_quiet(fect::estimand(fit, "att.cumu", "overall",
                                        cells = ~ event.time <= 5)),
               "does not take `cells`")
  expect_error(.fu_quiet(fect::estimand(fit, "att.cumu", "event.time",
                                        cells = ~ event.time <= 5)),
               "does not take `cells`")
  expect_error(.fu_quiet(fect::estimand(fit, "att.cumu", "event.time",
                                        window = c(1, 5))),
               "takes `window` only with")
  ## `window` with by = "overall" still works, as att.cumu() over the window
  w <- .fu_quiet(fect::estimand(fit, "att.cumu", "overall", window = c(1, 5)))
  ac <- .fu_quiet(fect::att.cumu(fit, period = c(1, 5)))
  expect_equal(w$estimate, unname(ac[nrow(ac), "catt"]), tolerance = 1e-12)
  ## messages name the cause, with no version or "this commit" wording
  msg <- function(expr) tryCatch({ .fu_quiet(expr); NA_character_ },
                                 error = function(e) conditionMessage(e))
  m <- c(
    cohort   = msg(fect::estimand(fit, "att", "cohort")),
    caltime  = msg(fect::estimand(fit, "att", "calendar.time")),
    column   = msg(fect::estimand(fit, "att", "region")),
    cumu_coh = msg(fect::estimand(fit, "att.cumu", "cohort")),
    aptt_all = msg(fect::estimand(fit, "aptt", "overall")),
    logatt   = msg(fect::estimand(fit, "log.att", "overall")),
    aptt_c   = msg(fect::estimand(fit, "aptt", "event.time", cells = ~ event.time <= 5)),
    log_c    = msg(fect::estimand(fit, "log.att", "event.time", cells = ~ event.time <= 5)),
    att_c    = msg(fect::estimand(fit, "att", "event.time", cells = ~ event.time <= 5)),
    nosims   = msg(fect::estimand(fit_nosims, "att", "overall"))
  )
  expect_false(anyNA(m))
  expect_false(any(grepl("v2.4.0|this commit|not yet", m)))
  expect_match(m[["cohort"]], "by = \"cohort\" is not available in this release", fixed = TRUE)
  expect_match(m[["column"]], "by = \"region\" is not available in this release", fixed = TRUE)
  expect_match(m[["aptt_all"]], "Use by = \"event.time\".", fixed = TRUE)
  expect_match(m[["aptt_c"]], "does not take `cells` or `window`", fixed = TRUE)
  expect_match(m[["att_c"]], "returns fit$est.att as it is, so it does not take `cells` or `window`",
               fixed = TRUE)
  expect_match(m[["nosims"]], "needs keep.sims = TRUE in fect() for its standard error",
               fixed = TRUE)
  ## direction = "off" by event time names direction, not `by`
  simdata <- .fu_data("simdata")
  fr <- .fu_quiet(fect::fect(Y ~ D + X1 + X2, data = simdata, index = c("id", "time"),
                             method = "fe", se = TRUE, nboots = 20,
                             keep.sims = TRUE, parallel = FALSE, seed = 1))
  expect_match(msg(fect::estimand(fr, "att", "event.time", direction = "off")),
               "does not take direction = \"off\"", fixed = TRUE)
})
