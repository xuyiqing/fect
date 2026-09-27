## ---------------------------------------------------------------
## fect 2.4.7, PR #151 review note N4: tests for the weighted paths of
## estimand() that W1 (db1076a) left without a committed test. They are
## the jackknife overall ATT over a window (the branch that recomputes
## each leave-one-out replicate), the jackknife APTT and log ATT, and the
## BCa intervals of a bootstrap fit, whose acceleration uses weighted
## leave-one-out values.
##
## Every reference is computed in this file from the definitions, from
## the fit's stored pieces (id, eff, Y.dat, D.dat, T.on, eff.boot,
## colnames.boot) and the weight column of the data. No fect function
## computes a reference. Every block fails on b1dded6, where these paths
## ignored the weights. Helpers prefixed .w1_.
##
## Data: turnout (47 states, 24 elections; 9 states adopt election-day
## registration at different times; the outcome is positive, as log ATT
## needs), with state weights w = 1, 2 or 3.
## ---------------------------------------------------------------

## turnout with a state-level weight column w.
.w1_data <- function() {
  e <- new.env()
  utils::data("turnout", package = "fect", envir = e)
  d <- e$turnout
  d$w <- 1 + (as.integer(factor(d$abb)) %% 3)
  d
}

## The two weighted fits, made once per file: gsynth without factors,
## aggregation weights W = "w". The jackknife has one replicate per
## state (47); the bootstrap 50 replicates with a fixed seed. Only a
## message about the F statistic is muffled.
.w1_fit <- local({
  cached <- list()
  function(vartype) {
    if (!is.null(cached[[vartype]])) return(cached[[vartype]])
    skip_on_cran()
    d <- .w1_data()
    fit <- if (vartype == "jackknife") {
      suppressMessages(fect::fect(
        turnout ~ policy_edr + policy_mail_in + policy_motor, data = d,
        index = c("abb", "year"), method = "gsynth", r = 0, CV = FALSE,
        force = "two-way", min.T0 = 5, se = TRUE, vartype = "jackknife",
        keep.sims = TRUE, parallel = FALSE, W = "w"))
    } else {
      suppressMessages(fect::fect(
        turnout ~ policy_edr + policy_mail_in + policy_motor, data = d,
        index = c("abb", "year"), method = "gsynth", r = 0, CV = FALSE,
        force = "two-way", min.T0 = 5, se = TRUE, vartype = "bootstrap",
        nboots = 50, seed = 2139, keep.sims = TRUE, parallel = FALSE,
        W = "w"))
    }
    cached[[vartype]] <<- fit
    fit
  }
})

## The data's weights as a matrix on the fit's panel (time by unit, like
## fit$Y.dat): column i holds the weight of state fit$id[i].
.w1_weights <- function(fit) {
  d <- .w1_data()
  w_state <- tapply(d$w, d$abb, function(v) v[1])
  matrix(unname(w_state[as.character(fit$id)]), nrow = nrow(fit$Y.dat),
         ncol = ncol(fit$Y.dat), byrow = TRUE)
}

## Treated cells after onset: D = 1 with an event time.
.w1_treated <- function(fit) {
  !is.na(fit$D.dat) & fit$D.dat == 1 & !is.na(fit$T.on)
}

## The three estimands over a set of cells, from the cells' effects, their
## observed outcomes y (so Y(0) = y - eff) and their weights w.
.w1_wmean <- function(v, w) sum(w * v) / sum(w)
.w1_att   <- function(eff, y, w) .w1_wmean(eff, w)
.w1_aptt  <- function(eff, y, w) .w1_wmean(eff, w) / .w1_wmean(y - eff, w)
.w1_latt  <- function(eff, y, w) .w1_wmean(log(y) - log(y - eff), w)

## Jackknife replicates. Replicate j is the fit without unit j: the
## columns of eff.boot[, , j] are the other units in panel order. Its
## value of the estimand uses those units' cells in `mask`, their observed
## outcomes and their weights, that is, the weights without unit j.
.w1_jack_reps <- function(fit, W, mask, stat) {
  vapply(seq_len(dim(fit$eff.boot)[3]), function(j) {
    m <- mask[, -j, drop = FALSE]
    stat(fit$eff.boot[, , j][m], fit$Y.dat[, -j, drop = FALSE][m],
         W[, -j, drop = FALSE][m])
  }, numeric(1))
}

## Tukey's jackknife SE from the n leave-one-out estimates theta:
## sqrt((n - 1) / n * sum((theta - mean(theta))^2)) (Efron and Tibshirani
## 1993, ch. 11).
.w1_jack_se <- function(theta) {
  n <- length(theta)
  sqrt((n - 1) / n * sum((theta - mean(theta))^2))
}

## Bootstrap replicates. Replicate b is the fit on the units
## colnames.boot[[b]], drawn with replacement; they are the columns of
## eff.boot[, , b]. Each drawn copy of a unit brings its cells in `mask`,
## its observed outcomes and its weight.
.w1_boot_reps <- function(fit, W, mask, stat) {
  vapply(seq_len(dim(fit$eff.boot)[3]), function(b) {
    u  <- as.integer(fit$colnames.boot[[b]])
    m  <- mask[, u, drop = FALSE]
    eb <- matrix(fit$eff.boot[, seq_along(u), b], nrow = nrow(mask))
    stat(eb[m], fit$Y.dat[, u, drop = FALSE][m], W[, u, drop = FALSE][m])
  }, numeric(1))
}

## Leave-one-out values for the BCa acceleration: the estimand over the
## group's cells with cell i left out (its effect and its weight leave
## together), the fitted model held fixed.
.w1_loo <- function(eff, y, w, stat) {
  vapply(seq_along(eff), function(i) stat(eff[-i], y[-i], w[-i]),
         numeric(1))
}

## BCa interval (Efron 1987; Efron and Tibshirani 1993, ch. 14) from the
## estimate, the bootstrap replicates and the leave-one-out values:
## influence values U_i = (n - 1) * (mean(loo) - loo_i), acceleration
## a = sum(U^3) / (6 * sum(U^2)^1.5), bias correction z0 = qnorm(share of
## replicates below the estimate). The interval is the replicates'
## quantiles (R's default, type 7) at pnorm(z0 + (z0 + z) / (1 - a *
## (z0 + z))) for z = qnorm(0.025) and qnorm(0.975).
.w1_bca <- function(est, boot, loo, level = 0.95) {
  n  <- length(loo)
  U  <- (n - 1) * (mean(loo) - loo)
  a  <- sum(U^3) / (6 * sum(U^2)^1.5)
  z0 <- stats::qnorm(mean(boot < est))
  z  <- stats::qnorm(c((1 - level) / 2, (1 + level) / 2))
  p  <- stats::pnorm(z0 + (z0 + z) / (1 - a * (z0 + z)))
  stats::quantile(boot, p, names = FALSE)
}


## -- jackknife: overall ATT over a window ------------------------------------

test_that("W1, jackknife: overall ATT over a window, weighted estimate and Tukey SE", {
  skip_on_cran()
  fit <- .w1_fit("jackknife")
  W   <- .w1_weights(fit)
  tr  <- .w1_treated(fit)
  N   <- ncol(fit$Y.dat)
  z   <- stats::qnorm(0.975)
  ## the fit keeps the data's weights, and replicate j is the fit without
  ## unit j (its columns are the other units, in panel order)
  expect_equal(unname(fit[["W.agg", exact = TRUE]]), W)
  expect_true(all(vapply(seq_len(N), function(j) {
    identical(as.integer(fit$colnames.boot[[j]]), seq_len(N)[-j])
  }, logical(1))))
  ## over all treated cells, the hand replicates are the fit's own
  ## leave-one-out overall ATTs (a check of the construction)
  expect_equal(.w1_jack_reps(fit, W, tr, .w1_att),
               as.vector(fit$att.avg.boot), tolerance = 1e-10)
  for (win in list(c(1, 5), c(4, 8))) {
    m   <- tr & fit$T.on >= win[1] & fit$T.on <= win[2]
    est <- .w1_att(fit$eff[m], fit$Y.dat[m], W[m])
    th  <- .w1_jack_reps(fit, W, m, .w1_att)
    se  <- .w1_jack_se(th)
    ov  <- fect::estimand(fit, "att", "overall", window = win,
                          ci.method = "normal")
    info <- paste0("window = c(", win[1], ", ", win[2], ")")
    ## the window keeps only part of the treated cells, so estimand()
    ## recomputes every replicate; every replicate has cells
    expect_equal(ov$n_cells, sum(m), info = info)
    expect_lt(sum(m), sum(tr))
    expect_true(all(is.finite(th)), info = info)
    ## weighted estimate, Tukey SE of the weighted replicates, normal
    ## interval. b1dded6 gave the unweighted mean and the SE of unweighted
    ## replicates: c(1, 5) -1.147 (SE 3.690), weighted -2.850 (SE 4.356);
    ## c(4, 8) 1.772 (SE 4.586), weighted 0.020 (SE 6.616)
    expect_equal(ov$estimate, est, tolerance = 1e-10, info = info)
    expect_equal(ov$se, se, tolerance = 1e-8, info = info)
    expect_equal(c(ov$ci.lo, ov$ci.hi), est + c(-z, z) * se,
                 tolerance = 1e-8, info = info)
  }
})


## -- jackknife: APTT by event time -------------------------------------------

test_that("W1, jackknife: APTT by event time, weighted estimate and Tukey SE", {
  skip_on_cran()
  fit <- .w1_fit("jackknife")
  W   <- .w1_weights(fit)
  tr  <- .w1_treated(fit)
  ets <- sort(unique(fit$T.on[tr]))
  ## per event time: weighted mean effect over weighted mean Y(0), in the
  ## estimate and in every leave-one-out replicate
  ref <- vapply(ets, function(e) {
    m  <- tr & fit$T.on == e
    th <- .w1_jack_reps(fit, W, m, .w1_aptt)
    c(n = sum(m), finite = all(is.finite(th)),
      est = .w1_aptt(fit$eff[m], fit$Y.dat[m], W[m]), se = .w1_jack_se(th))
  }, numeric(4))
  ap <- fect::estimand(fit, "aptt", "event.time", ci.method = "normal")
  z  <- stats::qnorm(0.975)
  ## b1dded6 (unweighted) at event time 2: -0.0108 (SE 0.0534); weighted
  ## -0.0365 (SE 0.0590)
  expect_equal(ap$event.time, ets)
  expect_equal(ap$n_cells, as.integer(ref["n", ]))
  expect_true(all(ref["finite", ] == 1))
  expect_equal(ap$estimate, ref["est", ], tolerance = 1e-10)
  expect_equal(ap$se, ref["se", ], tolerance = 1e-8)
  expect_equal(ap$ci.lo, ref["est", ] - z * ref["se", ], tolerance = 1e-8)
  expect_equal(ap$ci.hi, ref["est", ] + z * ref["se", ], tolerance = 1e-8)
})


## -- jackknife: log ATT by event time ----------------------------------------

test_that("W1, jackknife: log ATT by event time, weighted estimate and Tukey SE", {
  skip_on_cran()
  fit <- .w1_fit("jackknife")
  W   <- .w1_weights(fit)
  tr  <- .w1_treated(fit)
  ets <- sort(unique(fit$T.on[tr]))
  ## per event time: weighted mean of log(Y) - log(Y(0)), in the estimate
  ## and in every leave-one-out replicate. The outcome and every imputed
  ## Y(0) are positive, so no cell is dropped (a non-positive one would
  ## make its replicate non-finite)
  ref <- vapply(ets, function(e) {
    m  <- tr & fit$T.on == e
    th <- .w1_jack_reps(fit, W, m, .w1_latt)
    c(n = sum(m), finite = all(is.finite(th)),
      est = .w1_latt(fit$eff[m], fit$Y.dat[m], W[m]), se = .w1_jack_se(th))
  }, numeric(4))
  lg <- fect::estimand(fit, "log.att", "event.time", ci.method = "normal")
  z  <- stats::qnorm(0.975)
  ## b1dded6 (unweighted) at event time 2: -0.0128 (SE 0.0546); weighted
  ## -0.0396 (SE 0.0624)
  expect_equal(lg$event.time, ets)
  expect_equal(lg$n_cells, as.integer(ref["n", ]))
  expect_true(all(ref["finite", ] == 1))
  expect_equal(lg$estimate, ref["est", ], tolerance = 1e-10)
  expect_equal(lg$se, ref["se", ], tolerance = 1e-8)
  expect_equal(lg$ci.lo, ref["est", ] - z * ref["se", ], tolerance = 1e-8)
  expect_equal(lg$ci.hi, ref["est", ] + z * ref["se", ], tolerance = 1e-8)
})


## -- bootstrap: BCa interval of the overall ATT ------------------------------

test_that("W1, bootstrap: BCa interval of the overall ATT, weighted replicates and acceleration", {
  skip_on_cran()
  fit  <- .w1_fit("bootstrap")
  W    <- .w1_weights(fit)
  tr   <- .w1_treated(fit)
  est  <- .w1_att(fit$eff[tr], fit$Y.dat[tr], W[tr])
  boot <- .w1_boot_reps(fit, W, tr, .w1_att)
  loo  <- .w1_loo(fit$eff[tr], fit$Y.dat[tr], W[tr], .w1_att)
  ci   <- .w1_bca(est, boot, loo)
  ov   <- fect::estimand(fit, "att", "overall", ci.method = "bca")
  expect_equal(unname(fit[["W.agg", exact = TRUE]]), W)
  ## the hand replicates are the fit's own replicate overall ATTs (a check
  ## of the construction)
  expect_equal(boot, as.vector(fit$att.avg.boot), tolerance = 1e-10)
  ## the BCa interval covers the estimate, so estimand() keeps it (it
  ## falls back to the normal interval otherwise)
  expect_true(ci[1] < est && est < ci[2])
  ## weighted estimate, SD of the weighted replicates as the SE, and the
  ## BCa interval with the weighted acceleration. b1dded6 (unweighted):
  ## 0.826 in [-5.295, 4.628]; weighted: -0.502 (att.avg) in
  ## [-7.449, 4.063]
  expect_equal(ov$estimate, est, tolerance = 1e-10)
  expect_equal(ov$se, stats::sd(boot), tolerance = 1e-8)
  expect_equal(c(ov$ci.lo, ov$ci.hi), ci, tolerance = 1e-8)
  ## the check sees the acceleration's weights: with equal weights in the
  ## leave-one-out values (same estimate and replicates), the interval
  ## would move, here to [-7.355, 4.215]
  loo_eq <- .w1_loo(fit$eff[tr], fit$Y.dat[tr], rep(1, sum(tr)), .w1_att)
  expect_gt(max(abs(c(ov$ci.lo, ov$ci.hi) - .w1_bca(est, boot, loo_eq))),
            1e-6)
})


## -- bootstrap: BCa intervals of APTT and log ATT by event time --------------

test_that("W1, bootstrap: BCa intervals of APTT and log ATT by event time", {
  skip_on_cran()
  fit <- .w1_fit("bootstrap")
  W   <- .w1_weights(fit)
  tr  <- .w1_treated(fit)
  ets <- 1:5    # at least six treated states each: every replicate has cells
  for (type in c("aptt", "log.att")) {
    stat <- if (type == "aptt") .w1_aptt else .w1_latt
    ref <- vapply(ets, function(e) {
      m    <- tr & fit$T.on == e
      est  <- stat(fit$eff[m], fit$Y.dat[m], W[m])
      boot <- .w1_boot_reps(fit, W, m, stat)
      loo  <- .w1_loo(fit$eff[m], fit$Y.dat[m], W[m], stat)
      ci   <- .w1_bca(est, boot, loo)
      ## the same estimate and replicates with equal weights in the
      ## leave-one-out values
      ci_eq <- .w1_bca(est, boot, .w1_loo(fit$eff[m], fit$Y.dat[m],
                                          rep(1, sum(m)), stat))
      c(finite = all(is.finite(boot)), est = est, se = stats::sd(boot),
        lo = ci[1], hi = ci[2], lo_eq = ci_eq[1], hi_eq = ci_eq[2])
    }, numeric(7))
    out <- fect::estimand(fit, type, "event.time", ci.method = "bca")
    k   <- match(ets, out$event.time)
    ## b1dded6 (unweighted), APTT at event time 2: -0.0108 in
    ## [-0.1065, 0.0747]; weighted: -0.0365 in [-0.1456, 0.0289]
    expect_true(all(ref["finite", ] == 1), info = type)
    ## each BCa interval covers its estimate, so estimand() keeps it
    expect_true(all(ref["lo", ] < ref["est", ] & ref["est", ] < ref["hi", ]),
                info = type)
    expect_equal(out$estimate[k], ref["est", ], tolerance = 1e-10,
                 info = type)
    expect_equal(out$se[k], ref["se", ], tolerance = 1e-8, info = type)
    expect_equal(out$ci.lo[k], ref["lo", ], tolerance = 1e-8, info = type)
    expect_equal(out$ci.hi[k], ref["hi", ], tolerance = 1e-8, info = type)
    ## at every event time the interval would move with equal weights in
    ## the acceleration
    shift <- pmax(abs(out$ci.lo[k] - ref["lo_eq", ]),
                  abs(out$ci.hi[k] - ref["hi_eq", ]))
    expect_true(all(shift > 1e-6), info = type)
  }
})
