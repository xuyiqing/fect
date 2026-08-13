## ---------------------------------------------------------------
## Tests for estimand(test = "placebo" / "carryover").
##
## Issue #131 (ajunquera): pre-treatment APTT estimates for the
## alternative-estimand API. Generalized to all four shipped estimand
## types (att, aptt, log.att; att.cumu is intentionally disallowed).
## ---------------------------------------------------------------

suppressWarnings(data("simdata", package = "fect"))


## DGP helper duplicated from test-book-claims.R so this file is
## self-contained without cross-file source ordering.
.make_panel <- function(N = 40, TT = 20, T0 = 12, Ntr = 12,
                        tau = 3.0, seed = 9301, reversals = FALSE) {
    set.seed(seed)
    alpha_i <- rnorm(N, 0, 2)
    xi_t    <- rnorm(TT, 0, 1)
    D <- matrix(0L, TT, N)
    D[(T0 + 1):TT, 1:Ntr] <- 1L
    if (reversals) {
        for (i in 1:min(3, Ntr)) D[(TT - 1):TT, i] <- 0L
    }
    eps <- matrix(rnorm(N * TT, 0, 1), TT, N)
    Y <- outer(xi_t, rep(1, N)) +
         outer(rep(1, TT), alpha_i) +
         tau * D + eps
    data.frame(
        id   = rep(1:N, each = TT),
        time = rep(1:TT, N),
        Y    = as.vector(Y),
        D    = as.vector(D)
    )
}


.fit_placebo <- function(nboots = 50) {
    d <- .make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12, seed = 9301)
    suppressWarnings(suppressMessages(
        fect::fect(Y ~ D, data = d, index = c("id", "time"),
                   method = "fe", se = TRUE, nboots = nboots,
                   parallel = FALSE, keep.sims = TRUE,
                   placeboTest = TRUE, placebo.period = c(-2, 0),
                   CV = FALSE)
    ))
}

.fit_carryover <- function(nboots = 50) {
    d <- .make_panel(N = 40, TT = 20, T0 = 12, Ntr = 12,
                     seed = 9302, reversals = TRUE)
    suppressWarnings(suppressMessages(
        fect::fect(Y ~ D, data = d, index = c("id", "time"),
                   method = "fe", se = TRUE, nboots = nboots,
                   parallel = FALSE, keep.sims = TRUE,
                   carryoverTest = TRUE, carryover.period = c(1, 2),
                   CV = FALSE)
    ))
}


## -- PC.1  test argument is documented and validated ---------------

test_that("PC.1: estimand() rejects unknown test values", {

  skip_on_cran()

  fit <- .fit_placebo()
  expect_error(
    fect::estimand(fit, "att", "event.time", test = "bogus"),
    "should be one of"
  )
})


test_that("PC.2: test='placebo' errors when fit had placeboTest=FALSE", {

  skip_on_cran()

  fit_plain <- suppressWarnings(suppressMessages(
    fect::fect(Y ~ D, data = simdata, index = c("id", "time"),
               method = "fe", force = "two-way",
               se = TRUE, nboots = 50, parallel = FALSE,
               keep.sims = TRUE)
  ))
  expect_error(
    fect::estimand(fit_plain, "aptt", "event.time", test = "placebo"),
    "placeboTest = TRUE"
  )
  ## Also verify the message names the actionable refit hint.
  expect_error(
    fect::estimand(fit_plain, "aptt", "event.time", test = "placebo"),
    "Refit with"
  )
})


test_that("PC.3: test='carryover' errors when fit had carryoverTest=FALSE", {

  skip_on_cran()

  fit_plain <- suppressWarnings(suppressMessages(
    fect::fect(Y ~ D, data = simdata, index = c("id", "time"),
               method = "fe", force = "two-way",
               se = TRUE, nboots = 50, parallel = FALSE,
               keep.sims = TRUE)
  ))
  expect_error(
    fect::estimand(fit_plain, "aptt", "event.time", test = "carryover"),
    "carryoverTest = TRUE"
  )
  expect_error(
    fect::estimand(fit_plain, "aptt", "event.time", test = "carryover"),
    "Refit with"
  )
})


## -- PC.4  test='placebo' returns event-time series in placebo window

test_that("PC.4: placebo aptt returns rows in fit$placebo.period range", {

  skip_on_cran()

  fit <- .fit_placebo()
  est <- fect::estimand(fit, "aptt", "event.time", test = "placebo")

  expect_s3_class(est, "data.frame")
  expect_setequal(names(est),
                  c("event.time", "estimate", "se",
                    "ci.lo", "ci.hi", "n_cells", "vartype"))

  pp <- fit$placebo.period
  expect_true(all(est$event.time >= pp[1] & est$event.time <= pp[2]))
  expect_true(all(est$n_cells > 0L))
})


## -- PC.5  test='placebo' covers att (log.att is gated by v2.4.2 hard-error)

test_that("PC.5: placebo att returns rows in placebo window; placebo log.att hard-errors when DGP triggers cell-drop pathology", {

  skip_on_cran()

  fit <- .fit_placebo()

  est_att <- fect::estimand(fit, "att", "event.time", test = "placebo")

  pp <- fit$placebo.period
  expect_true(all(est_att$event.time >= pp[1] & est_att$event.time <= pp[2]))

  ## log.att on simdata triggers v2.4.2's point-level hard-error because
  ## simdata has Y <= 0 cells. This is the EXPECTED behavior --- the
  ## point-level check fires before any bootstrap work, with a clearer
  ## actionable message than the bootstrap-level pathology message.
  expect_error(
    fect::estimand(fit, "log.att", "event.time", test = "placebo"),
    "log\\.att requires Y > 0"
  )
})


## -- PC.6  byte-equality: placebo att rows match fit$est.att rows -----

test_that("PC.6: estimand(att, test=placebo) matches fit$est.att rows", {

  skip_on_cran()

  fit <- .fit_placebo()
  est <- fect::estimand(fit, "att", "event.time", test = "placebo")

  ## fit$est.att is per-event-time over all (treated post + masked
  ## placebo) cells. Restrict to placebo event times and compare.
  ea <- fit$est.att
  et <- as.numeric(rownames(ea))
  pp <- fit$placebo.period
  in_pp <- et >= pp[1] & et <= pp[2]

  ## Sort both by event.time for safety.
  est <- est[order(est$event.time), , drop = FALSE]
  ref <- ea[in_pp, , drop = FALSE]
  ref <- ref[order(as.numeric(rownames(ref))), , drop = FALSE]

  expect_equal(est$event.time,
               as.numeric(rownames(ref)))
  expect_equal(est$estimate,
               unname(ref[, "ATT"]),
               tolerance = 1e-10)
  expect_equal(est$n_cells,
               as.integer(unname(ref[, "count"])))
})


## -- PC.7  type='att.cumu' is rejected with non-none test ------------

test_that("PC.7: type='att.cumu' is incompatible with test != 'none'", {

  skip_on_cran()

  fit <- .fit_placebo()
  expect_error(
    fect::estimand(fit, "att.cumu", "event.time", test = "placebo"),
    "incompatible with test"
  )
})


## -- PC.8  carryover path works on a reversal panel ------------------

test_that("PC.8: carryover aptt returns rows in carryover window", {

  skip_on_cran()

  fit <- .fit_carryover()
  est <- fect::estimand(fit, "aptt", "event.time", test = "carryover")

  cp <- fit$carryover.period
  expect_true(all(est$event.time >= cp[1] & est$event.time <= cp[2]))
  expect_true(all(est$n_cells > 0L))
})


## -- PC.9  test='placebo' auto-overrides direction to 'on' -----------

test_that("PC.9: test='placebo' silently uses direction='on'", {

  skip_on_cran()

  fit <- .fit_placebo()
  ## Even when user passes direction='off', placebo path uses T.on
  ## (so we should get an event-time series, not error on missing T.off).
  est <- fect::estimand(fit, "aptt", "event.time",
                        test = "placebo", direction = "off")
  pp <- fit$placebo.period
  expect_true(all(est$event.time >= pp[1] & est$event.time <= pp[2]))
})


## -- PC.10 default vartype is sourced from fit, not from arg ---------

test_that("PC.10: vartype column reports method used at fit time", {

  skip_on_cran()

  fit <- .fit_placebo()
  est <- fect::estimand(fit, "aptt", "event.time", test = "placebo")
  expect_true(unique(est$vartype) %in%
              c("bootstrap", "jackknife", "parametric"))
})
