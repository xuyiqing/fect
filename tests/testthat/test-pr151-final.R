## ---------------------------------------------------------------
## fect 2.4.7: final fixes on PR #151 (b1dded6), group X: jackknife and
## parametric uncertainty of the cumulative ATT, aggregation weights in
## estimand() and effect(). Every test_that block below fails on b1dded6
## and passes after the fix. Self-contained: helpers prefixed .pf_.
## ---------------------------------------------------------------

.pf_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## One of fect's datasets, loaded into a local environment.
.pf_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "fect", envir = e)
  e[[name]]
}

## turnout, gsynth without factors: 47 states, staggered adoption, so the
## number of treated states changes over event time.
.pf_turnout <- function(vartype, ...) {
  turnout <- .pf_data("turnout")
  .pf_quiet(fect::fect(
    turnout ~ policy_edr + policy_mail_in + policy_motor, data = turnout,
    index = c("abb", "year"), method = "gsynth", r = 0, CV = FALSE,
    force = "two-way", min.T0 = 5, se = TRUE, vartype = vartype,
    keep.sims = TRUE, parallel = FALSE, ...))
}


## -- S1  parametric intervals of the cumulative ATT ------------------------

test_that("S1: effect() gives parametric fits normal intervals, as att.cumu() does", {
  skip_on_cran()
  fit <- .pf_turnout("parametric", nboots = 50, seed = 2139)
  acu <- .pf_quiet(fect::att.cumu(fit, period = c(1, 10)))
  eff <- .pf_quiet(fect::effect(fit, period = c(1, 10)))$effect.est.att
  z <- stats::qnorm(0.975)
  ## b1dded6: a t critical value with nboots - 1 degrees of freedom
  ## (k = 10: [-91.03, 134.03]; att.cumu() [-88.25, 131.25])
  expect_equal(unname(eff[, "CI.lower"]), unname(eff[, "ATT"] - z * eff[, "S.E."]),
               tolerance = 1e-10)
  expect_equal(unname(eff[, "CI.upper"]), unname(eff[, "ATT"] + z * eff[, "S.E."]),
               tolerance = 1e-10)
  expect_equal(unname(eff[, "p.value"]),
               unname(2 * stats::pnorm(-abs(eff[, "ATT"] / eff[, "S.E."]))),
               tolerance = 1e-10)
  expect_equal(unname(eff[, c("S.E.", "CI.lower", "CI.upper", "p.value")]),
               unname(acu[, c("S.E.", "CI.lower", "CI.upper", "p.value")]),
               tolerance = 1e-8)
  et <- .pf_quiet(fect::estimand(fit, "att.cumu", "event.time"))
  expect_equal(et$ci.lo[match(1:10, et$event.time)], unname(eff[, "CI.lower"]),
               tolerance = 1e-10)
})
