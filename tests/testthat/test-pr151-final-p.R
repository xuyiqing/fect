## ---------------------------------------------------------------
## fect 2.4.7: p-values with ci.method = "basic" under parametric
## inference. Parametric draws of an effect are simulated with no effect,
## so they are centered at zero; fect recenters them at the estimate for
## the basic interval, and the p-value is the one dual to that interval
## (#169). The test_that block below fails on 4faaf77 (and b1dded6) and
## passes after the fix. Self-contained: helpers prefixed .pfp_.
## ---------------------------------------------------------------

.pfp_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## One of fect's datasets, loaded into a local environment.
.pfp_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "fect", envir = e)
  e[[name]]
}

## Parametric draws of an effect are simulated with no effect, so they are
## centered at zero. fect recenters them at the estimate (the draws minus
## their mean, plus the estimate) and reflects them for the basic
## interval; the p-value is the one dual to that interval, on the same
## recentered draws (#169; before, twice the smaller share of the centered
## draws at or beyond the estimate, which agrees with the interval only up
## to 2 / nboots).
.pfp_shifted <- function(est, draws) draws - mean(draws, na.rm = TRUE) + est
.pfp_basic_p <- function(est, draws) {
  fect:::.pvalue_basic_dual(est, .pfp_shifted(est, draws))
}
## The basic interval at level 1 - a from the recentered draws.
.pfp_shifted_ci <- function(est, draws, a) {
  q <- stats::quantile(.pfp_shifted(est, draws), c(1 - a / 2, a / 2),
                       na.rm = TRUE, names = FALSE)
  c(2 * est - q[1], 2 * est - q[2])
}


## -- P1  parametric draws and ci.method = "basic" ---------------------------

test_that("P1: parametric p-values with ci.method = 'basic' compare the estimate with the no-effect draws", {
  skip_on_cran()
  sg <- .pfp_data("simgsynth")
  fit <- function(vt, ci) {
    .pfp_quiet(fect::fect(Y ~ D + X1 + X2, data = sg, index = c("id", "time"),
                          method = "gsynth", force = "two-way", r = 2,
                          CV = FALSE, se = TRUE, vartype = vt, nboots = 50,
                          seed = 1, ci.method = ci, keep.sims = TRUE,
                          parallel = FALSE, placeboTest = TRUE,
                          placebo.period = c(-2, 0)))
  }
  pb <- fit("parametric", "basic")
  rows <- which(!is.na(pb$att))
  ## 4faaf77 compared the zero-centered draws with zero: about 1 whatever
  ## the estimate (here est.avg 0.96 for an ATT of 5.5 with S.E. 0.27)
  expect_equal(unname(pb$est.avg[1, "p.value"]),
               .pfp_basic_p(pb$att.avg, c(pb$att.avg.boot)))
  expect_lt(pb$est.avg[1, "p.value"], 0.05)
  expect_equal(unname(pb$est.att[rows, "p.value"]),
               vapply(rows, function(i) {
                 .pfp_basic_p(pb$att[i], pb$att.boot[i, ])
               }, numeric(1)))
  expect_equal(unname(pb$est.placebo[1, "p.value"]),
               .pfp_basic_p(pb$est.placebo[1, 1], c(pb$att.placebo.boot)))
  ## each is below a exactly when the basic interval at level 1 - a, from
  ## the same recentered draws, excludes zero (#169; the counting rule
  ## disagreed at some levels)
  for (i in rows) {
    for (a in c(0.01, 0.05, 0.1, 0.5)) {
      ci <- .pfp_shifted_ci(pb$att[i], pb$att.boot[i, ], a)
      expect_identical(unname(pb$est.att[i, "p.value"] < a),
                       ci[1] > 0 || ci[2] < 0)
    }
  }
  ## the p-values go with the basic intervals: at event times 2 to 10 the
  ## effect is clear (ATT / S.E. above 3.6), and both reject no effect
  ## (4faaf77: p from 0.52 to 0.96 there)
  clear <- pb$est.att[as.character(2:10), , drop = FALSE]
  expect_true(all(clear[, "p.value"] < 0.05))
  expect_true(all(clear[, "CI.lower"] > 0))
  ## coefficients: the same rule on their draws recentered at the estimate
  ## (#169; they had kept the percentile p-value); the normal p-values
  ## are unchanged
  expect_equal(unname(pb$est.beta[, "p.value"]),
               vapply(seq_len(nrow(pb$beta.boot)), function(k) {
                 .pfp_basic_p(pb$est.beta[k, "Coef"], pb$beta.boot[k, ])
               }, numeric(1)))
  pn <- fit("parametric", "normal")
  expect_equal(unname(pn$est.att[, "p.value"]),
               unname((1 - stats::pnorm(abs(pn$est.att[, "ATT"] /
                                             pn$est.att[, "S.E."]))) * 2))
  ## bootstrap fits (2.4.7, #158): the p-value goes with the basic interval,
  ## below alpha exactly when the (1 - alpha) basic interval excludes zero
  bb <- fit("bootstrap", "basic")
  rows_b <- which(!is.na(bb$att))
  for (i in rows_b) {
    for (a in c(0.01, 0.05, 0.1, 0.5)) {
      q <- stats::quantile(bb$att.boot[i, ], c(1 - a / 2, a / 2),
                           na.rm = TRUE, names = FALSE)
      expect_identical(unname(bb$est.att[i, "p.value"] < a),
                       2 * bb$att[i] - q[1] > 0 || 2 * bb$att[i] - q[2] < 0)
    }
  }
})
