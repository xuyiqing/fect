## ---------------------------------------------------------------
## fect 2.4.7: p-values with ci.method = "basic" under parametric
## inference. Parametric draws of an effect are simulated with no effect,
## so they are centered at zero; the p-value compares the estimate with
## them. The test_that block below fails on 4faaf77 (and b1dded6) and
## passes after the fix. Self-contained: helpers prefixed .pfp_.
## ---------------------------------------------------------------

.pfp_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## One of fect's datasets, loaded into a local environment.
.pfp_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "fect", envir = e)
  e[[name]]
}

## Two-sided p-value from draws simulated under no effect: twice the
## smaller share of the centered draws at or beyond the estimate.
.pfp_null_p <- function(est, draws) {
  d <- draws - mean(draws, na.rm = TRUE)
  d <- d[!is.na(d)]
  min(1, 2 * min(mean(d >= est), mean(d <= est)))
}

## Two-sided p-value from draws centered at the estimate (bootstrap draws
## of an effect; coefficient draws under both schemes): twice the smaller
## share of the draws on either side of zero.
.pfp_pct_p <- function(draws) {
  d <- draws[!is.na(draws)]
  min(1, 2 * min(mean(d >= 0), mean(d <= 0)))
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
               .pfp_null_p(pb$att.avg, c(pb$att.avg.boot)))
  expect_lt(pb$est.avg[1, "p.value"], 0.05)
  expect_equal(unname(pb$est.att[rows, "p.value"]),
               vapply(rows, function(i) {
                 .pfp_null_p(pb$att[i], pb$att.boot[i, ])
               }, numeric(1)))
  expect_equal(unname(pb$est.placebo[1, "p.value"]),
               .pfp_null_p(pb$est.placebo[1, 1], c(pb$att.placebo.boot)))
  ## the p-values go with the basic intervals: at event times 2 to 10 the
  ## effect is clear (ATT / S.E. above 3.6), and both reject no effect
  ## (4faaf77: p from 0.52 to 0.96 there)
  clear <- pb$est.att[as.character(2:10), , drop = FALSE]
  expect_true(all(clear[, "p.value"] < 0.05))
  expect_true(all(clear[, "CI.lower"] > 0))
  ## unchanged: coefficients (their draws are centered at the estimate),
  ## the normal p-values, and bootstrap fits
  expect_equal(unname(pb$est.beta[, "p.value"]),
               apply(pb$beta.boot, 1, .pfp_pct_p))
  pn <- fit("parametric", "normal")
  expect_equal(unname(pn$est.att[, "p.value"]),
               unname((1 - stats::pnorm(abs(pn$est.att[, "ATT"] /
                                             pn$est.att[, "S.E."]))) * 2))
  bb <- fit("bootstrap", "basic")
  rows_b <- which(!is.na(bb$att))
  expect_equal(unname(bb$est.att[rows_b, "p.value"]),
               apply(bb$att.boot[rows_b, , drop = FALSE], 1, .pfp_pct_p))
})
