## ---------------------------------------------------------------
## fect #169: with ci.method = "basic", a parametric fit's p-values are
## dual to the basic interval printed beside them, as a bootstrap fit's
## are since #158. Fails on fad9f72 (fect 2.4.7 on dev), passes after.
## Self-contained: helpers prefixed .pb_. parallel = FALSE throughout.
## ---------------------------------------------------------------

.pb_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## One of fect's datasets, loaded into a local environment.
.pb_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "fect", envir = e)
  e[[name]]
}

## The basic interval at level 1 - a from draws b. With shift = TRUE the
## draws are first recentered at est (their mean moved to est), which is
## what .basic_ci_shifted_one() does for a parametric fit; quantile() type
## 7, the default, on the non-missing draws. And whether an interval
## excludes 0.
.pb_basic_ci <- function(est, b, a, shift = TRUE) {
  est <- unname(est)
  if (shift) b <- b - mean(b, na.rm = TRUE) + est
  q <- stats::quantile(b, c(1 - a / 2, a / 2), na.rm = TRUE, names = FALSE)
  c(2 * est - q[1], 2 * est - q[2])
}
.pb_excludes0 <- function(ci) ci[1] > 0 || ci[2] < 0

## Number of (row, level) pairs, over 999 levels from 0.001 to 0.999, where
## "p < level" and "the basic interval at level 1 - level excludes 0"
## disagree. est: estimates; p: p-values; boots: one row of draws per
## estimate.
.pb_levels <- seq(0.001, 0.999, by = 0.001)
.pb_disagreements <- function(est, p, boots, shift = TRUE) {
  sum(vapply(seq_along(est), function(k) {
    sum(vapply(.pb_levels, function(a) {
      (p[k] < a) != .pb_excludes0(.pb_basic_ci(est[k], boots[k, ], a, shift))
    }, logical(1)))
  }, integer(1)))
}

## The fits, each made once and reused across the blocks below.
## .pb_fit_issue(): the issue's fit. .pb_fit_x3(): the same with a pure
## noise covariate X3, whose coefficient is near 0, so that the
## coefficient p-value has room to disagree with the interval (in the
## issue's fit both coefficients have p = 0 under any rule).
.pb_memo <- new.env()
.pb_fit_issue <- function() {
  if (is.null(.pb_memo$issue)) {
    sim_gsynth <- .pb_data("sim_gsynth")
    .pb_memo$issue <- .pb_quiet(fect::fect(
      Y ~ D + X1 + X2, data = sim_gsynth, index = c("id", "time"),
      method = "gsynth", force = "two-way", r = 2, CV = FALSE, se = TRUE,
      vartype = "parametric", ci.method = "basic", nboots = 200, seed = 11,
      keep.sims = TRUE, parallel = FALSE))
  }
  .pb_memo$issue
}
.pb_fit_x3 <- function() {
  if (is.null(.pb_memo$x3)) {
    d <- .pb_data("sim_gsynth")
    set.seed(3)
    d$X3 <- stats::rnorm(nrow(d))
    .pb_memo$x3 <- .pb_quiet(fect::fect(
      Y ~ D + X1 + X2 + X3, data = d, index = c("id", "time"),
      method = "gsynth", force = "two-way", r = 2, CV = FALSE, se = TRUE,
      vartype = "parametric", ci.method = "basic", nboots = 200, seed = 11,
      keep.sims = TRUE, parallel = FALSE))
  }
  .pb_memo$x3
}


test_that("the issue's fit: at event time 1 the 95% interval excludes 0 and p is below 0.05", {
  skip_on_cran()
  fit <- .pb_fit_issue()
  row <- fit$est.att["1", ]
  ## the estimate, SE and interval are as in the issue: only the p-value changes
  expect_equal(unname(row["ATT"]), 1.23619259, tolerance = 1e-6)
  expect_equal(unname(row["S.E."]), 0.68209461, tolerance = 1e-6)
  expect_equal(unname(row[c("CI.lower", "CI.upper")]), c(0.05481349, 2.70937133),
               tolerance = 1e-6)
  expect_true(.pb_excludes0(row[c("CI.lower", "CI.upper")]))
  ## the base gives 0.05 (the counting rule); the dual rule gives 0.0420676
  expect_lt(unname(row["p.value"]), 0.05)
  expect_equal(unname(row["p.value"]), 0.0420676, tolerance = 1e-6)
  k <- which(rownames(fit$est.att) == "1")
  b <- fit$att.boot[k, ]
  th <- fit$est.att[k, "ATT"]
  expect_equal(unname(row["p.value"]),
               fect:::.pvalue_basic_dual(th, b - mean(b) + th), tolerance = 1e-12)
})


test_that("effects of a parametric fit: p < alpha exactly when the basic interval at level 1 - alpha excludes 0", {
  skip_on_cran()
  fit <- .pb_fit_issue()
  est <- fit$est.att[, "ATT"]
  p <- fit$est.att[, "p.value"]
  ab <- fit$att.boot
  expect_equal(nrow(ab), length(est))
  ## the fit's intervals are the shifted basic intervals of its stored draws
  for (k in seq_along(est)) {
    expect_equal(unname(fit$est.att[k, c("CI.lower", "CI.upper")]),
                 .pb_basic_ci(est[k], ab[k, ], 0.05), tolerance = 1e-10)
  }
  ## every p-value is the dual rule on the same shifted draws
  expect_equal(unname(p), vapply(seq_along(est), function(k) {
    fect:::.pvalue_basic_dual(est[k], ab[k, ] - mean(ab[k, ]) + est[k])
  }, numeric(1)), tolerance = 1e-12)
  ## no (row, level) pair disagrees over 999 levels (the base has 65)
  expect_identical(.pb_disagreements(est, p, ab), 0L)
  ## the overall ATT follows the same rule
  avg <- fit$est.avg[1, "ATT.avg"]
  expect_equal(unname(fit$est.avg[1, c("CI.lower", "CI.upper")]),
               .pb_basic_ci(avg, fit$att.avg.boot, 0.05), tolerance = 1e-10)
  expect_identical(.pb_disagreements(avg, fit$est.avg[1, "p.value"],
                                     matrix(fit$att.avg.boot, 1)), 0L)
})


test_that("coefficients of a parametric fit: the same duality for est.beta", {
  skip_on_cran()
  for (fit in list(.pb_fit_issue(), .pb_fit_x3())) {
    est <- fit$est.beta[, "Coef"]
    p <- fit$est.beta[, "p.value"]
    bb <- fit$beta.boot
    expect_equal(nrow(bb), length(est))
    for (k in seq_along(est)) {
      expect_equal(unname(fit$est.beta[k, c("CI.lower", "CI.upper")]),
                   .pb_basic_ci(est[k], bb[k, ], 0.05), tolerance = 1e-10)
    }
    expect_equal(unname(p), vapply(seq_along(est), function(k) {
      fect:::.pvalue_basic_dual(est[k], bb[k, ] - mean(bb[k, ]) + est[k])
    }, numeric(1)), tolerance = 1e-12)
    expect_identical(.pb_disagreements(est, p, bb), 0L)
  }
  ## the noise coefficient is the one with room to disagree: its interval
  ## covers 0 and its p-value is well inside (0, 1)
  x3 <- .pb_fit_x3()$est.beta["X3", ]
  expect_false(.pb_excludes0(x3[c("CI.lower", "CI.upper")]))
  expect_gt(unname(x3["p.value"]), 0.05)
  expect_lt(unname(x3["p.value"]), 0.95)
})


test_that("a bootstrap fit keeps its p-values: the dual rule on the draws as they are", {
  skip_on_cran()
  sim_gsynth <- .pb_data("sim_gsynth")
  fit <- .pb_quiet(fect::fect(
    Y ~ D + X1 + X2, data = sim_gsynth, index = c("id", "time"),
    method = "gsynth", force = "two-way", r = 2, CV = FALSE, se = TRUE,
    vartype = "bootstrap", ci.method = "basic", nboots = 50, seed = 5,
    keep.sims = TRUE, parallel = FALSE))
  for (slot in list(list(e = fit$est.att[, "ATT"], p = fit$est.att[, "p.value"],
                         ci = fit$est.att[, c("CI.lower", "CI.upper"), drop = FALSE],
                         b = fit$att.boot),
                    list(e = fit$est.beta[, "Coef"], p = fit$est.beta[, "p.value"],
                         ci = fit$est.beta[, c("CI.lower", "CI.upper"), drop = FALSE],
                         b = fit$beta.boot))) {
    for (k in seq_along(slot$e)) {
      expect_equal(unname(slot$ci[k, ]),
                   .pb_basic_ci(slot$e[k], slot$b[k, ], 0.05, shift = FALSE),
                   tolerance = 1e-10)
    }
    expect_equal(unname(slot$p), vapply(seq_along(slot$e), function(k) {
      fect:::.pvalue_basic_dual(slot$e[k], slot$b[k, ])
    }, numeric(1)), tolerance = 1e-12)
    expect_identical(.pb_disagreements(slot$e, slot$p, slot$b, shift = FALSE), 0L)
  }
})


test_that("missing draws: the parametric rule gives NA for all-NA draws and for fewer than two finite draws", {
  skip_on_cran()
  ## the parametric branch applies the dual rule to the draws recentered at
  ## the estimate; recentering all-NA draws leaves them NA
  shifted <- function(th, b) b - mean(b, na.rm = TRUE) + th
  expect_identical(fect:::.pvalue_basic_dual(1, shifted(1, rep(NA_real_, 5))), NA_real_)
  expect_identical(fect:::.pvalue_basic_dual(1, shifted(1, c(0.2, NA, NA, NaN))), NA_real_)
  expect_identical(fect:::.pvalue_basic_dual(NA_real_, shifted(NA_real_, c(0.2, 0.4))), NA_real_)
  two <- fect:::.pvalue_basic_dual(1, shifted(1, c(0.2, NA, 0.4)))
  expect_true(is.finite(two))
  expect_gte(two, 0)
  expect_lte(two, 1)
})
