## ---------------------------------------------------------------
## fect 2.4.6: fixes ported onto PR #151 (5afa708) from the home-mac run
## 2026-09-24-fix246-correctness: B10-cfe, B11, B12, B13 and M1. Every
## test_that block below fails on 5afa708 and passes after the port.
## Self-contained: the helpers are defined here (prefix .fix246p_).
## ---------------------------------------------------------------

.fix246p_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## One of fect's datasets, loaded into a local environment.
.fix246p_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "fect", envir = e)
  e[[name]]
}

## simgsynth: 50 units x 30 periods; units 101-105 treated from period 21,
## 45 never-treated units. Y ~ D + X1 + X2, two-way FE, no CV, seed 1.
.fix246p_nt_fit <- function(...) {
  d <- .fix246p_data("simgsynth")
  .fix246p_quiet(fect::fect(Y ~ D + X1 + X2, data = d,
                            index = c("id", "time"), force = "two-way",
                            CV = FALSE, parallel = FALSE, seed = 1,
                            nboots = 10, ...))
}


## -- B10-cfe  cfe + never-treated: replicates, loo SEs, dloo ----------------

test_that("B10-cfe: bootstrap and jackknife replicates of cfe + never-treated use the never-treated model", {
  skip_on_cran()
  c0 <- .fix246p_nt_fit(method = "cfe", time.component.from = "nevertreated",
                        r = 2, se = FALSE)
  for (vt in c("bootstrap", "jackknife")) {
    cc <- .fix246p_nt_fit(method = "cfe", time.component.from = "nevertreated",
                          r = 2, se = TRUE, vartype = vt)
    g <- .fix246p_nt_fit(method = "gsynth", r = 2, se = TRUE, vartype = vt)
    ## the point estimate is the se = FALSE fit's, as before
    expect_identical(cc$att.avg, c0$att.avg, info = vt)
    ## cfe without extra fixed effects is gsynth's model; the two solvers
    ## agree to about 5e-5. On 5afa708 the replicates were not-yet-treated
    ## cfe fits (S.E. off by up to 0.19).
    expect_lt(max(abs(cc$est.att[, "S.E."] - g$est.att[, "S.E."]),
                  na.rm = TRUE), 1e-3)
    expect_lt(max(abs(cc$att.avg.boot - g$att.avg.boot)), 1e-3)
  }
})

test_that("B10-cfe: leave-one-period-out SEs of cfe + never-treated use the never-treated model", {
  skip_on_cran()
  lo <- function(...) .fix246p_nt_fit(r = 2, se = TRUE, loo = TRUE, ...)
  g  <- lo(method = "gsynth")
  a  <- lo(method = "ife", time.component.from = "nevertreated")
  cc <- lo(method = "cfe", time.component.from = "nevertreated")
  ## ife + never-treated is relabelled gsynth (PR #151, B8e): identical
  expect_identical(a$pre.est.att, g$pre.est.att)
  ## ATT column: loo refits keep time.component.from (B8c, on 5afa708);
  ## S.E. column: the refits' replicates (0.16 off on 5afa708)
  cols <- c("ATT", "S.E.")
  expect_lt(max(abs(cc$pre.est.att[, cols] - g$pre.est.att[, cols]),
                na.rm = TRUE), 1e-3)
})

test_that("B10-cfe: dloo with never-treated fixed effects stops before the bootstrap", {
  skip_on_cran()
  ## test-dloo.R's panel: 8 periods; cohorts adopt at 4 and 6; 25 never treated
  dat <- withr::with_seed(42, {
    rows <- list()
    id <- 0
    for (g in c(4, 6, Inf)) for (u in seq_len(25)) {
      id <- id + 1
      y <- stats::rnorm(1) + 0.1 * (1:8) + stats::rnorm(8, 0, 0.4)
      d <- as.integer(!is.infinite(g) & (1:8) >= g)
      rows[[length(rows) + 1]] <- data.frame(id = id, time = 1:8,
                                             Y = y + 2 * d, D = d)
    }
    do.call(rbind, rows)
  })
  ## 5afa708: ran the bootstrap, then stopped with a message about `method`
  expect_error(
    .fix246p_quiet(fect::fect(Y ~ D, data = dat, index = c("id", "time"),
                              method = "fe", force = "two-way", dloo = TRUE,
                              se = TRUE, vartype = "bootstrap", nboots = 10,
                              time.component.from = "nevertreated",
                              parallel = FALSE, seed = 1)),
    "time.component.from")
})


## -- B11  cumulative ATT = running sum of the per-period ATTs ----------------

test_that("B11: att.cumu(), effect() and estimand('att.cumu') return the running sum", {
  skip_on_cran()
  turnout <- .fix246p_data("turnout")
  fit <- .fix246p_quiet(fect::fect(
    turnout ~ policy_edr + policy_mail_in + policy_motor, data = turnout,
    index = c("abb", "year"), method = "gsynth", r = 0, CV = FALSE,
    force = "two-way", min.T0 = 5, se = TRUE, vartype = "bootstrap",
    nboots = 50, seed = 2139, keep.sims = TRUE, parallel = FALSE))
  k   <- 1:10
  att <- fit$att[match(k, fit$time)]
  n   <- fit$count[match(k, fit$time)]
  running <- cumsum(att)
  ## treated counts vary over event time (9 8 6 6 6 3 3 3 3 3), so the two
  ## definitions differ: 21.501 (running sum) vs 8.259 (5afa708) at k = 10
  eff <- .fix246p_quiet(fect::effect(fit, period = c(1, 10)))
  acu <- .fix246p_quiet(fect::att.cumu(fit, period = c(1, 10)))
  est <- .fix246p_quiet(fect::estimand(fit, "att.cumu", "event.time"))
  ov5 <- .fix246p_quiet(fect::estimand(fit, "att.cumu", "overall",
                                       window = c(1, 5)))
  expect_equal(unname(eff$effect.est.avg), running, tolerance = 1e-10)
  expect_equal(unname(acu[, 3]), running, tolerance = 1e-10)
  expect_equal(est$estimate[1:10], running, tolerance = 1e-10)
  expect_equal(ov5$estimate, running[5], tolerance = 1e-10)
  ## weighted = TRUE keeps the count-weighted number
  acw <- .fix246p_quiet(fect::att.cumu(fit, period = c(1, 10),
                                       weighted = TRUE))
  expect_equal(unname(acw[, 3]), k * cumsum(att * n) / cumsum(n),
               tolerance = 1e-10)
  ## att.cumu() and effect() use the same replicate running sums
  ## (5afa708: 37.85 vs 27.54 at k = 10)
  expect_equal(unname(acu[2:10, "S.E."]),
               unname(eff$effect.est.att[2:10, "S.E."]), tolerance = 1e-10)
  expect_error(.fix246p_quiet(fect::att.cumu(fit, period = c(1, 10),
                                             weighted = NA)),
               "must be TRUE or FALSE")
})

test_that("B11: equal counts per event time give the old numbers; fits without SEs work", {
  skip_on_cran()
  simgsynth <- .fix246p_data("simgsynth")
  fit <- .fix246p_quiet(fect::fect(Y ~ D + X1 + X2, data = simgsynth,
                                   index = c("id", "time"), method = "gsynth",
                                   force = "two-way", r = 2, CV = FALSE,
                                   se = TRUE, vartype = "bootstrap",
                                   nboots = 20, seed = 3, keep.sims = TRUE,
                                   parallel = FALSE))
  k   <- 1:10
  att <- fit$att[match(k, fit$time)]
  n   <- fit$count[match(k, fit$time)]
  expect_true(all(n == 5))
  old <- k * cumsum(att * n) / cumsum(n)
  eff <- .fix246p_quiet(fect::effect(fit, period = c(1, 10)))
  acu <- .fix246p_quiet(fect::att.cumu(fit, period = c(1, 10)))
  expect_equal(unname(eff$effect.est.avg), old, tolerance = 1e-10)
  expect_equal(unname(acu[, 3]), old, tolerance = 1e-10)
  ## a fit without SEs (5afa708: "non-numeric argument to mathematical
  ## function")
  f0 <- .fix246p_quiet(fect::fect(Y ~ D + X1 + X2, data = simgsynth,
                                  index = c("id", "time"), method = "gsynth",
                                  force = "two-way", r = 2, CV = FALSE,
                                  se = FALSE, parallel = FALSE))
  run0 <- cumsum(f0$att[match(k, f0$time)])
  a0 <- .fix246p_quiet(fect::att.cumu(f0, period = c(1, 10)))
  a1 <- .fix246p_quiet(fect::att.cumu(f0, period = c(1, 10), weighted = TRUE))
  expect_equal(unname(a0[, 3]), run0, tolerance = 1e-10)
  expect_equal(unname(a1[, 3]), unname(a0[, 3]), tolerance = 1e-10)
  ov <- .fix246p_quiet(fect::estimand(f0, "att.cumu", "overall",
                                      window = c(1, 5)))
  expect_equal(ov$estimate, run0[5], tolerance = 1e-10)
  expect_true(is.na(ov$se))
})
