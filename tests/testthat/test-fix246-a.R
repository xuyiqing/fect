## ---------------------------------------------------------------
## fect 2.4.6 correctness fixes, group A (inference and bootstrap).
## Run 2026-09-24-fix246-correctness; one block per item (A1-A7).
## Every test here fails on dev @ 412d7ae and passes after the fix.
## ---------------------------------------------------------------

## A hand-built fit with 3 periods and 4 units. Units 1 and 4 are never
## treated; unit 2 is treated in period 3 and unit 3 from period 2. Its three
## bootstrap replicates have columns in resampled order (colnames.boot), so
## replicate column k is NOT unit k:
##   b = 1: units (3, 3, 1, 4)  unit 3 drawn twice, unit 2 not drawn
##   b = 2: units (2, 3, 4, 4)
##   b = 3: units (2, 2, 1, 1)  unit 2 drawn twice, unit 3 not drawn
.fix246_hand_fit <- function() {
  TT <- 3L
  N  <- 4L
  Y    <- matrix(c(10, 11, 12,  20, 21, 30,  40, 50, 52,  5, 6, 7), TT, N)
  D    <- matrix(c(0, 0, 0,  0, 0, 1,  0, 1, 1,  0, 0, 0), TT, N)
  T.on <- matrix(c(NA, NA, NA,  -1, 0, 1,  0, 1, 2,  NA, NA, NA), TT, N)
  eff  <- matrix(c(0, 0, 0,  0.1, -0.1, 8,  0.2, 9, 10,  0, 0, 0), TT, N)
  eb <- array(0, dim = c(TT, N, 3L))
  eb[, 1, 1] <- c(0.3, 7, 9)     # unit 3
  eb[, 2, 1] <- c(-0.2, 8, 11)   # unit 3 again
  eb[, 1, 2] <- c(0.1, 0, 6)     # unit 2
  eb[, 2, 2] <- c(0.4, 10, 12)   # unit 3
  eb[, 1, 3] <- c(0, 0.2, 9)     # unit 2
  eb[, 2, 3] <- c(0.1, -0.3, 5)  # unit 2 again
  structure(list(
    Y.dat = Y, D.dat = D, I.dat = matrix(1, TT, N), T.on = T.on, eff = eff,
    id = c("u1", "u2", "u3", "u4"), rawtime = 1:3, hasRevs = 0,
    vartype = "bootstrap", eff.boot = eb,
    colnames.boot = list(c(3L, 3L, 1L, 4L), c(2L, 3L, 4L, 4L),
                         c(2L, 2L, 1L, 1L))
  ), class = "fect")
}

## Staggered panel whose treated units are NOT the first columns: units 7-14
## of 20 are treated, half from period 9 and half from period 11.
.fix246_panel <- function(N = 20, TT = 16, tr = 7:14, T0 = c(9, 11),
                          seed = 246) {
  set.seed(seed)
  F <- stats::rnorm(TT)
  L <- stats::rnorm(N)
  d <- expand.grid(time = seq_len(TT), id = seq_len(N))[, c("id", "time")]
  start <- rep(Inf, N)
  start[tr] <- rep(T0, length.out = length(tr))
  d$D <- as.numeric(d$time >= start[d$id])
  d$Y <- 5 + stats::rnorm(N)[d$id] + stats::rnorm(TT)[d$time] +
    F[d$time] * L[d$id] + 2 * d$D + stats::rnorm(nrow(d), sd = 0.5)
  d
}

.fix246_fit <- function(d, nboots = 40, keep.sims = TRUE, ...) {
  suppressWarnings(suppressMessages(
    fect::fect(Y ~ D, data = d, index = c("id", "time"), method = "gsynth",
               r = 1, CV = FALSE, se = TRUE, nboots = nboots, parallel = FALSE,
               keep.sims = keep.sims, seed = 1, ...)
  ))
}

## Evaluate `expr`, collecting (and muffling) its messages and warnings.
.fix246_capture <- function(expr) {
  msgs  <- character()
  warns <- character()
  value <- withCallingHandlers(
    expr,
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, messages = msgs, warnings = warns)
}


## -- A1  post-hoc estimands read each replicate's own cells ----------

test_that("A1: estimand(att, overall) reads each replicate through colnames.boot", {
  fit <- .fix246_hand_fit()
  ## Treated cells: (t3, u2), (t2, u3), (t3, u3). Replicate means over the
  ## copies of these cells that each replicate drew:
  att_b <- c(mean(c(7, 9, 8, 11)), mean(c(6, 10, 12)), mean(c(9, 5)))
  est <- fect::estimand(fit, "att", "overall", ci.method = "normal")
  expect_equal(est$estimate, 9)
  expect_equal(est$se, stats::sd(att_b), tolerance = 1e-12)
  expect_equal(est$n_cells, 3L)
})

test_that("A1: placebo event-time, aptt and log.att use the replicate's cells", {
  fit <- .fix246_hand_fit()

  ## placebo cells: T.on in [-1, 0] -> (t1, u2) at -1; (t2, u2), (t1, u3) at 0
  fit_pl <- fit
  fit_pl$placeboTest <- TRUE
  fit_pl$placebo.period <- c(-1, 0)
  pl <- fect::estimand(fit_pl, "att", "event.time", test = "placebo",
                       ci.method = "normal")
  expect_equal(pl$event.time, c(-1, 0))
  ## event time -1: b1 has no unit 2 (missing draw)
  expect_equal(pl$se[1], stats::sd(c(0.1, mean(c(0, 0.1)))), tolerance = 1e-12)
  expect_equal(pl$se[2],
               stats::sd(c(mean(c(0.3, -0.2)), mean(c(0, 0.4)),
                           mean(c(0.2, -0.3)))),
               tolerance = 1e-12)

  ## aptt = mean(eff_b) / mean(Y - eff_b) over each replicate's cells.
  ## event time 1: (t3, u2) Y = 30 and (t2, u3) Y = 50
  aptt1 <- c(mean(c(7, 8)) / mean(c(50 - 7, 50 - 8)),
             mean(c(6, 10)) / mean(c(30 - 6, 50 - 10)),
             mean(c(9, 5)) / mean(c(30 - 9, 30 - 5)))
  ## event time 2: (t3, u3) Y = 52; replicate 3 drew no unit 3
  aptt2 <- c(mean(c(9, 11)) / mean(c(52 - 9, 52 - 11)), 12 / (52 - 12))
  ap <- fect::estimand(fit, "aptt", "event.time", ci.method = "normal")
  expect_equal(ap$event.time, c(1, 2))
  expect_equal(ap$estimate, c(8.5 / mean(c(22, 41)), 10 / 42),
               tolerance = 1e-12)
  expect_equal(ap$se, c(stats::sd(aptt1), stats::sd(aptt2)), tolerance = 1e-12)

  la1 <- c(mean(log(50) - log(c(43, 42))),
           mean(log(c(30, 50)) - log(c(24, 40))),
           mean(log(30) - log(c(21, 25))))
  la2 <- c(mean(log(52) - log(c(43, 41))), log(52) - log(40))
  la <- fect::estimand(fit, "log.att", "event.time", ci.method = "normal")
  expect_equal(la$se, c(stats::sd(la1), stats::sd(la2)), tolerance = 1e-12)
})

test_that("A1: log.att's unstable-cell share counts replicates that contain the cell", {
  fit <- .fix246_hand_fit()
  ## Replicate 2's draw for (t3, u3) makes Y0_b = 52 - 60 < 0. Only
  ## replicates 1 and 2 contain that cell, so the share is 1 of 2.
  fit$eff.boot[3, 2, 2] <- 60
  expect_error(
    fect::estimand(fit, "log.att", "event.time", ci.method = "normal"),
    "1 of 2 bootstrap replicates"
  )
})

test_that("A1: imputed_outcomes(replicates = TRUE) returns each replicate's cells", {
  fit <- .fix246_hand_fit()
  po_b <- fect::imputed_outcomes(fit, replicates = TRUE)
  expect_equal(po_b$replicate, c(1, 1, 1, 1, 2, 2, 2, 3, 3))
  expect_equal(po_b$id, c("u3", "u3", "u3", "u3", "u2", "u3", "u3", "u2", "u2"))
  expect_equal(po_b$time, c(2, 3, 2, 3, 3, 2, 3, 3, 3))
  expect_equal(po_b$eff, c(7, 9, 8, 11, 6, 10, 12, 9, 5))
  expect_equal(po_b$Y0_hat, po_b$Y_obs - po_b$eff)
  expect_equal(unname(c(tapply(po_b$eff, po_b$replicate, mean))),
               c(mean(c(7, 9, 8, 11)), mean(c(6, 10, 12)), mean(c(9, 5))))
})

test_that("A1: a cells formula sees the variables of the function that built it", {
  fit <- .fix246_hand_fit()
  pick <- function(fit, k_fix246) {
    fect::imputed_outcomes(fit, cells = ~ event.time <= k_fix246)
  }
  expect_equal(nrow(pick(fit, 1)), 2L)
  overall <- function(fit, k_fix246) {
    fect::estimand(fit, "att", "overall", ci.method = "normal",
                   cells = ~ event.time <= k_fix246)
  }
  expect_equal(overall(fit, 1)$estimate, mean(c(8, 9)))
})

test_that("A1: estimand SEs equal fect's own SEs when treated units are not first", {
  skip_on_cran()
  d <- .fix246_panel()
  for (vt in c("bootstrap", "parametric")) {
    fit <- .fix246_fit(d, vartype = vt, placeboTest = TRUE,
                       placebo.period = c(-2, 0))
    est <- fect::estimand(fit, "att", "overall")
    expect_equal(est$se, unname(fit$est.avg[1, "S.E."]), tolerance = 1e-10,
                 info = vt)
    pl <- fect::estimand(fit, "att", "event.time", test = "placebo")
    expect_equal(pl$se,
                 unname(fit$est.att[as.character(pl$event.time), "S.E."]),
                 tolerance = 1e-10, info = vt)
    ## each replicate's rows average to that replicate's att.avg
    po_b <- fect::imputed_outcomes(fit, replicates = TRUE)
    expect_equal(unname(c(tapply(po_b$eff, po_b$replicate, mean))),
                 unname(c(fit$att.avg.boot)), tolerance = 1e-10, info = vt)
  }
})


## -- A2  effect() / att.cumu() use the fit's stored vartype --------------

test_that("A2: effect() reads the stored vartype when vartype was passed as a variable", {
  skip_on_cran()
  d <- .fix246_panel()
  vt <- "parametric"
  fit <- suppressWarnings(suppressMessages(
    fect::fect(Y ~ D, data = d, index = c("id", "time"), method = "gsynth",
               r = 1, CV = FALSE, se = TRUE, vartype = vt, nboots = 40,
               parallel = FALSE, keep.sims = TRUE, seed = 1)
  ))
  expect_identical(fit$vartype, "parametric")
  expect_true(is.name(fit$call$vartype))
  M <- suppressMessages(fect::effect(fit, cumu = TRUE, plot = FALSE))$effect.est.att
  ## parametric draws are centred at 0; the CI must be centred on the estimate
  expect_equal(unname((M[, "CI.lower"] + M[, "CI.upper"]) / 2),
               unname(M[, "ATT"]), tolerance = 1e-8)
  expect_true(all(M[, "CI.lower"] < M[, "ATT"] & M[, "ATT"] < M[, "CI.upper"]))
  ## the effect is 2 per period, so the cumulative effect is clearly nonzero
  expect_true(all(M[, "p.value"] < 0.01))
  est <- fect::estimand(fit, "att.cumu", "event.time", ci.method = "percentile")
  expect_equal(est$ci.lo, unname(M[, "CI.lower"]))
})

test_that("A2: att.cumu() gives parametric fits a CI around the cumulative effect", {
  skip_on_cran()
  d <- .fix246_panel()
  fit <- .fix246_fit(d, vartype = "parametric")
  cm <- suppressMessages(fect::att.cumu(fit, period = c(1, 5), plot = FALSE))
  last <- cm[nrow(cm), ]
  z <- stats::qnorm(0.975)
  expect_equal(unname(last["CI.lower"]), unname(last["catt"] - z * last["S.E."]))
  expect_equal(unname(last["CI.upper"]), unname(last["catt"] + z * last["S.E."]))
  expect_equal(unname(last["p.value"]),
               unname(2 * stats::pnorm(-abs(last["catt"] / last["S.E."]))))
  ov <- fect::estimand(fit, "att.cumu", "overall", window = c(1, 5))
  expect_equal(ov$ci.lo, unname(last["CI.lower"]))
  expect_true(ov$ci.lo < ov$estimate && ov$estimate < ov$ci.hi)
})

test_that("A2: effect() works for a single unit and for fits with one treated unit", {
  skip_on_cran()
  d <- .fix246_panel()
  fit <- .fix246_fit(d, vartype = "bootstrap")
  one <- suppressMessages(fect::effect(fit, cumu = TRUE, id = 8, plot = FALSE))
  expect_true(is.matrix(one$effect.est.att))
  expect_true(all(is.finite(one$effect.est.att[, "ATT"])))

  d1 <- .fix246_panel(N = 12, tr = 5, T0 = 9)
  fit1 <- .fix246_fit(d1, vartype = "parametric")
  eff1 <- suppressMessages(fect::effect(fit1, cumu = TRUE, plot = FALSE))
  expect_true(all(is.finite(eff1$effect.est.att[, "S.E."])))
})


## -- A3  normalize = TRUE and the parametric bootstrap --------------------

test_that("A3: parametric SEs are the same with and without normalize = TRUE", {
  skip_on_cran()
  d <- .fix246_panel()
  f0 <- .fix246_fit(d, vartype = "parametric", normalize = FALSE)
  f1 <- .fix246_fit(d, vartype = "parametric", normalize = TRUE)
  ## before 2.4.6 the normalized SEs were multiplied by sd(Y)
  expect_gt(stats::sd(d$Y), 1.5)
  expect_equal(f1$est.avg[1, "S.E."], f0$est.avg[1, "S.E."], tolerance = 1e-4)
  expect_equal(f1$est.att[, "S.E."], f0$est.att[, "S.E."], tolerance = 1e-4)
})
