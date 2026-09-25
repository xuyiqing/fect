## ---------------------------------------------------------------
## fix246 section B: CV, weights in CV, implied weights, plots, guards,
## never-treated se = TRUE, cumulative ATT (items B1-B13).
##
## Regression tests for the 2.4.6 correctness run. Each block fails on
## fect dev @ 412d7ae and passes after the fix. Fits are small and serial
## (parallel = FALSE); helpers use pkg::fn for non-base functions.
## ---------------------------------------------------------------

.fb_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## Run `expr`, returning its value (or the error) and every message text.
.fb_msgs <- function(expr) {
  m <- character(0)
  res <- tryCatch(
    withCallingHandlers(
      suppressWarnings(expr),
      message = function(cnd) {
        m <<- c(m, conditionMessage(cnd))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) e
  )
  list(res = res, msgs = m)
}

.fb_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "fect", envir = e)
  e[[name]]
}


## -- B1  criterion = "pc" selects the r with the lowest PC -----------------

test_that("B1: criterion = 'pc' picks argmin PC for gsynth and cfe + never-treated", {
  skip_on_cran()
  turnout <- .fb_data("turnout")
  for (m in c("gsynth", "cfe")) {
    args <- list(turnout ~ policy_edr + policy_mail_in + policy_motor,
                 data = turnout, index = c("abb", "year"), method = m,
                 force = "two-way", CV = TRUE, r = c(0, 5), criterion = "pc",
                 se = FALSE, parallel = FALSE, seed = 1)
    if (m == "cfe") args$time.component.from <- "nevertreated"
    fit <- .fb_quiet(do.call(fect::fect, args))
    cv <- fit$CV.out
    ## PC is lowest at r = 4 on turnout; the MSPE rules pick r = 2
    expect_equal(unname(cv[which.min(cv[, "PC"]), "r"]), 4, info = m)
    expect_equal(as.numeric(fit$r.cv), 4, info = m)
  }
})

test_that("B1: fect's IFE CV table under 'pc' has the PC columns, correctly labelled", {
  skip_on_cran()
  simdata <- .fb_data("simdata")
  fit <- .fb_quiet(fect::fect(Y ~ D + X1 + X2, data = simdata,
                              index = c("id", "time"), method = "ife",
                              CV = TRUE, r = c(0, 3), criterion = "pc",
                              se = FALSE, parallel = FALSE, seed = 1))
  cv <- fit$CV.out.ife
  expect_identical(colnames(cv), c("r", "sigma2", "IC", "PC", "MSPTATT", "MSE"))
  ## the two last columns hold the values computed for them (finite), and
  ## no MSPE-family column carries MSPTATT / MSE values
  expect_true(all(is.finite(cv[, "MSPTATT"])))
  expect_true(all(is.finite(cv[, "MSE"])))
  expect_true(all(cv[, "MSE"] > cv[, "MSPTATT"]))
  ## fect's own IFE path already selects argmin PC
  expect_equal(as.numeric(fit$r.cv), unname(cv[which.min(cv[, "PC"]), "r"]))
})


## -- B2  the CV range for r is capped by the number of never-treated units --

test_that("B2: CV with few never-treated units searches at most Nco - 1 factors", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  d4 <- simgsynth[simgsynth$id %in% c(101:105, 106:109), ]  # 4 never-treated
  for (m in c("gsynth", "ife", "cfe")) {
    args <- list(Y ~ D + X1 + X2, data = d4, index = c("id", "time"),
                 method = m, force = "two-way", CV = TRUE, r = c(0, 5),
                 se = FALSE, parallel = FALSE, seed = 1)
    if (m != "gsynth") args$time.component.from <- "nevertreated"
    out <- .fb_msgs(do.call(fect::fect, args))
    ## 412d7ae: "Mat::head_cols(): size out of bounds"
    expect_false(inherits(out$res, "error"), info = m)
    expect_true(any(grepl(
      "With 4 never-treated units at most 3 factor(s) can be estimated",
      out$msgs, fixed = TRUE)), info = m)
    if (!inherits(out$res, "error")) {
      expect_equal(unname(out$res$CV.out[, "r"]), 0:3, info = m)
      expect_true(out$res$r.cv <= 3, info = m)
    }
  }
})

test_that("B2: a feasible CV range is not changed by the cap", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  out <- .fb_msgs(fect::fect(Y ~ D + X1 + X2, data = simgsynth,
                             index = c("id", "time"), method = "gsynth",
                             force = "two-way", CV = TRUE, r = c(0, 5),
                             se = FALSE, parallel = FALSE, seed = 1))
  expect_false(any(grepl("factor(s) can be estimated", out$msgs, fixed = TRUE)))
  expect_equal(unname(out$res$CV.out[, "r"]), 0:5)
})

test_that("B2: panel_factor() extracts at most min(T, N) factors (C++ guard)", {
  set.seed(1)
  E <- matrix(stats::rnorm(20), 5, 4)
  ## 412d7ae: "Mat::head_cols(): size out of bounds"
  pf <- fect:::panel_factor(E, 6L)
  expect_equal(dim(pf$factor), c(5L, 4L))
  expect_equal(dim(pf$lambda), c(4L, 4L))
  expect_equal(dim(pf$FE), c(5L, 4L))
  ## r within range: unchanged meaning
  pf2 <- fect:::panel_factor(E, 2L)
  expect_equal(dim(pf2$factor), c(5L, 2L))
  pf0 <- fect:::panel_factor(E, 0L)
  expect_equal(dim(pf0$factor), c(5L, 0L))
  expect_equal(pf0$FE, matrix(0, 5, 4))
})


## -- B3  the "CV skipped" message states the real cause ----------------------

test_that("B3: a single candidate r = 0 is reported as such, not as too few records", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  for (m in c("gsynth", "cfe")) {
    args <- list(Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
                 method = m, force = "two-way", CV = TRUE, r = 0,
                 se = FALSE, parallel = FALSE, seed = 1)
    if (m == "cfe") args$time.component.from <- "nevertreated"
    out <- .fb_msgs(do.call(fect::fect, args))
    expect_false(inherits(out$res, "error"), info = m)
    expect_true(any(grepl("Only one candidate number of factors (r = 0) was given",
                          out$msgs, fixed = TRUE)), info = m)
    expect_false(any(grepl("records of treated units are too few", out$msgs,
                           fixed = TRUE)), info = m)
    ## which r is searched does not change
    expect_equal(as.numeric(out$res$r.cv), 0, info = m)
  }
})

test_that("B3: too few pre-treatment records keep their message; one never-treated unit gets its own", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  ## two pre-treatment periods, two-way FE: no factor can be estimated
  dT <- simgsynth[simgsynth$time >= 19, ]
  out <- .fb_msgs(fect::fect(Y ~ D + X1 + X2, data = dT, index = c("id", "time"),
                             method = "gsynth", force = "two-way", CV = TRUE,
                             r = c(0, 5), min.T0 = 2, se = FALSE,
                             parallel = FALSE, seed = 1))
  expect_true(any(grepl("records of treated units are too few", out$msgs,
                        fixed = TRUE)))
  expect_false(any(grepl("Only one candidate", out$msgs, fixed = TRUE)))
  ## one never-treated unit: no factor can be estimated (412d7ae: "Not a matrix.")
  d1 <- simgsynth[simgsynth$id %in% c(101:105, 106), ]
  out1 <- .fb_msgs(fect::fect(Y ~ D, data = d1, index = c("id", "time"),
                              method = "gsynth", force = "two-way", CV = TRUE,
                              r = c(0, 5), se = FALSE, parallel = FALSE,
                              seed = 1))
  expect_false(inherits(out1$res, "error"))
  expect_true(any(grepl("with 1 never-treated unit no factor can be estimated",
                        out1$msgs, fixed = TRUE)))
  expect_false(any(grepl("records of treated units are too few", out1$msgs,
                         fixed = TRUE)))
})


## -- B4  W.agg never enters the fit or its CV scoring ------------------------

## simgsynth with unit weights and 60 control cells dropped at random, so a
## fit weight would change the control-panel fit.
.fb_unbalanced_weighted <- function() {
  d0 <- .fb_data("simgsynth")
  set.seed(7)
  ids <- unique(d0$id)
  w <- stats::setNames(stats::runif(length(ids), 0.2, 5), ids)
  d0$w <- w[as.character(d0$id)]
  ctrl_rows <- which(!(d0$id %in% 101:105))
  set.seed(3)
  d0[-sample(ctrl_rows, 60), ]
}

test_that("B4: W.agg alone leaves the model fit and its CV unweighted (gsynth, cfe)", {
  skip_on_cran()
  d <- .fb_unbalanced_weighted()
  for (m in c("gsynth", "cfe")) {
    fit <- function(...) .fb_quiet(fect::fect(
      Y ~ D + X1 + X2, data = d, index = c("id", "time"), method = m,
      force = "two-way", time.component.from = "nevertreated", se = FALSE,
      parallel = FALSE, seed = 1, ...))
    u1 <- fit(r = c(2, 2), CV = TRUE)
    a1 <- fit(r = c(2, 2), CV = TRUE, W.agg = "w")
    e1 <- fit(r = c(2, 2), CV = TRUE, W.est = "w")
    ## 412d7ae: the W.agg fit equals the W.est fit (beta diff 0.015)
    expect_equal(a1$beta, u1$beta, tolerance = 1e-10, info = m)
    expect_equal(a1$eff, u1$eff, tolerance = 1e-10, info = m)
    expect_equal(a1$CV.out, u1$CV.out, tolerance = 1e-10, info = m)
    ## a fit weight still acts on this panel
    expect_gt(max(abs(e1$beta - u1$beta)), 1e-3)
    ## 412d7ae, cfe: "subscript out of bounds" even with CV = FALSE
    a0 <- fit(r = 2, CV = FALSE, W.agg = "w")
    u0 <- fit(r = 2, CV = FALSE)
    expect_equal(a0$eff, u0$eff, tolerance = 1e-10, info = m)
  }
})
