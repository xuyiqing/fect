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


## -- B5  wgt.implied: Xu (2017)'s weights, Nco x Ntr in both modes ----------

test_that("B5: wgt.implied solves t(Lco) W = t(Ltr) and rebuilds the counterfactual", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  for (m in c("gsynth", "cfe")) {
    args <- list(Y ~ D, data = simgsynth, index = c("id", "time"), method = m,
                 force = "none", r = 2, CV = FALSE, se = FALSE,
                 parallel = FALSE)
    if (m == "cfe") args$time.component.from <- "nevertreated"
    fit <- .fb_quiet(do.call(fect::fect, args))
    W   <- fit$wgt.implied
    tr  <- which(colSums(fit$D.dat) > 0)
    co  <- which(colSums(fit$D.dat) == 0)
    expect_equal(dim(W), c(length(co), length(tr)), info = m)
    expect_identical(rownames(W), as.character(fit$id[co]), info = m)
    expect_identical(colnames(W), as.character(fit$id[tr]), info = m)
    ## the defining equation (412d7ae misses by 40)
    expect_equal(unname(t(fit$lambda.co) %*% W), unname(t(fit$lambda.tr)),
                 tolerance = 1e-10, info = m)
    ## no fixed effects, no covariates: the controls' outcomes weighted by W
    ## rebuild the treated counterfactuals (412d7ae misses by 111)
    mu <- if (is.null(fit$mu)) 0 else as.numeric(fit$mu)
    Y0.tr <- (fit$Y.dat - fit$eff)[, tr, drop = FALSE]
    rebuilt <- mu + (fit$Y.dat[, co, drop = FALSE] - mu) %*% W
    expect_equal(unname(rebuilt), unname(Y0.tr), tolerance = 1e-8, info = m)
  }
})

test_that("B5: simplex wgt.implied is Nco x Ntr with columns on the simplex", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  fit <- .fb_quiet(fect::fect(Y ~ D, data = simgsynth, index = c("id", "time"),
                              method = "ife", time.component.from = "nevertreated",
                              force = "two-way", r = 2, CV = FALSE, se = FALSE,
                              parallel = FALSE, loading.bound = "simplex",
                              gamma.loading = 1))
  expect_identical(fit$loading.bound, "simplex")
  W <- fit$wgt.implied
  expect_equal(dim(W), c(45L, 5L))   # 412d7ae: 5 x 45
  expect_equal(unname(colSums(W)), rep(1, 5), tolerance = 1e-6)
  expect_true(all(W >= -1e-10))
  expect_equal(unname(t(fit$lambda.co) %*% W), unname(t(fit$lambda.tr)),
               tolerance = 1e-8)
})


## -- B6  plots with se = FALSE and with a non-numeric time index ------------

.fb_builds <- function(p) {
  inherits(p, "ggplot") &&
    !inherits(tryCatch(ggplot2::ggplot_build(p), error = function(e) e), "error")
}

test_that("B6: counterfactual plots work for fits without SEs", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  fit <- .fb_quiet(fect::fect(Y ~ D + X1 + X2, data = simgsynth,
                              index = c("id", "time"), method = "gsynth",
                              force = "two-way", r = 2, CV = FALSE,
                              se = FALSE, parallel = FALSE))
  expect_null(fit[["vartype"]])   # se = FALSE fits store no vartype
  ## 412d7ae: "argument is of length zero"
  expect_true(.fb_builds(.fb_quiet(plot(fit, type = "counterfactual"))))
  expect_true(.fb_builds(.fb_quiet(plot(fit, type = "ct", id = 101))))
  turnout <- .fb_data("turnout")
  fit2 <- .fb_quiet(fect::fect(turnout ~ policy_edr + policy_mail_in + policy_motor,
                               data = turnout, index = c("abb", "year"),
                               method = "gsynth", force = "two-way", r = 1,
                               CV = FALSE, se = FALSE, parallel = FALSE))
  expect_true(.fb_builds(.fb_quiet(plot(fit2, type = "ct"))))
  expect_true(.fb_builds(.fb_quiet(plot(fit2, type = "ct", id = "CT"))))
})

test_that("B6: the factors plot works with a Date time index", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  simgsynth$date <- as.Date("2020-01-01") + 7 * (simgsynth$time - 1)
  fit <- .fb_quiet(fect::fect(Y ~ D + X1 + X2, data = simgsynth,
                              index = c("id", "date"), method = "gsynth",
                              force = "two-way", r = 2, CV = FALSE,
                              se = FALSE, parallel = FALSE))
  expect_s3_class(fit$rawtime, "Date")
  ## 412d7ae: "object 'T.b' not found"
  p <- .fb_quiet(plot(fit, type = "factors"))
  expect_true(.fb_builds(p))
  labs <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x$get_labels()
  expect_true("2020-01-01" %in% labs)
})


## -- B7  dropping periods without controls must leave treated cells --------

## 50 units, 62 months; 10 units treated from 2020-01. `all_post = TRUE`: the
## controls are unobserved in every post-treatment month; otherwise in two.
.fb_panel_057 <- function(all_post) {
  months <- as.integer(format(seq(as.Date("2015-08-01"), as.Date("2020-09-01"),
                                  by = "month"), "%Y%m"))
  TT <- length(months)
  set.seed(57)
  N <- 50; Ntr <- 10; start <- which(months == 202001)
  d <- expand.grid(t = 1:TT, u = 1:N)
  d$ID <- sprintf("u%02d", d$u)
  d$time <- months[d$t]
  d$D <- as.numeric(d$u <= Ntr & d$t >= start)
  f <- cumsum(stats::rnorm(TT)); l <- stats::rnorm(N)
  d$Y <- 20 + 3 * stats::rnorm(N)[d$u] + f[d$t] * l[d$u] - 4 * d$D +
    stats::rnorm(nrow(d))
  gone <- if (all_post) d$t >= start else d$t %in% c(start + 2, start + 5)
  d[!(gone & d$u > Ntr), c("ID", "time", "Y", "D")]
}

test_that("B7: no treated cell left after dropping control-free periods stops clearly", {
  skip_on_cran()
  d <- .fb_panel_057(all_post = TRUE)
  for (m in c("gsynth", "fe")) {
    ## 412d7ae: "non-numeric argument to binary operator" (gsynth),
    ## "missing value where TRUE/FALSE needed" (fe)
    expect_error(
      .fb_quiet(fect::fect(Y ~ D, data = d, index = c("ID", "time"), method = m,
                           force = "two-way", r = 1, CV = FALSE, se = FALSE,
                           parallel = FALSE)),
      "No treated observations remain after dropping the periods", info = m)
  }
  ## partial drops keep running
  dp <- .fb_panel_057(all_post = FALSE)
  fit <- .fb_quiet(fect::fect(Y ~ D, data = dp, index = c("ID", "time"),
                              method = "gsynth", force = "two-way", r = 1,
                              CV = FALSE, se = FALSE, parallel = FALSE))
  expect_true(is.finite(fit$att.avg))
  expect_equal(length(fit$rawtime), 60)   # two months dropped
})


## -- B8  loading.bound = "simplex" is honoured on every path ---------------

.fb_simplex_cols <- function(W) {
  is.matrix(W) && isTRUE(all.equal(unname(colSums(W)), rep(1, ncol(W)),
                                   tolerance = 1e-6)) && all(W >= -1e-10)
}

test_that("B8: simplex is kept for gsynth without SEs and in every CV path", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  run <- function(...) .fb_quiet(fect::fect(
    Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
    force = "two-way", se = FALSE, parallel = FALSE, seed = 1,
    loading.bound = "simplex", ...))
  fits <- list(
    gsynth_noCV = run(method = "gsynth", r = 2, CV = FALSE),
    gsynth_CV   = run(method = "gsynth", r = c(2, 2), CV = TRUE),
    ife_nt_CV   = run(method = "ife", time.component.from = "nevertreated",
                      r = c(2, 2), CV = TRUE))
  for (nm in names(fits)) {
    ## 412d7ae: loading.bound "none" and unbounded weights
    expect_identical(fits[[nm]]$loading.bound, "simplex", info = nm)
    expect_true(.fb_simplex_cols(fits[[nm]]$wgt.implied), info = nm)
  }
})

test_that("B8: every refit of a simplex fit with SEs uses the bound and the fit's gamma", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  orig <- get("fect_nevertreated", envir = asNamespace("fect"))
  rec <- new.env()
  settings <- list(
    list(CV = TRUE, r = c(2, 2), vartype = "bootstrap", gamma = NULL),
    list(CV = TRUE, r = c(2, 2), vartype = "parametric", gamma = NULL),
    list(CV = FALSE, r = 2, vartype = "parametric", gamma = 1))
  for (s in settings) {
    rec$lb <- character(0)
    rec$g <- numeric(0)
    fit <- testthat::with_mocked_bindings(
      .fb_quiet(fect::fect(Y ~ D, data = simgsynth, index = c("id", "time"),
                           method = "gsynth", force = "two-way", r = s$r,
                           CV = s$CV, se = TRUE, vartype = s$vartype,
                           nboots = 10, loading.bound = "simplex",
                           gamma.loading = s$gamma, parallel = FALSE,
                           seed = 1)),
      fect_nevertreated = function(..., loading.bound = "none",
                                   gamma.loading = NULL) {
        rec$lb <- c(rec$lb, loading.bound)
        rec$g <- c(rec$g, if (is.null(gamma.loading)) NA_real_ else gamma.loading)
        orig(..., loading.bound = loading.bound, gamma.loading = gamma.loading)
      },
      .package = "fect")
    info <- paste(s$vartype, "CV =", s$CV)
    expect_identical(fit$loading.bound, "simplex", info = info)
    ## the main fit plus 10 (bootstrap) or 20 (parametric) refits
    expect_gte(length(rec$lb), 11)
    ## 412d7ae: "none" in the CV main fit and in the parametric refits
    expect_true(all(rec$lb == "simplex"), info = info)
    ## the refits reuse the main fit's gamma (no CV of gamma per replicate)
    expect_true(all(rec$g[-1] == fit$gamma.loading), info = info)
    expect_true(all(is.finite(fit$att.avg.boot)), info = info)
  }
})

test_that("B8: simplex with no factor in the model says it has no effect", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  out <- .fb_msgs(fect::fect(Y ~ D, data = simgsynth, index = c("id", "time"),
                             method = "fe", time.component.from = "nevertreated",
                             se = FALSE, parallel = FALSE,
                             loading.bound = "simplex"))
  expect_false(inherits(out$res, "error"))
  expect_true(any(grepl("loading.bound = \"simplex\" has no effect because the selected number of factors is 0",
                        out$msgs, fixed = TRUE)))
  out0 <- .fb_msgs(fect::fect(Y ~ D, data = simgsynth, index = c("id", "time"),
                              method = "gsynth", r = 0, CV = FALSE, se = FALSE,
                              parallel = FALSE))
  expect_false(any(grepl("has no effect", out0$msgs, fixed = TRUE)))
})


## -- B9  CV folds do not depend on se; remove.id whenever units are removed --

test_that("B9: with parallel = FALSE, se = TRUE draws the CV folds of se = FALSE", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  for (m in c("gsynth", "ife")) {
    run <- function(se) .fb_quiet(fect::fect(
      Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"), method = m,
      force = "two-way", CV = TRUE, r = c(0, 4), se = se, nboots = 5,
      parallel = FALSE, seed = 12))
    f0 <- run(FALSE)
    f1 <- run(TRUE)
    slot <- if (m == "gsynth") "CV.out" else "CV.out.ife"
    ## 412d7ae: se = TRUE used the folds of seed + 1 (CV.out differs by 1.25)
    expect_identical(f1[[slot]], f0[[slot]], info = m)
    expect_identical(f1$r.cv, f0$r.cv, info = m)
    expect_identical(f1$att.avg, f0$att.avg, info = m)
  }
})

test_that("B9: remove.id lists the removed units even when the first unit is kept", {
  simgsynth <- .fb_data("simgsynth")
  ## unit 103 keeps only 3 pre-treatment periods (< min.T0 = 5)
  d2 <- simgsynth[!(simgsynth$id == 103 & simgsynth$time <= 17), ]
  fit <- .fb_quiet(fect::fect(Y ~ D, data = d2, index = c("id", "time"),
                              method = "gsynth", force = "two-way", r = 0,
                              CV = FALSE, se = FALSE, parallel = FALSE))
  expect_false(103 %in% fit$id)
  ## 412d7ae: remove.id was set only when unit 1 (here 101) was removed
  expect_equal(as.numeric(fit[["remove.id"]]), 103)
})
