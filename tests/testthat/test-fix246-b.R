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


## -- B10  se = TRUE uses the never-treated estimator for ife / fe / cfe ------

.fb_nt_fit <- function(...) {
  simgsynth <- .fb_data("simgsynth")
  .fb_quiet(fect::fect(Y ~ D + X1 + X2, data = simgsynth,
                       index = c("id", "time"), force = "two-way",
                       CV = FALSE, parallel = FALSE, seed = 1, nboots = 10,
                       ...))
}

test_that("B10: ife + never-treated with se = TRUE is the gsynth fit (all vartypes)", {
  skip_on_cran()
  f0 <- .fb_nt_fit(method = "ife", time.component.from = "nevertreated",
                   r = 2, se = FALSE)
  for (vt in c("bootstrap", "parametric", "jackknife")) {
    a <- .fb_nt_fit(method = "ife", time.component.from = "nevertreated",
                    r = 2, se = TRUE, vartype = vt)
    g <- .fb_nt_fit(method = "gsynth", r = 2, se = TRUE, vartype = vt)
    ## 412d7ae: 5.5788 (the not-yet-treated fit) vs 5.5433 without SEs
    expect_identical(a$att.avg, f0$att.avg, info = vt)
    for (s in c("att.avg", "att", "est.avg", "est.att", "est.beta",
                "att.avg.boot")) {
      expect_identical(a[[s]], g[[s]], info = paste(vt, s))
    }
    expect_identical(a$method, "gsynth", info = vt)
  }
  ## method = "fe" arrives as "ife" with r = 0
  fe <- .fb_nt_fit(method = "fe", time.component.from = "nevertreated",
                   se = TRUE, vartype = "bootstrap")
  g0 <- .fb_nt_fit(method = "gsynth", r = 0, se = TRUE, vartype = "bootstrap")
  for (s in c("att.avg", "att", "est.avg", "est.att", "att.avg.boot")) {
    expect_identical(fe[[s]], g0[[s]], info = paste("fe", s))
  }
})

test_that("B10: cfe + never-treated replicates use the never-treated estimator", {
  skip_on_cran()
  c0 <- .fb_nt_fit(method = "cfe", time.component.from = "nevertreated",
                   r = 2, se = FALSE)
  for (vt in c("bootstrap", "jackknife")) {
    cc <- .fb_nt_fit(method = "cfe", time.component.from = "nevertreated",
                     r = 2, se = TRUE, vartype = vt)
    g <- .fb_nt_fit(method = "gsynth", r = 2, se = TRUE, vartype = vt)
    expect_identical(cc$att.avg, c0$att.avg, info = vt)
    ## same model as gsynth (no extra FE); the solvers agree to ~5e-5.
    ## 412d7ae: the replicates were not-yet-treated CFE fits (off by 0.1-0.2)
    expect_lt(max(abs(cc$est.att[, "S.E."] - g$est.att[, "S.E."]), na.rm = TRUE),
              1e-3)
    expect_lt(max(abs(cc$att.avg.boot - g$att.avg.boot)), 1e-3)
  }
})

test_that("B10: the leave-one-period-out refits use the never-treated estimator", {
  skip_on_cran()
  lo <- function(...) .fb_nt_fit(r = 2, se = TRUE, loo = TRUE, ...)
  g  <- lo(method = "gsynth")
  a  <- lo(method = "ife", time.component.from = "nevertreated")
  cc <- lo(method = "cfe", time.component.from = "nevertreated")
  expect_identical(a$pre.est.att, g$pre.est.att)
  ## 412d7ae: the cfe loo refits were not-yet-treated CFE fits
  expect_lt(max(abs(cc$pre.est.att - g$pre.est.att), na.rm = TRUE), 1e-3)
})

test_that("B10: dloo with never-treated fixed effects stops before the bootstrap", {
  set.seed(42)
  rows <- list(); id <- 0
  for (g in c(4, 6, Inf)) for (u in seq_len(25)) {
    id <- id + 1
    y <- stats::rnorm(1) + 0.1 * (1:8) + stats::rnorm(8, 0, 0.4)
    d <- as.integer(!is.infinite(g) & (1:8) >= g)
    rows[[length(rows) + 1]] <- data.frame(id = id, time = 1:8, Y = y + 2 * d, D = d)
  }
  dat <- do.call(rbind, rows)
  ## 412d7ae: ran, silently ignoring time.component.from
  expect_error(
    suppressMessages(fect::fect(Y ~ D, data = dat, index = c("id", "time"),
                                method = "fe", force = "two-way", dloo = TRUE,
                                se = TRUE, vartype = "bootstrap", nboots = 10,
                                time.component.from = "nevertreated",
                                parallel = FALSE)),
    "time.component.from")
})


## -- B11  cumulative ATT = running sum of the per-period ATTs ----------------

test_that("B11: effect(), att.cumu() and estimand('att.cumu') return the running sum", {
  skip_on_cran()
  turnout <- .fb_data("turnout")
  fit <- .fb_quiet(fect::fect(
    turnout ~ policy_edr + policy_mail_in + policy_motor, data = turnout,
    index = c("abb", "year"), method = "gsynth", r = 0, CV = FALSE,
    force = "two-way", min.T0 = 5, se = TRUE, vartype = "bootstrap",
    nboots = 50, seed = 2139, keep.sims = TRUE, parallel = FALSE))
  k   <- 1:10
  att <- fit$att[match(k, fit$time)]
  n   <- fit$count[match(k, fit$time)]
  running <- cumsum(att)
  ## the treated counts vary over event time (9 8 6 6 6 3 3 3 3 3), so the
  ## two definitions differ: 21.501 vs 8.259 at k = 10
  eff <- .fb_quiet(fect::effect(fit, period = c(1, 10)))
  acu <- .fb_quiet(fect::att.cumu(fit, period = c(1, 10)))
  est <- .fb_quiet(fect::estimand(fit, "att.cumu", "event.time"))
  ov5 <- .fb_quiet(fect::estimand(fit, "att.cumu", "overall", window = c(1, 5)))
  expect_equal(unname(eff$effect.est.avg), running, tolerance = 1e-10)
  expect_equal(unname(acu[, 3]), running, tolerance = 1e-10)
  expect_equal(est$estimate[1:10], running, tolerance = 1e-10)
  expect_equal(ov5$estimate, running[5], tolerance = 1e-10)
  ## weighted = TRUE keeps the old count-weighted number
  acw <- .fb_quiet(fect::att.cumu(fit, period = c(1, 10), weighted = TRUE))
  expect_equal(unname(acw[, 3]), k * cumsum(att * n) / cumsum(n),
               tolerance = 1e-10)
  ## att.cumu() and effect() share the replicate running sums
  expect_equal(unname(acu[2:10, "S.E."]),
               unname(eff$effect.est.att[2:10, "S.E."]), tolerance = 1e-10)
  expect_error(.fb_quiet(fect::att.cumu(fit, period = c(1, 10), weighted = NA)),
               "must be TRUE or FALSE")
})

test_that("B11: equal counts per event time give the old numbers; no-SE fits work", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  fit <- .fb_quiet(fect::fect(Y ~ D + X1 + X2, data = simgsynth,
                              index = c("id", "time"), method = "gsynth",
                              force = "two-way", r = 2, CV = FALSE, se = TRUE,
                              vartype = "bootstrap", nboots = 20, seed = 3,
                              keep.sims = TRUE, parallel = FALSE))
  k   <- 1:10
  att <- fit$att[match(k, fit$time)]
  n   <- fit$count[match(k, fit$time)]
  expect_true(all(n == 5))
  old <- k * cumsum(att * n) / cumsum(n)
  eff <- .fb_quiet(fect::effect(fit, period = c(1, 10)))
  acu <- .fb_quiet(fect::att.cumu(fit, period = c(1, 10)))
  expect_equal(unname(eff$effect.est.avg), old, tolerance = 1e-10)
  expect_equal(unname(acu[, 3]), old, tolerance = 1e-10)
  ## fits without SEs (412d7ae: "non-numeric argument to mathematical function")
  f0 <- .fb_quiet(fect::fect(Y ~ D + X1 + X2, data = simgsynth,
                             index = c("id", "time"), method = "gsynth",
                             force = "two-way", r = 2, CV = FALSE, se = FALSE,
                             parallel = FALSE))
  a0 <- .fb_quiet(fect::att.cumu(f0, period = c(1, 10)))
  a1 <- .fb_quiet(fect::att.cumu(f0, period = c(1, 10), weighted = TRUE))
  expect_equal(unname(a0[, 3]), cumsum(f0$att[match(k, f0$time)]),
               tolerance = 1e-10)
  expect_equal(unname(a1[, 3]), unname(a0[, 3]), tolerance = 1e-10)
})


## -- B12  cfe + never-treated with one treated unit -------------------------

test_that("B12: cfe + never-treated runs with a single treated unit and equals gsynth", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  d <- simgsynth[!(simgsynth$id %in% 102:105), ]   # unit 101 + 45 controls
  y1   <- d$Y[d$id == 101][order(d$time[d$id == 101])]
  cbar <- tapply(d$Y[d$id != 101], d$time[d$id != 101], mean)
  hand <- mean(y1[21:30] - cbar[21:30]) - mean(y1[1:20] - cbar[1:20])
  run <- function(fml, method, force = "two-way", ...) {
    args <- list(fml, data = d, index = c("id", "time"), method = method,
                 force = force, r = 0, CV = FALSE, parallel = FALSE, seed = 1,
                 ...)
    if (method == "cfe") args$time.component.from <- "nevertreated"
    .fb_quiet(do.call(fect::fect, args))
  }
  ## 412d7ae: "'x' must be an array of at least two dimensions"
  f  <- run(Y ~ D, "cfe", se = FALSE)
  g  <- run(Y ~ D, "gsynth", se = FALSE)
  expect_equal(f$att.avg, hand, tolerance = 1e-10)
  expect_equal(f$att.avg, g$att.avg, tolerance = 1e-10)
  fx <- run(Y ~ D + X1 + X2, "cfe", se = FALSE)
  gx <- run(Y ~ D + X1 + X2, "gsynth", se = FALSE)
  expect_equal(fx$att.avg, gx$att.avg, tolerance = 1e-8)
  fu <- run(Y ~ D, "cfe", force = "unit", se = FALSE)
  gu <- run(Y ~ D, "gsynth", force = "unit", se = FALSE)
  expect_equal(fu$att.avg, gu$att.avg, tolerance = 1e-8)
  ## with SEs (and B10's never-treated replicates): the gsynth SE
  fb <- run(Y ~ D, "cfe", se = TRUE, vartype = "bootstrap", nboots = 20)
  gb <- run(Y ~ D, "gsynth", se = TRUE, vartype = "bootstrap", nboots = 20)
  expect_true(is.finite(fb$est.avg[1, "S.E."]))
  expect_equal(fb$est.avg[1, "S.E."], gb$est.avg[1, "S.E."], tolerance = 1e-6)
})


## -- B13  imputed_outcomes() reports the aggregation weights used -----------

test_that("B13: imputed_outcomes() reports the weights of a weighted fit", {
  simgsynth <- .fb_data("simgsynth")
  simgsynth$w2 <- 0.5 + (as.integer(factor(simgsynth$id)) %% 7) / 4
  for (role in c("W", "W.agg", "W.est")) {
    args <- list(Y ~ D, data = simgsynth, index = c("id", "time"),
                 method = "gsynth", force = "two-way", r = 2, CV = FALSE,
                 se = FALSE, parallel = FALSE)
    args[[role]] <- "w2"
    fit <- .fb_quiet(do.call(fect::fect, args))
    po  <- fect::imputed_outcomes(fit)
    w_data <- simgsynth$w2[match(paste(po$id, po$time),
                                 paste(simgsynth$id, simgsynth$time))]
    ## W.est weights the fit only; the aggregation stays unweighted
    expected <- if (role == "W.est") rep(1, nrow(po)) else w_data
    ## 412d7ae: fit$W.agg partially matched W.agg.col, giving W.agg = NA
    expect_equal(po$W.agg, expected, info = role)
    ## the reported weights are the ones the fit's ATT used
    expect_equal(sum(po$eff * po$W.agg) / sum(po$W.agg), fit$att.avg,
                 tolerance = 1e-10, info = role)
    if (role != "W.est") {
      expect_equal(dim(fit[["W.agg"]]), dim(fit$Y.dat), info = role)
      ## a fit object from before 2.4.6 has no W.agg slot
      old <- fit
      old[["W.agg"]] <- NULL
      expect_equal(fect::imputed_outcomes(old)$W.agg, expected, info = role)
    }
  }
})


## -- B8b  leave-one-out refits keep loading.bound (leader-added) ------------

test_that("B8b: the leave-one-period-out refits of a simplex fit use the bound", {
  skip_on_cran()
  simgsynth <- .fb_data("simgsynth")
  orig <- get("fect_boot", envir = asNamespace("fect"))
  rec <- new.env()
  loo_fit <- function(lb) {
    rec$lb <- character(0)
    rec$g <- numeric(0)
    testthat::with_mocked_bindings(
      .fb_quiet(fect::fect(Y ~ D, data = simgsynth, index = c("id", "time"),
                           method = "gsynth", force = "two-way", r = 2,
                           CV = FALSE, se = TRUE, nboots = 5, loo = TRUE,
                           loading.bound = lb,
                           gamma.loading = if (lb == "simplex") 1 else NULL,
                           parallel = FALSE, seed = 1)),
      fect_boot = function(..., loading.bound = "none", gamma.loading = NULL) {
        rec$lb <- c(rec$lb, loading.bound)
        rec$g <- c(rec$g, if (is.null(gamma.loading)) NA_real_ else gamma.loading)
        orig(..., loading.bound = loading.bound, gamma.loading = gamma.loading)
      },
      .package = "fect")
  }
  fs <- loo_fit("simplex")
  lb_s <- rec$lb
  g_s <- rec$g
  fn <- loo_fit("none")
  lb_n <- rec$lb
  ## the main fit plus one refit per pre-treatment period
  expect_gt(length(lb_s), 10)
  ## da20c2f / 412d7ae: every loo refit got loading.bound = "none"
  expect_true(all(lb_s == "simplex"))
  expect_true(all(g_s == 1))
  ## so the loo placebo estimates now differ from the unbounded fit's (before,
  ## the two fits shared the same unbounded refits and were identical)
  expect_identical(fs$loading.bound, "simplex")
  expect_gt(max(abs(fs$pre.est.att[, "ATT"] - fn$pre.est.att[, "ATT"])), 1e-4)
  ## a fit without the bound passes "none" to every refit, as before
  expect_true(all(lb_n == "none"))
  expect_true(all(is.finite(fs$pre.est.att[, "ATT"])))
})
