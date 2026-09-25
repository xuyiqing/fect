## ---------------------------------------------------------------
## fect 2.4.6 correctness fixes, group B (CV, weights in CV, implied
## weights, plots, guards). Run 2026-09-24-fix246-correctness; one block
## (or more) per item (B1-B9). Every test here fails on dev @ 412d7ae
## and passes after the fix.
## ---------------------------------------------------------------

.fix246b_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## Collect the messages of `expr` (warnings are muffled); errors propagate.
.fix246b_messages <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(
    suppressWarnings(expr),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(value = val, messages = msgs)
}

## Balanced panel with `Nco` never-treated units and `Ntr` units treated
## from period T0 + 1 on; two factors, two covariates.
.fix246b_panel <- function(Nco, Ntr = 3, TT = 30, T0 = 20, seed = 32) {
  set.seed(seed)
  N <- Ntr + Nco
  d <- expand.grid(time = seq_len(TT), id = seq_len(N))[, c("id", "time")]
  F <- matrix(stats::rnorm(TT * 2), TT, 2)
  L <- matrix(stats::rnorm(N * 2), N, 2)
  d$X1 <- stats::rnorm(nrow(d))
  d$X2 <- stats::rnorm(nrow(d))
  d$D <- as.numeric(d$id <= Ntr & d$time > T0)
  d$Y <- 1 + stats::rnorm(N)[d$id] + stats::rnorm(TT)[d$time] +
    rowSums(F[d$time, ] * L[d$id, ]) + 0.5 * d$X1 + 0.3 * d$X2 + 2 * d$D +
    stats::rnorm(nrow(d), sd = 0.3)
  d
}


## -- B1  criterion = "pc" -------------------------------------------------

test_that("B1: criterion = 'pc' picks the lowest-PC r for gsynth and never-treated CFE", {
  skip_on_cran()
  data(turnout, package = "fect")
  ## On turnout the PC is lowest at r = 4 while the MSPE rule picks r = 2.
  for (m in c("gsynth", "cfe")) {
    fit <- .fix246b_quiet(fect::fect(
      turnout ~ policy_edr + policy_mail_in + policy_motor, data = turnout,
      index = c("abb", "year"), method = m, time.component.from = "nevertreated",
      force = "two-way", CV = TRUE, r = c(0, 5), k = 5, criterion = "pc",
      se = FALSE, parallel = FALSE, seed = 1
    ))
    r.pc <- unname(fit$CV.out[which.min(fit$CV.out[, "PC"]), "r"])
    expect_equal(r.pc, 4, info = m)
    expect_equal(unname(fit$r.cv), r.pc, info = m)
    ## the final fit is estimated at the selected r
    expect_equal(ncol(fit$factor), r.pc, info = m)
  }
})

test_that("B1: fect's IFE CV table under criterion = 'pc' has the right columns", {
  skip_on_cran()
  data(turnout, package = "fect")
  fit <- .fix246b_quiet(fect::fect(
    turnout ~ policy_edr + policy_mail_in + policy_motor, data = turnout,
    index = c("abb", "year"), method = "ife", force = "two-way", CV = TRUE,
    r = c(0, 5), k = 5, criterion = "pc", se = FALSE, parallel = FALSE, seed = 1
  ))
  tab <- fit$CV.out.ife
  expect_identical(colnames(tab), c("r", "sigma2", "IC", "PC", "MSPTATT", "MSE"))
  ## every row holds its own statistics (no sentinel or NA left behind)
  expect_true(all(is.finite(tab)))
  expect_true(all(tab[, "PC"] < 1e19))
  expect_equal(unname(fit$r.cv), unname(tab[which.min(tab[, "PC"]), "r"]))
})


## -- B2  CV range for r capped by the number of never-treated units ------

test_that("B2: the CV range for r is capped at Nco - 1 with a message", {
  skip_on_cran()
  d <- .fix246b_panel(Nco = 4)
  for (m in c("gsynth", "cfe")) {
    ## r = c(0, 5) with 4 never-treated units crashed in panel_factor()
    ## ("Mat::head_cols(): size out of bounds")
    out <- .fix246b_messages(fect::fect(
      Y ~ D, data = d, index = c("id", "time"), method = m,
      time.component.from = "nevertreated", force = "two-way", CV = TRUE,
      r = c(0, 5), k = 5, se = FALSE, parallel = FALSE, seed = 1
    ))
    expect_true(any(grepl(
      "With 4 never-treated units at most 3 factor(s) can be estimated; cross-validation searches r = 0 to 3.",
      out$messages, fixed = TRUE)), info = m)
    expect_equal(unname(out$value$CV.out[, "r"]), 0:3, info = m)
    expect_true(out$value$r.cv <= 3, info = m)
    ## a feasible range is searched as given, without the message
    ok <- .fix246b_messages(fect::fect(
      Y ~ D, data = d, index = c("id", "time"), method = m,
      time.component.from = "nevertreated", force = "two-way", CV = TRUE,
      r = c(0, 3), k = 5, se = FALSE, parallel = FALSE, seed = 1
    ))
    expect_false(any(grepl("never-treated units at most", ok$messages)), info = m)
    expect_equal(unname(ok$value$CV.out[, "r"]), 0:3, info = m)
  }
})

test_that("B2: panel_factor() clamps r at min(T, N)", {
  set.seed(2)
  E <- matrix(stats::rnorm(20), 5, 4)   # T = 5 > N = 4
  big <- fect:::panel_factor(E, 6L)
  full <- fect:::panel_factor(E, 4L)
  expect_equal(dim(big$factor), c(5L, 4L))
  expect_equal(dim(big$lambda), c(4L, 4L))
  expect_equal(big$FE, full$FE)
  Et <- t(E)                             # T = 4 < N = 5
  bigt <- fect:::panel_factor(Et, 9L)
  expect_equal(dim(bigt$factor), c(4L, 4L))
  expect_equal(bigt$FE, fect:::panel_factor(Et, 4L)$FE)
  zero <- fect:::panel_factor(E, 0L)
  expect_equal(dim(zero$factor), c(5L, 0L))
  expect_equal(zero$FE, matrix(0, 5, 4))
})


## -- B3  the reason cross-validation is skipped -------------------------

test_that("B3: r = 0 with CV = TRUE says there is only one candidate r", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  for (m in c("gsynth", "cfe")) {
    out <- .fix246b_messages(fect::fect(
      Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"), method = m,
      time.component.from = "nevertreated", force = "two-way", CV = TRUE,
      r = 0, se = FALSE, parallel = FALSE
    ))
    expect_true(any(grepl(
      "Only one candidate number of factors (r = 0) was given, so cross-validation is skipped and r.cv = 0.",
      out$messages, fixed = TRUE)), info = m)
    expect_false(any(grepl("too few", out$messages)), info = m)
    expect_equal(unname(out$value$r.cv), 0, info = m)
  }
  ## too few pre-treatment records keep the old message: two-way FE with two
  ## pre-treatment periods leaves no room for a factor
  d <- .fix246b_panel(Nco = 10, TT = 10, T0 = 2)
  few <- .fix246b_messages(fect::fect(
    Y ~ D, data = d, index = c("id", "time"), method = "gsynth",
    force = "two-way", CV = TRUE, r = c(0, 2), min.T0 = 2, k = 5, se = FALSE,
    parallel = FALSE, seed = 1
  ))
  expect_true(any(grepl("pre-treatment records of treated units are too few",
                        few$messages, fixed = TRUE)))
  expect_false(any(grepl("Only one candidate", few$messages)))
})

test_that("B3: one never-treated unit skips CV with its own message", {
  skip_on_cran()
  ## Before the B2 cap this range crashed ("Not a matrix.").
  d <- .fix246b_panel(Nco = 1)
  out <- .fix246b_messages(fect::fect(
    Y ~ D, data = d, index = c("id", "time"), method = "gsynth",
    force = "two-way", CV = TRUE, r = c(0, 1), k = 5, se = FALSE,
    parallel = FALSE, seed = 1
  ))
  expect_true(any(grepl(
    "With 1 never-treated unit no factor can be estimated, so cross-validation is skipped and r.cv = 0.",
    out$messages, fixed = TRUE)), info = paste(out$messages, collapse = " | "))
  expect_false(any(grepl("too few", out$messages)))
  expect_equal(unname(out$value$r.cv), 0)
  expect_true(is.finite(out$value$att.avg))
})


## -- B4  W.agg stays out of the never-treated CV fits and scores ---------

## simgsynth with unit weights and 60 control cells dropped (unbalanced
## controls, so a fit weight would change the estimates).
.fix246b_weighted_panel <- function() {
  data(simgsynth, package = "fect")
  d <- simgsynth
  set.seed(7)
  ids <- unique(d$id)
  w <- stats::setNames(stats::runif(length(ids), 0.2, 5), ids)
  d$w <- w[as.character(d$id)]
  tr.ids <- unique(d$id[d$D == 1])
  ctrl.rows <- which(!(d$id %in% tr.ids))
  set.seed(3)
  d[-sample(ctrl.rows, 60), ]
}

test_that("B4: W.agg alone never enters the gsynth / never-treated CFE fit under CV", {
  skip_on_cran()
  d <- .fix246b_weighted_panel()
  for (m in c("gsynth", "cfe")) {
    fit <- function(...) .fix246b_quiet(fect::fect(
      Y ~ D + X1 + X2, data = d, index = c("id", "time"), method = m,
      time.component.from = "nevertreated", force = "two-way", k = 5,
      se = FALSE, parallel = FALSE, seed = 1, ...
    ))
    none <- fit(CV = TRUE, r = c(0, 3))
    agg  <- fit(CV = TRUE, r = c(0, 3), W.agg = "w")
    ## same fit, same CV table, same r
    expect_identical(agg$CV.out, none$CV.out, info = m)
    expect_identical(agg$r.cv, none$r.cv, info = m)
    expect_identical(agg$beta, none$beta, info = m)
    expect_identical(agg$eff, none$eff, info = m)
    ## the weights still enter the aggregation
    expect_false(isTRUE(all.equal(agg$att.avg, none$att.avg)), info = m)
    ## a fixed r without CV (CFE crashed here: "subscript out of bounds")
    none2 <- fit(CV = FALSE, r = 2)
    agg2  <- fit(CV = FALSE, r = 2, W.agg = "w")
    expect_identical(agg2$beta, none2$beta, info = m)
    expect_identical(agg2$eff, none2$eff, info = m)
    ## W.est still weights the fit, with or without CV
    est1 <- fit(CV = TRUE, r = c(2, 2), W.est = "w")
    est0 <- fit(CV = FALSE, r = 2, W.est = "w")
    expect_equal(est1$beta, est0$beta, tolerance = 1e-8, info = m)
    expect_false(isTRUE(all.equal(est0$beta, none2$beta)), info = m)
  }
})


## -- B5  wgt.implied: Xu (2017)'s formula, Nco x Ntr ----------------------

test_that("B5: wgt.implied rebuilds the counterfactual (gsynth and never-treated CFE)", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  for (m in c("gsynth", "cfe")) {
    fit <- .fix246b_quiet(fect::fect(
      Y ~ D, data = simgsynth, index = c("id", "time"), method = m,
      time.component.from = "nevertreated", force = "none", r = 2, CV = FALSE,
      se = FALSE, parallel = FALSE
    ))
    tr <- which(colSums(fit$D.dat) > 0)
    co <- which(colSums(fit$D.dat) == 0)
    W <- fit$wgt.implied
    expect_equal(dim(W), c(length(co), length(tr)), info = m)
    expect_identical(dimnames(W), list(as.character(fit$id[co]),
                                       as.character(fit$id[tr])), info = m)
    Lco <- unname(as.matrix(fit$lambda.co))
    Ltr <- unname(as.matrix(fit$lambda.tr))
    expect_equal(unname(t(Lco) %*% W), t(Ltr), tolerance = 1e-8, info = m)
    ## no covariates, no fixed effects: the weights rebuild Y(0) - mu exactly
    rebuilt <- (fit$Y.dat[, co] - fit$mu) %*% unname(W)
    expect_equal(unname(rebuilt), unname(fit$Y.ct[, tr] - fit$mu),
                 tolerance = 1e-8, info = m)
  }
})

test_that("B5: simplex weights have the same Nco x Ntr orientation", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  fit <- .fix246b_quiet(fect::fect(
    Y ~ D, data = simgsynth, index = c("id", "time"), method = "ife",
    time.component.from = "nevertreated", force = "two-way", r = 2,
    CV = FALSE, se = FALSE, parallel = FALSE, loading.bound = "simplex",
    gamma.loading = 1
  ))
  expect_identical(fit$loading.bound, "simplex")
  tr <- which(colSums(fit$D.dat) > 0)
  co <- which(colSums(fit$D.dat) == 0)
  W <- fit$wgt.implied
  expect_equal(dim(W), c(length(co), length(tr)))
  expect_identical(colnames(W), as.character(fit$id[tr]))
  expect_equal(unname(colSums(W)), rep(1, length(tr)), tolerance = 1e-6)
  expect_true(all(W >= -1e-10))
  for (i in seq_along(tr)) {
    expect_equal(as.numeric(fit$lambda.tr[i, ]),
                 as.numeric(crossprod(fit$lambda.co, W[, i])), tolerance = 1e-6)
  }
})


## -- B6  plots of fits without SEs, and factors with a Date index -------

test_that("B6: counterfactual plots work when se = FALSE (vartype is NULL)", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  fit <- .fix246b_quiet(fect::fect(
    Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
    method = "gsynth", force = "two-way", r = 2, CV = FALSE, se = FALSE,
    parallel = FALSE
  ))
  expect_null(fit$vartype)
  p1 <- .fix246b_quiet(plot(fit, type = "counterfactual"))
  expect_s3_class(p1, "ggplot")
  p2 <- .fix246b_quiet(plot(fit, type = "ct", id = 101))
  expect_s3_class(p2, "ggplot")
})

test_that("B6: the factors plot labels a Date time index", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  d <- simgsynth
  d$date <- as.Date("2000-01-01") + (d$time - 1) * 31
  fit <- .fix246b_quiet(fect::fect(
    Y ~ D + X1 + X2, data = d, index = c("id", "date"), method = "gsynth",
    force = "two-way", r = 2, CV = FALSE, se = FALSE, parallel = FALSE
  ))
  p <- .fix246b_quiet(plot(fit, type = "factors"))
  expect_s3_class(p, "ggplot")
  ## 30 periods: labels thinned to every 2nd period
  b <- ggplot2::ggplot_build(p)
  labs <- b$layout$panel_params[[1]]$x$get_labels()
  expect_equal(length(labs), 15)
  expect_equal(as.character(labs[1]), as.character(fit$rawtime[1]))
})


## -- B7  no treated observation left after dropping control-less periods --

test_that("B7: dropping every post-treatment period stops with a plain message", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  tr.ids <- unique(simgsynth$id[simgsynth$D == 1])
  start <- min(simgsynth$time[simgsynth$D == 1])
  ## the controls are not observed once treatment starts
  d <- simgsynth[!(simgsynth$time >= start & !(simgsynth$id %in% tr.ids)), ]
  for (m in c("gsynth", "ife", "fe", "mc")) {
    expect_error(
      .fix246b_quiet(fect::fect(
        Y ~ D, data = d, index = c("id", "time"), method = m,
        force = "two-way", r = 1, CV = FALSE,
        lambda = if (m == "mc") 0.1 else NULL, se = FALSE, parallel = FALSE
      )),
      "No treated observations remain after dropping the periods in which no unit is under control \\(21, 22, 23",
      info = m
    )
  }
  ## dropping only some post-treatment periods keeps running
  d2 <- simgsynth[!(simgsynth$time %in% c(start + 2, start + 3) &
                      !(simgsynth$id %in% tr.ids)), ]
  fit <- .fix246b_quiet(fect::fect(
    Y ~ D, data = d2, index = c("id", "time"), method = "gsynth",
    force = "two-way", r = 1, CV = FALSE, se = FALSE, parallel = FALSE
  ))
  expect_equal(nrow(fit$Y.dat), 28)
  expect_true(is.finite(fit$att.avg))
})


## -- B8  loading.bound = "simplex" honoured in every path ----------------

test_that("B8: simplex is honoured with CV and for method = 'gsynth' without SEs", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  fit <- function(...) .fix246b_quiet(fect::fect(
    Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
    force = "two-way", k = 5, parallel = FALSE, seed = 1,
    loading.bound = "simplex", gamma.loading = 1, ...
  ))
  ## the path that already honoured it: ife + nevertreated, fixed r
  ref <- fit(method = "ife", time.component.from = "nevertreated", r = 2,
             CV = FALSE, se = FALSE)
  expect_identical(ref$loading.bound, "simplex")
  cases <- list(
    gsynth_fixed = list(method = "gsynth", r = 2, CV = FALSE, se = FALSE),
    gsynth_cv    = list(method = "gsynth", r = c(0, 3), CV = TRUE, se = FALSE),
    ifent_cv     = list(method = "ife", time.component.from = "nevertreated",
                        r = c(0, 3), CV = TRUE, se = FALSE),
    gsynth_cv_se = list(method = "gsynth", r = c(0, 3), CV = TRUE, se = TRUE,
                        nboots = 10)
  )
  for (nm in names(cases)) {
    f <- do.call(fit, cases[[nm]])
    expect_identical(f$loading.bound, "simplex", info = nm)
    expect_equal(f$gamma.loading, 1, info = nm)
    expect_equal(unname(f$r.cv), 2, info = nm)
    expect_equal(f$att.avg, ref$att.avg, tolerance = 1e-8, info = nm)
    expect_equal(unname(colSums(f$wgt.implied)), rep(1, 5), tolerance = 1e-6,
                 info = nm)
  }
})

test_that("B8: parametric replicates bound the loadings too", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  ## impute_Y0() forwards the bound to fect_nevertreated()
  d <- simgsynth
  Y <- matrix(d$Y, 30, 50); D <- matrix(d$D, 30, 50)
  I <- matrix(1, 30, 50); II <- I; II[D == 1] <- 0
  T.on <- matrix(NA_real_, 30, 50); T.on[, 1:5] <- seq_len(30) - 20
  res <- .fix246b_quiet(fect:::impute_Y0(
    method = "gsynth", predictive = "nevertreated", Y = Y, X = NULL, D = D,
    W = NULL, I = I, II = II, T.on = T.on, tuning = 2, boot = 0, force = 3,
    hasRevs = 0,
    tol = 1e-5, loading.bound = "simplex", gamma.loading = 1
  ))
  expect_identical(res$loading.bound, "simplex")
  expect_equal(unname(colSums(res$wgt.implied)), rep(1, 5), tolerance = 1e-6)
  ## every fect_nevertreated() call of a parametric fit gets the bound: the
  ## main fit, the simulated errors (draw.error) and the replicates. trace()
  ## reads only loading.bound; a wrapper that evaluates list(...) would force
  ## arguments the function never uses (time.on.W does not exist without
  ## weights), and the replicates would fail inside try().
  seen <- new.env()
  seen$lb <- character(0)
  suppressMessages(trace(
    "fect_nevertreated", where = asNamespace("fect"), print = FALSE,
    tracer = bquote(assign("lb", c(get("lb", envir = .(seen)), loading.bound),
                           envir = .(seen)))
  ))
  withr::defer(suppressMessages(
    untrace("fect_nevertreated", where = asNamespace("fect"))
  ))
  fit <- .fix246b_quiet(fect::fect(
    Y ~ D, data = simgsynth, index = c("id", "time"), method = "gsynth",
    force = "two-way", r = 2, CV = FALSE, se = TRUE, vartype = "parametric",
    nboots = 5, parallel = FALSE, seed = 1, loading.bound = "simplex",
    gamma.loading = 1
  ))
  expect_identical(fit$loading.bound, "simplex")
  expect_gte(length(seen$lb), 11)   # 1 main fit + 5 error draws + 5 replicates
  expect_true(all(seen$lb == "simplex"))
})

test_that("B8: simplex with no factors says it has no effect", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  out <- .fix246b_messages(fect::fect(
    Y ~ D, data = simgsynth, index = c("id", "time"), method = "gsynth",
    force = "two-way", r = 0, CV = FALSE, se = FALSE, parallel = FALSE,
    loading.bound = "simplex", gamma.loading = 1
  ))
  expect_true(any(grepl(
    "loading.bound = \"simplex\" has no effect because the selected number of factors is 0.",
    out$messages, fixed = TRUE)))
})


## -- B9  CV folds do not depend on se; remove.id --------------------------

## N = 30, T = 16, 10 treated units adopting at t = 11..13, two factors, one
## covariate (the panel of test-cv-seed.R).
.fix246b_seed_panel <- function() {
  withr::with_seed(2, {
    N <- 30; TT <- 16; Ntr <- 10
    F <- matrix(stats::rnorm(TT * 2), TT, 2)
    L <- matrix(stats::rnorm(N * 2), N, 2); L[, 2] <- 0.5 * L[, 2]
    a <- stats::rnorm(N); xi <- stats::rnorm(TT)
    X <- matrix(stats::rnorm(N * TT), TT, N)
    T0 <- rep(Inf, N); T0[seq_len(Ntr)] <- sample(11:13, Ntr, replace = TRUE)
    D <- sapply(seq_len(N), function(i) as.integer(seq_len(TT) >= T0[i]))
    Y <- outer(xi, rep(1, N)) + outer(rep(1, TT), a) + F %*% t(L) + 0.5 * X +
      2 * D + matrix(stats::rnorm(N * TT), TT, N)
    data.frame(id = rep(seq_len(N), each = TT), time = rep(seq_len(TT), N),
               Y = c(Y), D = c(D), X = c(X))
  })
}

test_that("B9: with parallel = FALSE, se = TRUE draws the same CV folds as se = FALSE", {
  skip_on_cran()
  d <- .fix246b_seed_panel()
  for (m in c("ife", "gsynth")) {
    tab <- if (m == "ife") "CV.out.ife" else "CV.out"
    fit <- function(...) .fix246b_quiet(fect::fect(
      Y ~ D + X, data = d, index = c("id", "time"), method = m, CV = TRUE,
      r = c(0, 3), k = 5, seed = 7, parallel = FALSE, ...
    ))
    nose <- fit(se = FALSE)
    boot <- fit(se = TRUE, nboots = 5)
    expect_identical(boot[[tab]], nose[[tab]], info = m)
    expect_identical(boot$r.cv, nose$r.cv, info = m)
    expect_identical(boot$att.avg, nose$att.avg, info = m)
  }
})

test_that("B9: remove.id lists the removed units even when unit 1 stays", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  ## control unit 120 is observed in 3 periods only (< min.T0 = 5)
  d <- simgsynth[!(simgsynth$id == 120 & simgsynth$time > 3), ]
  fit <- .fix246b_quiet(fect::fect(
    Y ~ D, data = d, index = c("id", "time"), method = "gsynth", r = 1,
    CV = FALSE, se = FALSE, parallel = FALSE
  ))
  expect_equal(ncol(fit$Y.dat), 49)
  expect_equal(fit$remove.id, 120)
})
