## vartype = "conformal": pooled information-matched Jackknife+ (Liu, v2).
## See R/conformal.R. Small panels; every fit runs with parallel = FALSE.

.cf_sim <- function(J = 20, Ntr = 1, TT = 20, T0 = 15, r = 1, seed = 1,
                    adopt = NULL, tau = 1) {
  set.seed(seed)
  N <- J + Ntr
  F <- matrix(rnorm(TT * r), TT, r)
  L <- matrix(rnorm(N * r), N, r)
  Y0 <- outer(rnorm(TT), rep(1, N)) + outer(rep(1, TT), rnorm(N)) +
    F %*% t(L) + matrix(rnorm(TT * N), TT, N)
  if (is.null(adopt)) adopt <- rep(T0 + 1, Ntr)
  D <- matrix(0, TT, N)
  for (i in seq_len(Ntr)) D[seq_len(TT) >= adopt[i], i] <- 1
  Y <- Y0 + tau * D
  data.frame(id = rep(seq_len(N), each = TT), time = rep(seq_len(TT), N),
             Y = c(Y), D = c(D))
}

.cf_fit <- function(df, ...) {
  suppressMessages(fect(Y ~ D, data = df, index = c("id", "time"),
                        se = TRUE, vartype = "conformal", parallel = FALSE,
                        CV = FALSE, ...))
}

## ---- calibration ---------------------------------------------------------

test_that("q_gamma matches hand-computed values", {
  q <- fect:::.conformal_q
  ## (2q - 1)_+ / N <= gamma, q <= floor(N / 2)
  expect_equal(q(5, 0.05), 0L)   # 2q - 1 <= 0.25
  expect_equal(q(5, 0.10), 0L)   # 2q - 1 <= 0.5
  expect_equal(q(10, 0.05), 0L)  # 2q - 1 <= 0.5
  expect_equal(q(10, 0.10), 1L)  # 2q - 1 <= 1
  expect_equal(q(20, 0.05), 1L)  # 2q - 1 <= 1
  expect_equal(q(20, 0.10), 1L)  # 2q - 1 <= 2
  expect_equal(q(40, 0.05), 1L)  # 2q - 1 <= 2
  expect_equal(q(40, 0.10), 2L)  # 2q - 1 <= 4
  expect_equal(q(31, 0.10), 2L)  # 2q - 1 <= 3.1
  expect_equal(q(61, 0.10), 3L)  # 2q - 1 <= 6.1
  expect_equal(q(4, 0.99), 2L)   # capped at floor(N / 2)
})

test_that("the interval is the real line when gamma < 1/N", {
  set.seed(2)
  M0 <- rnorm(9)
  R <- abs(rnorm(9))
  v <- fect:::.conformal_interval(M0, R, y = 0.3, gamma = 0.05)  # N = 10
  expect_equal(unname(v["q"]), 0)
  expect_equal(unname(v[c("CI.lower", "CI.upper")]), c(-Inf, Inf))
  expect_true(v["p.value"] > 0.05)
})

test_that("interval endpoints are the order statistics and the p-value is dual", {
  set.seed(3)
  J <- 39
  M0 <- rnorm(J)
  R <- abs(rnorm(J))
  a <- sort(M0 - R)
  b <- sort(M0 + R)
  for (gam in c(0.05, 0.1, 0.2)) {
    v <- fect:::.conformal_interval(M0, R, y = 0, gamma = gam)
    q <- fect:::.conformal_q(J + 1, gam)
    expect_equal(unname(v[c("pred.lower", "pred.upper")]),
                 c(a[q], b[J + 1 - q]))
    expect_equal(unname(v[c("CI.lower", "CI.upper")]),
                 c(0 - b[J + 1 - q], 0 - a[q]))
  }
  ## duality: p(c = 0) <= gamma exactly when 0 is outside the effect interval
  for (y in seq(-4, 4, by = 0.05)) {
    for (gam in c(0.05, 0.1)) {
      v <- fect:::.conformal_interval(M0, R, y = y, gamma = gam)
      out <- 0 < v["CI.lower"] || 0 > v["CI.upper"]
      expect_identical(unname(v["p.value"] <= gam), unname(out))
    }
  }
})

## ---- gsynth folds against a manual reimplementation ----------------------

test_that("gsynth fold predictions equal factor estimation on B_j plus projection", {
  df <- .cf_sim(J = 12, TT = 20, T0 = 15, r = 1, seed = 11)
  fit <- .cf_fit(df, method = "gsynth", r = 1, keep.sims = TRUE,
                 tol = 1e-10)
  Y <- matrix(df$Y, 20, 13)
  co <- 2:13
  pre <- 1:15
  post <- 16:20
  manual_pred <- function(B, y) {
    YB <- Y[, B, drop = FALSE]
    mu <- mean(YB)
    xi <- rowMeans(YB) - mu
    al <- colMeans(YB) - mu
    Yd <- YB - mu - outer(xi, rep(1, ncol(YB))) - outer(rep(1, 20), al)
    Fh <- svd(Yd)$u[, 1, drop = FALSE]
    U <- y - mu - xi
    cf <- stats::lm.fit(cbind(1, Fh[pre, , drop = FALSE]), U[pre])$coefficients
    mu + xi + cbind(1, Fh) %*% cf
  }
  M0.man <- E.man <- matrix(NA_real_, length(post), length(co))
  for (jj in seq_along(co)) {
    B <- setdiff(co, co[jj])
    M0.man[, jj] <- manual_pred(B, Y[, 1])[post]
    E.man[, jj] <- Y[post, co[jj]] - manual_pred(B, Y[, co[jj]])[post]
  }
  folds <- fit$conformal$folds
  expect_equal(folds$cell$time, post)
  expect_equal(unname(folds$M0), M0.man, tolerance = 1e-6)
  expect_equal(unname(folds$E), E.man, tolerance = 1e-6)

  ## the reported date-specific intervals are built from these folds
  for (s in seq_along(post)) {
    v <- fect:::.conformal_interval(M0.man[s, ], abs(E.man[s, ]),
                                    Y[post[s], 1], 0.05)
    u <- fit$conformal$unit[s, ]
    expect_equal(c(u$CI.lower, u$CI.upper),
                 unname(v[c("CI.lower", "CI.upper")]), tolerance = 1e-6)
    expect_equal(u$N, 13)
  }
  ## one treated unit: the event-time interval is the unit's interval
  expect_equal(unname(fit$est.att[as.character(1:5), "CI.lower"]),
               fit$conformal$unit$CI.lower)
})

## ---- point estimates and outputs -----------------------------------------

test_that("point estimates equal the ordinary fit's", {
  df <- .cf_sim(J = 20, Ntr = 3, seed = 4)
  for (m in c("gsynth", "ife")) {
    cf <- .cf_fit(df, method = m, r = 1)
    pl <- suppressMessages(fect(Y ~ D, data = df, index = c("id", "time"),
                                method = m, r = 1, CV = FALSE, se = FALSE))
    expect_equal(cf$att, pl$att)
    expect_equal(cf$att.avg, pl$att.avg)
    expect_equal(cf$eff, pl$eff)
    expect_equal(unname(cf$est.att[, "ATT"]), unname(pl$att))
    expect_equal(unname(cf$est.avg[1, "ATT.avg"]), pl$att.avg)
    expect_true(all(is.na(cf$est.att[, "S.E."])))
  }
})

test_that("ife and mc give finite intervals", {
  df <- .cf_sim(J = 25, Ntr = 1, seed = 5)
  fi <- .cf_fit(df, method = "ife", r = 1)
  fm <- .cf_fit(df, method = "mc", lambda = 0.01)
  for (f in list(fi, fm)) {
    post <- f$est.att[f$time >= 1, , drop = FALSE]
    expect_true(all(is.finite(post[, "CI.lower"])))
    expect_true(all(is.finite(post[, "CI.upper"])))
    expect_true(all(post[, "CI.lower"] <= post[, "CI.upper"]))
    expect_true(all(is.finite(f$conformal$unit$CI.lower)))
    expect_equal(f$conformal$n.failed.fits, 0)
    expect_identical(f$conformal$status, "ok")
  }
})

test_that("staggered adoption runs and the average interval covers its estimate", {
  df <- .cf_sim(J = 30, Ntr = 6, TT = 20, seed = 6,
                adopt = c(13, 13, 15, 15, 17, 17))
  fit <- .cf_fit(df, method = "gsynth", r = 1)
  expect_true(is.finite(fit$est.avg[1, "CI.lower"]))
  expect_true(fit$est.avg[1, "CI.lower"] <= fit$est.avg[1, "ATT.avg"])
  expect_true(fit$est.avg[1, "ATT.avg"] <= fit$est.avg[1, "CI.upper"])
  expect_true(fit$est.avg.unit[1, "CI.lower"] <= fit$est.avg.unit[1, "ATT.avg.unit"])
  expect_true(fit$est.avg.unit[1, "ATT.avg.unit"] <= fit$est.avg.unit[1, "CI.upper"])
  ## every post-treatment event time has an interval, pre-treatment ones none
  post <- fit$time >= 1
  expect_true(all(is.finite(fit$est.att[post, "CI.lower"])))
  expect_true(all(is.na(fit$est.att[!post, "CI.lower"])))
  ## per-unit table: one row per treated post-treatment cell
  expect_equal(nrow(fit$conformal$unit), sum(df$D))
  expect_equal(fit$conformal$unit$eff,
               fit$eff[cbind(match(fit$conformal$unit$time, rownames(fit$eff)),
                             match(fit$conformal$unit$unit, colnames(fit$eff)))])
  ## the inner band (att.bound) is narrower than the main interval
  w <- fit$est.att[post, "CI.upper"] - fit$est.att[post, "CI.lower"]
  w90 <- fit$att.bound[post, "CI.upper"] - fit$att.bound[post, "CI.lower"]
  expect_true(all(w90 <= w + 1e-12))
  ## ife with staggered adoption
  fi <- .cf_fit(df, method = "ife", r = 1)
  expect_true(all(is.finite(fi$est.att[post, "CI.lower"])))
})

test_that("normalize = TRUE gives the same intervals", {
  df <- .cf_sim(J = 20, seed = 7)
  df$Y <- df$Y * 50
  f1 <- .cf_fit(df, method = "gsynth", r = 1, tol = 1e-10)
  f2 <- .cf_fit(df, method = "gsynth", r = 1, tol = 1e-10, normalize = TRUE)
  expect_equal(f1$est.att, f2$est.att, tolerance = 1e-6)
  expect_equal(f1$est.avg, f2$est.avg, tolerance = 1e-6)
  m1 <- .cf_fit(df, method = "mc", lambda = 0.5, tol = 1e-8)
  m2 <- .cf_fit(df, method = "mc", lambda = 0.5, tol = 1e-8, normalize = TRUE)
  expect_equal(m1$est.att, m2$est.att, tolerance = 1e-4)
})

test_that("controls missing the treated pre-period window are dropped and counted", {
  df <- .cf_sim(J = 25, seed = 8)
  ## two controls miss a pre-treatment date of the treated unit
  df <- df[!(df$id %in% c(5, 9) & df$time == 4), ]
  fit <- .cf_fit(df, method = "gsynth", r = 1)
  expect_equal(fit$conformal$n.controls, 25)
  expect_equal(fit$conformal$n.controls.ineligible, 2)
  expect_true(all(fit$conformal$unit$N == 24))
  ## a control missing one post-treatment date drops out at that date only
  df2 <- .cf_sim(J = 25, seed = 8)
  df2 <- df2[!(df2$id == 7 & df2$time == 18), ]
  fit2 <- .cf_fit(df2, method = "gsynth", r = 1)
  N2 <- fit2$conformal$unit$N
  expect_equal(N2[fit2$conformal$unit$time == 18], 25)
  expect_true(all(N2[fit2$conformal$unit$time != 18] == 26))
})

test_that("too few controls give the whole real line", {
  df <- .cf_sim(J = 8, seed = 9)  # N = 9 < 1 / 0.05
  fit <- .cf_fit(df, method = "gsynth", r = 1)
  post <- fit$time >= 1
  expect_true(all(fit$est.att[post, "CI.lower"] == -Inf))
  expect_true(all(fit$est.att[post, "CI.upper"] == Inf))
  expect_match(fit$conformal$status, "infinite")
})

## ---- unsupported options ----------------------------------------------------

test_that("unsupported methods and options stop with a message", {
  df <- .cf_sim(J = 15, seed = 10)
  expect_error(.cf_fit(df, method = "fe"), "method = \"ife\" with r = 0")
  expect_error(.cf_fit(df, method = "cfe", r = 0), "gsynth")
  expect_error(.cf_fit(df, method = "both"), "gsynth")
  df$w <- 1
  expect_error(.cf_fit(df, method = "gsynth", r = 1, W = "w"), "weights")
  ## reversals
  dr <- df
  dr$D[dr$id == 1 & dr$time == 20] <- 0
  expect_error(.cf_fit(dr, method = "ife", r = 1), "absorbing")
  ## no never-treated units (ife with not-yet-treated controls still fits)
  dn <- .cf_sim(J = 1, Ntr = 6, TT = 20, seed = 10,
                adopt = c(8, 9, 10, 11, 12, 13))
  expect_error(.cf_fit(dn, method = "ife", r = 0), "never-treated")
})

test_that("covariates and cross-validated r are held fixed in the folds", {
  df <- .cf_sim(J = 25, seed = 13)
  set.seed(13)
  df$x <- rnorm(nrow(df))
  df$Y <- df$Y + 0.5 * df$x
  for (m in c("gsynth", "ife")) {
    fit <- suppressMessages(fect(Y ~ D + x, data = df, index = c("id", "time"),
                                 method = m, r = 1, CV = FALSE, se = TRUE,
                                 vartype = "conformal", parallel = FALSE))
    expect_true(all(is.finite(fit$conformal$unit$CI.lower)))
    expect_true(is.na(fit$est.beta[1, "S.E."]))
    expect_equal(unname(fit$est.beta[1, "Coef"]), unname(fit$beta[1, 1]))
  }
  cv <- suppressMessages(fect(Y ~ D, data = df, index = c("id", "time"),
                              method = "gsynth", r = c(0, 3), CV = TRUE,
                              se = TRUE, vartype = "conformal",
                              parallel = FALSE, seed = 1))
  expect_equal(cv$conformal$r, as.numeric(cv$r.cv))
  expect_true(all(is.finite(cv$conformal$unit$CI.lower)))
})

## ---- print, plot and post-fit functions ---------------------------------------

test_that("print and plot handle conformal fits", {
  df <- .cf_sim(J = 25, Ntr = 2, seed = 12)
  fit <- .cf_fit(df, method = "gsynth", r = 1)
  out <- utils::capture.output(print(fit))
  expect_true(any(grepl("conformal", out)))
  pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_s3_class(plot(fit), "ggplot")
  expect_s3_class(plot(fit, type = "counterfactual"), "ggplot")
  expect_s3_class(suppressMessages(plot(fit, type = "calendar")), "ggplot")
  expect_s3_class(suppressMessages(plot(fit, type = "gap", id = 1)), "ggplot")
  expect_error(plot(fit, type = "equiv"), "conformal")
  expect_error(suppressMessages(estimand(fit, "att", "overall")), "conformal")
  expect_error(suppressMessages(att.cumu(fit)), "conformal")
})

## ---- rule D: failed folds in .conformal_interval ---------------------------

test_that("a failed fold keeps N and q and never shortens the interval", {
  set.seed(21)
  J <- 39
  M0 <- rnorm(J)
  R <- abs(rnorm(J))
  f5 <- seq_len(J) == 5
  for (gam in c(0.1, 0.2)) {
    base <- fect:::.conformal_interval(M0, R, y = 0.2, gamma = gam)
    one <- fect:::.conformal_interval(M0, R, y = 0.2, gamma = gam, failed = f5)
    expect_equal(unname(one[c("N", "q")]), unname(base[c("N", "q")]))
    expect_equal(unname(one["n.failed"]), 1)
    expect_true(one["pred.lower"] <= base["pred.lower"])
    expect_true(one["pred.upper"] >= base["pred.upper"])
    expect_true(one["CI.lower"] <= base["CI.lower"])
    expect_true(one["CI.upper"] >= base["CI.upper"])
    ## values behind a failed slot are ignored, even NA
    M0.na <- M0
    R.na <- R
    M0.na[5] <- NA
    R.na[5] <- NA
    one.na <- fect:::.conformal_interval(M0.na, R.na, y = 0.2, gamma = gam,
                                         failed = f5)
    expect_identical(one.na, one)
  }
})

test_that("q or more failed folds give the whole real line, fewer do not", {
  set.seed(22)
  J <- 39  # N = 40, gamma = 0.1 -> q = 2
  M0 <- rnorm(J)
  R <- abs(rnorm(J))
  expect_equal(fect:::.conformal_q(J + 1, 0.1), 2L)
  f1 <- seq_len(J) == 3
  f2 <- seq_len(J) %in% c(3, 30)
  v1 <- fect:::.conformal_interval(M0, R, y = 0, gamma = 0.1, failed = f1)
  v2 <- fect:::.conformal_interval(M0, R, y = 0, gamma = 0.1, failed = f2)
  expect_true(is.finite(v1["pred.lower"]) && is.finite(v1["pred.upper"]))
  ## one -Inf in a: the q = 2 smallest value is the smallest finite a
  expect_equal(unname(v1["pred.lower"]), sort((M0 - R)[!f1])[1])
  ## two failed slots: the q-th smallest a is -Inf, the (N - q)-th b is +Inf
  expect_equal(unname(v2[c("pred.lower", "pred.upper")]), c(-Inf, Inf))
  expect_equal(unname(v2[c("CI.lower", "CI.upper")]), c(-Inf, Inf))
  expect_equal(unname(v2[c("N", "q", "n.failed")]), c(40, 2, 2))
})

test_that("failed slots keep the p-value dual to the interval", {
  set.seed(23)
  J <- 39
  M0 <- rnorm(J)
  R <- abs(rnorm(J))
  fl <- seq_len(J) == 11
  for (y in seq(-4, 4, by = 0.05)) {
    for (gam in c(0.05, 0.1, 0.2)) {
      v <- fect:::.conformal_interval(M0, R, y = y, gamma = gam, failed = fl)
      out <- 0 < v["CI.lower"] || 0 > v["CI.upper"]
      expect_identical(unname(v["p.value"] <= gam), unname(out))
    }
  }
  ## a failed slot is no tail violation: three of four slots violate
  v <- fect:::.conformal_interval(c(0, 0, 0, 0), c(1, 1, 1, 1), y = 5,
                                  gamma = 0.5,
                                  failed = c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(unname(v["p.value"]), min(1, (2 * (5 - 3) - 1) / 5))
})

test_that("failed = NULL reproduces the no-failure result; pool size sets N", {
  set.seed(24)
  M0 <- rnorm(30)
  R <- abs(rnorm(30))
  v0 <- fect:::.conformal_interval(M0, R, y = 0.1, gamma = 0.1)
  v1 <- fect:::.conformal_interval(M0, R, y = 0.1, gamma = 0.1,
                                   failed = rep(FALSE, 30))
  expect_identical(v0, v1)
  expect_equal(unname(v0["n.failed"]), 0)
  expect_equal(names(v0), c("pred.lower", "pred.upper", "CI.lower",
                            "CI.upper", "p.value", "N", "q", "n.failed"))
  ## a control removed by the caller (out of pool) reduces N; a failed one
  ## does not
  v.drop <- fect:::.conformal_interval(M0[-1], R[-1], y = 0.1, gamma = 0.1)
  v.fail <- fect:::.conformal_interval(M0, R, y = 0.1, gamma = 0.1,
                                       failed = seq_len(30) == 1)
  expect_equal(unname(v.drop["N"]), 30)
  expect_equal(unname(v.fail["N"]), 31)
  ## non-finite inputs count as failed
  M0.bad <- M0
  M0.bad[2] <- NA
  expect_equal(unname(fect:::.conformal_interval(M0.bad, R, y = 0.1,
                                                 gamma = 0.1)["n.failed"]), 1)
  expect_error(fect:::.conformal_interval(M0, R, 0.1, 0.1, failed = TRUE),
               "same length")
  expect_error(fect:::.conformal_interval(M0, R[-1], 0.1, 0.1), "same length")
})
