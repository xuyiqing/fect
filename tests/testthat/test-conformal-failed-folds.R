## Rule D for failed folds in vartype = "conformal" (tester, from test-spec.md).
## Failure injection: mock fect:::.conformal_component_fit (single component fit).

.cff_sim <- function(J = 20, Ntr = 1, TT = 20, T0 = 15, r = 1, seed = 1,
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

.cff_fit <- function(df, ...) {
  suppressMessages(fect(Y ~ D, data = df, index = c("id", "time"),
                        se = TRUE, vartype = "conformal", parallel = FALSE,
                        CV = FALSE, ...))
}

## optional log of observed values (CFF_LOG=<file>)
.cff_log <- function(label, ...) {
  f <- Sys.getenv("CFF_LOG", "")
  if (nzchar(f)) {
    cat(label, ": ", paste(format(c(...), digits = 10), collapse = " | "),
        "\n", file = f, append = TRUE, sep = "")
  }
  invisible(NULL)
}

## fit with an injected failure; returns fit, collected warnings, call count
## fail.calls: component-fit call numbers that return NULL
## nan.calls / nan.fun: call numbers whose returned matrix is modified
## all.fail: every call returns NULL
.cff_run <- function(df, fail.calls = integer(0), nan.calls = integer(0),
                     nan.fun = NULL, all.fail = FALSE, ...) {
  orig <- fect:::.conformal_component_fit
  st <- new.env()
  st$n <- 0L
  mock <- function(...) {
    st$n <- st$n + 1L
    if (all.fail || st$n %in% fail.calls) return(NULL)
    m <- orig(...)
    if (st$n %in% nan.calls) m <- nan.fun(m)
    m
  }
  testthat::local_mocked_bindings(.conformal_component_fit = mock,
                                  .package = "fect")
  ws <- character(0)
  fit <- withCallingHandlers(
    .cff_fit(df, ...),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, ws = ws, calls = st$n)
}

.cff_nwarn <- function(r) sum(grepl("failed", r$ws))

## interval of `new` contains interval of `ref` (NA in ref skipped)
.cff_contains <- function(new, ref, lo = "CI.lower", hi = "CI.upper") {
  keep <- !is.na(ref[[lo]]) & !is.na(ref[[hi]])
  isTRUE(all(new[[lo]][keep] <= ref[[lo]][keep])) &&
    isTRUE(all(new[[hi]][keep] >= ref[[hi]][keep]))
}

.cff_post <- function(fit) {
  ea <- fit$est.att
  as.data.frame(ea[as.numeric(rownames(ea)) >= 1 & !is.na(ea[, "CI.lower"]), ,
                   drop = FALSE])
}

## section 6 invariants
.cff_invariants <- function(fit, alpha) {
  cf <- fit$conformal
  u <- cf$unit
  a <- cf$att[!is.na(cf$att$N), ]
  v <- cf$avg
  for (tab in list(u, a, v)) {
    expect_equal(tab$N, cf$n.controls - tab$n.excluded + 1)
    expect_true(all(tab$n.failed >= 0 & tab$n.failed <= tab$N - 1))
    expect_equal(tab$q, vapply(tab$N, fect:::.conformal_q, integer(1),
                               gamma = alpha) * 1)
  }
  expect_gte(cf$n.failed.fits, cf$n.failed.folds)
  expect_equal(cf$n.intervals.affected, sum(u$n.failed > 0))
}

## point estimates and pre-treatment rows are not touched by a failure
.cff_same_points <- function(fit, ref) {
  expect_identical(fit$est.att[, "ATT"], ref$est.att[, "ATT"])
  expect_identical(fit$est.avg[1, "ATT.avg"], ref$est.avg[1, "ATT.avg"])
  expect_identical(fit$eff, ref$eff)
  expect_identical(fit$att, ref$att)
  pre <- as.numeric(rownames(fit$est.att)) < 1
  expect_true(all(is.na(fit$est.att[pre, "CI.lower"])))
  expect_true(all(is.na(fit$est.att[pre, "CI.upper"])))
}

## ---- legacy e99f2d2 .conformal_interval (oracle) ---------------------------
.cff_legacy <- function(M0, R, y, gamma) {
  ok <- is.finite(M0) & is.finite(R)
  M0 <- M0[ok]
  R <- R[ok]
  J <- length(M0)
  N <- J + 1L
  q <- if (J >= 1L) fect:::.conformal_q(N, gamma) else 0L
  out <- c(pred.lower = -Inf, pred.upper = Inf,
           CI.lower = -Inf, CI.upper = Inf,
           p.value = 1, N = N, q = q)
  if (J < 1L || !is.finite(y)) {
    out[c("p.value")] <- NA_real_
    if (!is.finite(y)) out[c("CI.lower", "CI.upper")] <- NA_real_
    return(out)
  }
  a <- M0 - R
  b <- M0 + R
  m <- max(sum(y < a), sum(y > b))
  out["p.value"] <- min(1, (2 * (N - m) - 1) / N)
  if (q >= 1L) {
    a.s <- sort(a)
    b.s <- sort(b)
    lo <- a.s[q]
    hi <- b.s[N - q]
    out[c("pred.lower", "pred.upper")] <- c(lo, hi)
    out[c("CI.lower", "CI.upper")] <- c(y - hi, y - lo)
  }
  out
}

## =========================================================================
## 3. Unit scenarios for .conformal_interval
## =========================================================================

.U <- function(y, gamma, failed = NULL, M0 = rep(0, 9), R = 1:9) {
  fect:::.conformal_interval(M0, R, y, gamma, failed = failed)
}
.F1 <- c(TRUE, rep(FALSE, 8))
.nm <- c("pred.lower", "pred.upper", "CI.lower", "CI.upper", "p.value", "N",
         "q", "n.failed")

test_that("A0: names and order of the return vector", {
  v <- .U(0.5, 0.3)
  expect_identical(names(v), .nm)
})

test_that("A1: no failure", {
  v <- .U(0.5, 0.3)
  expect_equal(unname(v), c(-8, 8, -7.5, 8.5, 1, 10, 2, 0))
})

test_that("A2: one failure keeps N and q, widens", {
  v <- .U(0.5, 0.3, .F1)
  expect_equal(unname(v), c(-9, 9, -8.5, 9.5, 1, 10, 2, 1))
  v0 <- .U(0.5, 0.3)
  expect_true(v["pred.lower"] <= v0["pred.lower"] &&
                v["pred.upper"] >= v0["pred.upper"])
  expect_true(v["CI.lower"] <= v0["CI.lower"] && v["CI.upper"] >= v0["CI.upper"])
})

test_that("A3: values behind a failed slot are ignored; non-finite is failed", {
  ref <- .U(0.5, 0.3, .F1)
  M0 <- rep(0, 9); R <- 1:9
  M0[1] <- NA; R <- as.numeric(R); R[1] <- NA
  expect_equal(.U(0.5, 0.3, .F1, M0 = M0, R = R), ref)
  M0b <- rep(0, 9); M0b[1] <- NA
  expect_equal(.U(0.5, 0.3, NULL, M0 = M0b), ref)
  Rb <- as.numeric(1:9); Rb[1] <- Inf
  expect_equal(.U(0.5, 0.3, NULL, R = Rb), ref)
})

test_that("A4/A5: q failures give the whole real line", {
  f2 <- c(TRUE, TRUE, rep(FALSE, 7))
  v <- .U(0.5, 0.3, f2)
  expect_equal(unname(v), c(-Inf, Inf, -Inf, Inf, 1, 10, 2, 2))
  f89 <- (1:9) %in% c(8, 9)
  expect_equal(.U(0.5, 0.3, f89), v)
})

test_that("A6: q = 1 and one failure", {
  v <- .U(0.5, 0.1, .F1)
  expect_equal(unname(v[c("q", "N")]), c(1, 10))
  expect_equal(unname(v[c("CI.lower", "CI.upper")]), c(-Inf, Inf))
  v0 <- .U(0.5, 0.1)
  expect_equal(unname(v0[c("pred.lower", "pred.upper", "CI.lower", "CI.upper")]),
               c(-9, 9, -8.5, 9.5))
})

test_that("A7: p-value with a failure does not decrease", {
  v <- .U(5.5, 0.3)
  expect_equal(unname(v["p.value"]), 0.9)
  expect_equal(unname(v[c("CI.lower", "CI.upper")]), c(-2.5, 13.5))
  vf <- .U(5.5, 0.3, .F1)
  expect_equal(unname(vf["p.value"]), 1)
  expect_equal(unname(vf[c("CI.lower", "CI.upper")]), c(-3.5, 14.5))
  expect_gte(vf[["p.value"]], v[["p.value"]])
})

test_that("A8: boundary p <= gamma with 0 outside the CI", {
  v <- .U(9.5, 0.3)
  expect_equal(unname(v["p.value"]), 0.1)
  expect_equal(unname(v[c("CI.lower", "CI.upper")]), c(1.5, 17.5))
  expect_true(v[["p.value"]] <= 0.3 && v[["CI.lower"]] > 0)
  vf <- .U(9.5, 0.3, .F1)
  expect_equal(unname(vf["p.value"]), 0.3)
  expect_equal(unname(vf[c("CI.lower", "CI.upper")]), c(0.5, 18.5))
  expect_true(vf[["p.value"]] <= 0.3)
  expect_true(vf[["CI.lower"]] > 0)  # 0 outside
})

test_that("A9: every slot failed", {
  v <- fect:::.conformal_interval(rep(0, 5), rep(1, 5), 0.2, 0.5,
                                  failed = rep(TRUE, 5))
  expect_equal(unname(v), c(-Inf, Inf, -Inf, Inf, 1, 6, 2, 5))
})

test_that("A10: empty pool as before", {
  v <- fect:::.conformal_interval(numeric(0), numeric(0), 1, 0.1)
  expect_equal(unname(v[c("N", "q", "n.failed")]), c(1, 0, 0))
  expect_equal(unname(v[c("pred.lower", "pred.upper", "CI.lower", "CI.upper")]),
               c(-Inf, Inf, -Inf, Inf))
  expect_true(is.na(v[["p.value"]]))
})

test_that("A11: y = NA", {
  v <- fect:::.conformal_interval(c(0, 1, 2), c(1, 1, 1), NA_real_, 0.3)
  expect_true(all(is.na(v[c("CI.lower", "CI.upper", "p.value")])))
  expect_equal(v[["N"]], 4)
  l <- .cff_legacy(c(0, 1, 2), c(1, 1, 1), NA_real_, 0.3)
  expect_identical(v[1:7], l)
})

test_that("A12: failed of the wrong length is an error", {
  expect_error(fect:::.conformal_interval(rep(0, 9), 1:9, 0.5, 0.3,
                                          failed = c(TRUE, FALSE)))
  expect_error(fect:::.conformal_interval(rep(0, 9), 1:9, 0.5, 0.3,
                                          failed = rep(FALSE, 10)))
})

test_that("P-a: failed sets keep N and q, never narrow, match the oracle", {
  set.seed(3)
  J <- 39
  M0 <- rnorm(J)
  R <- abs(rnorm(J))
  for (gam in c(0.05, 0.1, 0.2)) {
    ref <- fect:::.conformal_interval(M0, R, y = 0.1, gamma = gam)
    for (k in 0:4) {
      for (rep in 1:5) {
        fl <- seq_len(J) %in% sample(J, k)
        v <- fect:::.conformal_interval(M0, R, y = 0.1, gamma = gam, failed = fl)
        N <- J + 1
        q <- fect:::.conformal_q(N, gam)
        expect_equal(unname(v[c("N", "q", "n.failed")]), c(N, q, k))
        expect_equal(unname(v["N"]), unname(ref["N"]))
        expect_equal(unname(v["q"]), unname(ref["q"]))
        expect_true(v[["pred.lower"]] <= ref[["pred.lower"]])
        expect_true(v[["pred.upper"]] >= ref[["pred.upper"]])
        a <- M0 - R; b <- M0 + R
        if (q >= 1) {
          expect_equal(v[["pred.lower"]], sort(replace(a, fl, -Inf))[q])
          expect_equal(v[["pred.upper"]], sort(replace(b, fl, Inf))[N - q])
        }
      }
    }
  }
})

test_that("P-b: duality and p-value formula over a grid", {
  set.seed(3)
  J <- 39
  M0 <- rnorm(J)
  R <- abs(rnorm(J))
  nviol <- 0
  for (k in 0:3) {
    fl <- seq_len(J) %in% sample(J, k)
    for (y in seq(-4, 4, by = 0.05)) {
      for (gam in c(0.05, 0.1)) {
        v <- fect:::.conformal_interval(M0, R, y = y, gamma = gam, failed = fl)
        out <- 0 < v[["CI.lower"]] || 0 > v[["CI.upper"]]
        if (!identical(v[["p.value"]] <= gam, out)) nviol <- nviol + 1
        a <- (M0 - R)[!fl]; b <- (M0 + R)[!fl]
        N <- J + 1
        m <- max(sum(y < a), sum(y > b))
        if (!isTRUE(all.equal(v[["p.value"]], min(1, (2 * (N - m) - 1) / N),
                              tolerance = 0))) nviol <- nviol + 1
      }
    }
  }
  expect_equal(nviol, 0)
})

test_that("P-c: legacy equivalence on all-finite inputs (identical)", {
  set.seed(11)
  for (i in 1:50) {
    J <- sample(1:60, 1)
    M0 <- rnorm(J); R <- abs(rnorm(J)); y <- rnorm(1) * 2
    gam <- sample(c(0.02, 0.05, 0.1, 0.2, 0.5), 1)
    new <- fect:::.conformal_interval(M0, R, y, gam)
    old <- .cff_legacy(M0, R, y, gam)
    expect_identical(new[1:7], old)
    expect_identical(unname(new[8]), 0)
  }
})

## =========================================================================
## 4. Integration scenarios
## =========================================================================

test_that("S1: gsynth, one treated unit, failures in the pool", {
  skip_on_cran()
  df <- .cff_sim(J = 40, Ntr = 1, TT = 20, T0 = 15, r = 1, seed = 21, tau = 6)
  fit_args <- list(method = "gsynth", r = 1, alpha = 0.1, keep.sims = TRUE)
  run <- function(df, ...) do.call(.cff_run, c(list(df), list(...), fit_args))
  ref <- do.call(.cff_fit, c(list(df), fit_args))
  cf0 <- ref$conformal
  ## S1-0
  expect_equal(unique(cf0$unit$N), 41)
  expect_equal(unique(cf0$unit$q), 2)
  expect_true(all(cf0$unit$n.failed == 0) && all(cf0$unit$n.excluded == 0))
  expect_equal(cf0$n.failed.fits, 0)
  expect_equal(cf0$n.failed.folds, 0)
  expect_equal(cf0$n.intervals.affected, 0)
  expect_identical(cf0$status, "ok")
  pt <- run(df)
  expect_identical(pt$fit$est.att, ref$est.att)
  expect_identical(pt$fit$est.avg, ref$est.avg)
  expect_identical(pt$fit$conformal$unit, ref$conformal$unit)
  expect_length(pt$ws, 0)
  expect_equal(pt$calls, 40)

  Mx <- cf0$folds$M0; Ex <- cf0$folds$E
  yv <- cf0$unit$Y
  oracle <- function(M, E, y, failcols, N = 41, q = 2) {
    a <- M - abs(E); b <- M + abs(E)
    a[failcols] <- -Inf; b[failcols] <- Inf
    m <- max(sum(y < a), sum(y > b))
    c(CI.lower = y - sort(b)[N - q], CI.upper = y - sort(a)[q],
      p.value = min(1, (2 * (N - m) - 1) / N))
  }

  ## S1-1
  r1 <- run(df, fail.calls = 7)
  f1 <- r1$fit; u1 <- f1$conformal$unit
  expect_equal(.cff_nwarn(r1), 1)
  expect_true(grepl("widened", r1$ws[grepl("failed", r1$ws)]))
  expect_match(r1$ws[grepl("failed", r1$ws)],
               "1 fold\\(s\\) failed \\(1 fit\\(s\\) returned no prediction\\)")
  expect_true(all(u1$n.failed == 1) && all(u1$n.excluded == 0))
  expect_true(all(u1$N == 41) && all(u1$q == 2))
  expect_equal(f1$conformal$n.failed.fits, 1)
  expect_equal(f1$conformal$n.failed.folds, 1)
  expect_equal(f1$conformal$n.intervals.affected, 5)
  expect_match(f1$conformal$status, "failed")
  expect_true(.cff_contains(u1, cf0$unit))
  expect_true(.cff_contains(.cff_post(f1), .cff_post(ref)))
  expect_true(.cff_contains(as.data.frame(f1$est.avg), as.data.frame(ref$est.avg)))
  for (s in 1:5) {
    o <- oracle(Mx[s, ], Ex[s, ], yv[s], 7)
    expect_equal(unname(c(u1$CI.lower[s], u1$CI.upper[s], u1$p.value[s])),
                 unname(o), tolerance = 1e-10)
    .cff_log(paste0("S1-1 unit row ", s, " CI.lower/CI.upper/p ref|new"),
             cf0$unit$CI.lower[s], u1$CI.lower[s], cf0$unit$CI.upper[s],
             u1$CI.upper[s], cf0$unit$p.value[s], u1$p.value[s])
  }
  Mbar <- colMeans(Mx); Ebar <- colMeans(Ex)
  oa <- oracle(Mbar, Ebar, mean(yv), 7)
  expect_equal(unname(c(f1$est.avg[1, "CI.lower"], f1$est.avg[1, "CI.upper"],
                        f1$est.avg[1, "p.value"])), unname(oa), tolerance = 1e-10)
  expect_true(all(is.na(f1$conformal$folds$M0[, 7])))
  expect_true(all(is.na(f1$conformal$folds$E[, 7])))
  expect_equal(f1$conformal$folds$M0[, -7], Mx[, -7])
  expect_equal(f1$conformal$folds$E[, -7], Ex[, -7])
  expect_equal(colnames(f1$conformal$folds$M0), colnames(Mx))
  expect_true(all(u1$p.value >= cf0$unit$p.value))
  .cff_same_points(f1, ref)
  .cff_invariants(f1, 0.1)
  .cff_invariants(ref, 0.1)
  expect_equal(f1$conformal$att$n.failed[!is.na(f1$conformal$att$N)], rep(1, 5))
  expect_equal(f1$conformal$avg$n.failed, c(1, 1))
  .cff_log("S1-1 est.avg CI ref|new", ref$est.avg[1, "CI.lower"], f1$est.avg[1, "CI.lower"],
           ref$est.avg[1, "CI.upper"], f1$est.avg[1, "CI.upper"])

  ## S1-2
  r2 <- run(df, fail.calls = c(7, 20))
  f2 <- r2$fit; u2 <- f2$conformal$unit
  expect_equal(.cff_nwarn(r2), 1)
  expect_true(all(u2$CI.lower == -Inf) && all(u2$CI.upper == Inf))
  expect_true(all(u2$n.failed == 2) && all(u2$N == 41) && all(u2$q == 2))
  p2 <- .cff_post(f2)
  expect_true(all(p2[, "CI.lower"] == -Inf) && all(p2[, "CI.upper"] == Inf))
  expect_equal(unname(f2$est.avg[1, c("CI.lower", "CI.upper")]), c(-Inf, Inf))
  expect_true(all(is.finite(u2$p.value)) && all(u2$p.value >= cf0$unit$p.value))
  expect_equal(f2$conformal$n.failed.folds, 2)
  .cff_same_points(f2, ref)
  .cff_invariants(f2, 0.1)

  ## S1-3
  df3 <- df
  ref3 <- .cff_fit(df3, method = "gsynth", r = 1, alpha = 0.05, keep.sims = TRUE)
  r3 <- .cff_run(df3, fail.calls = 7, method = "gsynth", r = 1, alpha = 0.05,
                 keep.sims = TRUE)
  u3 <- r3$fit$conformal$unit
  expect_equal(unique(ref3$conformal$unit$q), 1)
  expect_true(all(u3$CI.lower == -Inf) && all(u3$CI.upper == Inf))
  expect_true(all(u3$n.failed == 1) && all(u3$q == 1))
  expect_equal(.cff_nwarn(r3), 1)

  ## S1-4: fit not NULL but non-finite in the control's query column
  r4 <- run(df, nan.calls = 7, nan.fun = function(m) { m[, ncol(m)] <- NaN; m })
  u4 <- r4$fit$conformal$unit
  expect_equal(.cff_nwarn(r4), 1)
  expect_true(all(u4$n.failed == 1) && all(u4$N == 41))
  expect_equal(u4$CI.lower, u1$CI.lower)
  expect_equal(u4$CI.upper, u1$CI.upper)
  expect_match(r4$fit$conformal$status, "failed")
  expect_equal(r4$fit$conformal$n.failed.folds, 1)
  expect_match(r4$ws[grepl("failed", r4$ws)],
               "^vartype = \"conformal\": 1 fold\\(s\\) failed")
  expect_match(r4$fit$conformal$status, "1 failed fold\\(s\\)")
  .cff_log("S1-4 n.failed.fits/folds/affected", r4$fit$conformal$n.failed.fits,
           r4$fit$conformal$n.failed.folds, r4$fit$conformal$n.intervals.affected)

  ## S1-5: non-finite in unused pre-treatment rows
  r5 <- run(df, nan.calls = 7,
            nan.fun = function(m) { m[1:15, ncol(m)] <- NaN; m })
  expect_length(r5$ws, 0)
  expect_true(all(r5$fit$conformal$unit$n.failed == 0))
  expect_identical(r5$fit$est.att, ref$est.att)
  expect_identical(r5$fit$est.avg, ref$est.avg)
  expect_identical(r5$fit$conformal$unit, ref$conformal$unit)

  ## S1-6: duality in fits, tau = 6 and tau = 0
  dual <- function(f, alpha = 0.1) {
    u <- f$conformal$unit
    expect_identical(u$p.value <= alpha, (0 < u$CI.lower | 0 > u$CI.upper))
    ep <- .cff_post(f)
    expect_identical(unname(ep[, "p.value"] <= alpha),
                     unname(0 < ep[, "CI.lower"] | 0 > ep[, "CI.upper"]))
    e <- f$est.avg
    expect_identical(unname(e[1, "p.value"] <= alpha),
                     unname(0 < e[1, "CI.lower"] | 0 > e[1, "CI.upper"]))
  }
  for (f in list(ref, f1, f2)) dual(f)
  df0 <- .cff_sim(J = 40, Ntr = 1, TT = 20, T0 = 15, r = 1, seed = 21, tau = 0)
  ref0 <- do.call(.cff_fit, c(list(df0), fit_args))
  dual(ref0)
  dual(run(df0, fail.calls = 7)$fit)
  dual(run(df0, fail.calls = c(7, 20))$fit)

  ## S1-7 print and plot
  out <- capture.output(print(f1))
  expect_true(any(grepl("failed", out)))
  p <- plot(f1)
  expect_s3_class(p, "ggplot")
  p <- plot(f1, type = "gap")
  expect_s3_class(p, "ggplot")
  pp <- tryCatch(plot(f2), error = function(e) e)
  .cff_log("S1-7 plot(S1-2) class", class(pp)[1])
  pp <- tryCatch(plot(f2, type = "gap"), error = function(e) e)
  .cff_log("S1-7 plot(S1-2, gap) class", class(pp)[1])
})

test_that("S2: gsynth, unbalanced controls plus failures", {
  skip_on_cran()
  df <- .cff_sim(J = 25, seed = 8)
  drop <- (df$id %in% c(5, 9) & df$time == 4) | (df$id == 7 & df$time == 18)
  df <- df[!drop, ]
  args <- list(method = "gsynth", r = 1, alpha = 0.2, keep.sims = TRUE)
  ref <- do.call(.cff_fit, c(list(df), args))
  run <- function(...) do.call(.cff_run, c(list(df), list(...), args))
  cf <- ref$conformal; u <- cf$unit
  is18 <- u$time == 18
  ## S2-0
  expect_true(all(u$N[!is18] == 24) && all(u$n.excluded[!is18] == 2))
  expect_equal(u$N[is18], 23); expect_equal(u$n.excluded[is18], 3)
  expect_true(all(u$n.failed == 0))
  expect_equal(cf$n.controls, 25); expect_equal(cf$n.controls.ineligible, 2)
  ## S2-1
  r1 <- run(fail.calls = 11)
  u1 <- r1$fit$conformal$unit
  expect_equal(.cff_nwarn(r1), 1)
  expect_equal(u1$N, u$N); expect_equal(u1$q, u$q)
  expect_equal(u1$n.excluded, u$n.excluded)
  expect_true(all(u1$n.failed == 1))
  expect_true(.cff_contains(u1, u))
  expect_true(.cff_contains(.cff_post(r1$fit), .cff_post(ref)))
  .cff_invariants(r1$fit, 0.2)
  .cff_same_points(r1$fit, ref)
  ## S2-2: failed control is ineligible
  r2 <- run(fail.calls = 4)
  f2 <- r2$fit; u2 <- f2$conformal$unit
  expect_equal(.cff_nwarn(r2), 1)
  expect_equal(f2$conformal$n.failed.fits, 1)
  expect_true(all(u2$n.failed == 0))
  expect_equal(u2$N, u$N)
  for (cn in c("CI.lower", "CI.upper", "p.value", "N", "q")) {
    expect_identical(u2[[cn]], u[[cn]])
  }
  for (cn in c("CI.lower", "CI.upper", "p.value")) {
    expect_identical(f2$est.att[, cn], ref$est.att[, cn])
  }
  for (cn in c("N", "q")) {
    expect_identical(f2$conformal$att[[cn]], ref$conformal$att[[cn]])
  }
  expect_equal(f2$conformal$n.intervals.affected, 0)
  .cff_log("S2-2 status", f2$conformal$status)
  ## S2-3: failed control out of the pool at one date only
  r3 <- run(fail.calls = 6)
  f3 <- r3$fit; u3 <- f3$conformal$unit
  expect_equal(.cff_nwarn(r3), 1)
  expect_equal(u3$n.failed[is18], 0); expect_equal(u3$N[is18], 23)
  expect_identical(u3$CI.lower[is18], u$CI.lower[is18])
  expect_identical(u3$CI.upper[is18], u$CI.upper[is18])
  expect_true(all(u3$n.failed[!is18] == 1)); expect_true(all(u3$N[!is18] == 24))
  expect_true(.cff_contains(u3, u))
  expect_equal(f3$conformal$n.intervals.affected, 4)
  .cff_invariants(f3, 0.2)
})

test_that("S3: ife, two treated units adopting together", {
  skip_on_cran()
  df <- .cff_sim(J = 40, Ntr = 2, TT = 20, T0 = 15, r = 1, seed = 22, tau = 6)
  args <- list(method = "ife", r = 1, alpha = 0.1, keep.sims = TRUE)
  ref <- do.call(.cff_fit, c(list(df), args))
  run <- function(...) do.call(.cff_run, c(list(df), list(...), args))
  u0 <- ref$conformal$unit
  expect_true(all(u0$N == 41) && all(u0$q == 2) && all(u0$n.failed == 0))
  agg_ok <- function(f, nf) {
    cf <- f$conformal
    expect_true(all(cf$att$n.failed[!is.na(cf$att$N)] == nf))
    expect_true(all(cf$avg$n.failed == nf))
  }
  ## S3-1 / S3-2
  for (cl in c(14, 13)) {
    r <- run(fail.calls = cl)
    f <- r$fit; u <- f$conformal$unit
    bad.unit <- if (cl == 14) 2 else 1
    expect_equal(f$conformal$n.failed.fits, 1)
    expect_true(all(u$n.failed[u$unit == bad.unit] == 1))
    expect_true(all(u$n.failed[u$unit != bad.unit] == 0))
    agg_ok(f, 1)
    expect_true(all(u$N == 41))
    expect_equal(r$calls, 3 * 40)
    expect_true(.cff_contains(u, u0))
    expect_true(.cff_contains(.cff_post(f), .cff_post(ref)))
    expect_true(.cff_contains(as.data.frame(f$est.avg), as.data.frame(ref$est.avg)))
    expect_true(.cff_contains(as.data.frame(f$est.avg.unit),
                              as.data.frame(ref$est.avg.unit)))
    expect_equal(.cff_nwarn(r), 1)
    .cff_invariants(f, 0.1)
    .cff_same_points(f, ref)
    .cff_log(paste0("S3 fail.call ", cl, " unit n.failed by unit"),
             tapply(u$n.failed, u$unit, mean))
  }
  ## S3-3: control's query
  r <- run(fail.calls = 15)
  f <- r$fit; u <- f$conformal$unit
  expect_true(all(u$n.failed == 1)); agg_ok(f, 1)
  expect_equal(f$conformal$n.failed.fits, 1)
  expect_true(.cff_contains(u, u0))
  expect_equal(.cff_nwarn(r), 1)
  .cff_invariants(f, 0.1)
  ## S3-4
  r <- run(fail.calls = c(13, 14))
  f <- r$fit; u <- f$conformal$unit
  expect_equal(f$conformal$n.failed.fits, 2)
  expect_equal(f$conformal$n.failed.folds, 1)
  expect_true(all(u$n.failed == 1)); agg_ok(f, 1)
  expect_true(.cff_contains(u, u0))
  expect_equal(.cff_nwarn(r), 1)
  expect_match(r$ws[grepl("failed", r$ws)],
               "1 fold\\(s\\) failed \\(2 fit\\(s\\) returned no prediction\\)")
  .cff_invariants(f, 0.1)
})

test_that("S3-5: mc, control's query fails", {
  skip_on_cran()
  df <- .cff_sim(J = 25, Ntr = 1, seed = 5)
  args <- list(method = "mc", lambda = 0.01, alpha = 0.2, keep.sims = TRUE)
  ref <- do.call(.cff_fit, c(list(df), args))
  r <- do.call(.cff_run, c(list(df), list(fail.calls = 6), args))
  u0 <- ref$conformal$unit; u <- r$fit$conformal$unit
  expect_equal(unique(u0$N), 26); expect_equal(unique(u0$q), 3)
  expect_true(all(u$n.failed == 1))
  expect_true(all(r$fit$conformal$att$n.failed[!is.na(r$fit$conformal$att$N)] == 1))
  expect_true(all(r$fit$conformal$avg$n.failed == 1))
  expect_true(.cff_contains(u, u0))
  expect_equal(.cff_nwarn(r), 1)
  expect_equal(r$calls, 2 * 25)
  .cff_invariants(r$fit, 0.2)
})

## ---- S4 staggered ---------------------------------------------------------
.cff_s4_pool <- function(df, adopt, TT = 20) {
  tr.ids <- seq_along(adopt)
  ids <- sort(unique(df$id))
  co.ids <- setdiff(ids, tr.ids)
  Iobs <- matrix(FALSE, TT, max(ids))
  Iobs[cbind(df$time, df$id)] <- TRUE
  pool_cell <- function(i, t) {
    pre <- seq_len(adopt[i] - 1)
    co.ids[vapply(co.ids, function(j) all(Iobs[c(pre, t), j]), logical(1))]
  }
  list(co.ids = co.ids, pool_cell = pool_cell)
}

test_that("S4: gsynth staggered adoption, exclusion and failure mix", {
  skip_on_cran()
  adopt <- c(13, 13, 15, 15, 17, 17)
  df0 <- .cff_sim(J = 30, Ntr = 6, TT = 20, seed = 6, adopt = adopt)
  args <- list(method = "gsynth", r = 1, alpha = 0.2, keep.sims = TRUE)
  check_run <- function(df, fail.calls, failed.ids, label) {
    ref <- do.call(.cff_fit, c(list(df), args))
    r <- do.call(.cff_run, c(list(df), list(fail.calls = fail.calls), args))
    f <- r$fit; cf <- f$conformal; u <- cf$unit
    P <- .cff_s4_pool(df, adopt)
    Nco <- length(P$co.ids)
    expect_equal(cf$n.controls, Nco)
    ## unit rows
    for (s in seq_len(nrow(u))) {
      pool <- P$pool_cell(u$unit[s], u$time[s])
      expect_equal(u$n.excluded[s], Nco - length(pool))
      expect_equal(u$N[s], length(pool) + 1)
      expect_equal(u$n.failed[s], sum(failed.ids %in% pool))
    }
    ## aggregates: event time e
    att <- cf$att[!is.na(cf$att$N), ]
    cells <- u
    for (l in seq_len(nrow(att))) {
      e <- att$event.time[l]
      sel <- which(cells$event.time == e)
      pool <- Reduce(intersect, lapply(sel, function(s)
        P$pool_cell(cells$unit[s], cells$time[s])))
      expect_equal(att$n.excluded[l], Nco - length(pool), info = paste("e", e))
      expect_equal(att$N[l], length(pool) + 1)
      expect_equal(att$n.failed[l], sum(failed.ids %in% pool), info = paste("e", e))
    }
    poolall <- Reduce(intersect, lapply(seq_len(nrow(cells)), function(s)
      P$pool_cell(cells$unit[s], cells$time[s])))
    for (k in 1:2) {
      expect_equal(cf$avg$n.excluded[k], Nco - length(poolall))
      expect_equal(cf$avg$N[k], length(poolall) + 1)
      expect_equal(cf$avg$n.failed[k], sum(failed.ids %in% poolall))
    }
    ## containment
    expect_true(.cff_contains(u, ref$conformal$unit))
    expect_true(.cff_contains(.cff_post(f), .cff_post(ref)))
    expect_true(.cff_contains(as.data.frame(f$est.avg), as.data.frame(ref$est.avg)))
    expect_true(.cff_contains(as.data.frame(f$est.avg.unit),
                              as.data.frame(ref$est.avg.unit)))
    cal <- as.data.frame(f$est.eff.calendar); cal0 <- as.data.frame(ref$est.eff.calendar)
    expect_true(.cff_contains(cal, cal0))
    expect_true(sum(!is.na(cal0$CI.lower)) > 0)
    expect_equal(.cff_nwarn(r), 1)
    .cff_invariants(f, 0.2)
    .cff_same_points(f, ref)
    .cff_log(paste(label, "n.failed.fits/folds/affected/status"),
             cf$n.failed.fits, cf$n.failed.folds, cf$n.intervals.affected, cf$status)
    .cff_log(paste(label, "unit n.excluded range, n.failed range"),
             range(u$n.excluded), range(u$n.failed))
    .cff_log(paste(label, "att n.excluded by e"), att$n.excluded)
    .cff_log(paste(label, "att n.failed by e"), att$n.failed)
    .cff_log(paste(label, "avg n.excluded/n.failed"), cf$avg$n.excluded, cf$avg$n.failed)
    list(ref = ref, r = r)
  }
  ## S4-1
  o <- check_run(df0, fail.calls = 4, failed.ids = 10, "S4-1")
  cf <- o$r$fit$conformal
  expect_true(all(cf$unit$n.failed == 1) && all(cf$unit$n.excluded == 0))
  expect_true(all(cf$unit$N == 31))
  expect_equal(cf$n.failed.fits, 1)
  ## S4-2
  df2 <- df0[!(df0$id == 12 & df0$time == 15), ]
  o <- check_run(df2, fail.calls = 4, failed.ids = 10, "S4-2")
  cf <- o$r$fit$conformal
  expect_true(any(cf$unit$n.excluded > 0) && any(cf$unit$n.excluded == 0))
  att <- cf$att[!is.na(cf$att$N), ]
  expect_true(all(att$n.excluded[att$event.time %in% 1:4] == 1))
  expect_true(all(att$n.excluded[att$event.time >= 5] == 0))
  expect_true(all(cf$avg$n.excluded == 1))
  ## S4-3
  o <- check_run(df2, fail.calls = 6, failed.ids = 12, "S4-3")
  cf <- o$r$fit$conformal
  u <- cf$unit
  att <- cf$att[!is.na(cf$att$N), ]
  expect_true(all(att$n.failed[att$event.time %in% 1:4] == 0))
  expect_true(all(att$n.failed[att$event.time >= 5] == 1))
  expect_true(all(cf$avg$n.failed == 0))
  ok.rows <- u$unit <= 4 & u$time != 15
  expect_true(all(u$n.failed[ok.rows] == 1 & u$n.excluded[ok.rows] == 0))
  expect_equal(cf$n.failed.fits, 1)
})

## ---- S5 regression --------------------------------------------------------
.cff_configs <- function() {
  c5 <- .cff_sim(J = 25, seed = 8)
  c5 <- c5[!((c5$id %in% c(5, 9) & c5$time == 4) | (c5$id == 7 & c5$time == 18)), ]
  c7 <- .cff_sim(J = 20, seed = 7); c7$Y <- c7$Y * 50
  list(
    C1 = list(df = .cff_sim(J = 20, Ntr = 1, seed = 1),
              args = list(method = "gsynth", r = 1, tol = 1e-5,
                          max.iteration = 5000)),
    C2 = list(df = .cff_sim(J = 30, Ntr = 6, TT = 20, seed = 6,
                            adopt = c(13, 13, 15, 15, 17, 17)),
              args = list(method = "gsynth", r = 1, tol = 1e-5,
                          max.iteration = 5000)),
    C3 = list(df = .cff_sim(J = 25, Ntr = 2, seed = 5),
              args = list(method = "ife", r = 1, tol = 1e-5,
                          max.iteration = 5000)),
    C4 = list(df = .cff_sim(J = 25, Ntr = 1, seed = 5),
              args = list(method = "mc", lambda = 0.01, tol = 1e-5,
                          max.iteration = 5000)),
    C5 = list(df = c5, args = list(method = "gsynth", r = 1, tol = 1e-5,
                          max.iteration = 5000)),
    C6 = list(df = .cff_sim(J = 8, seed = 9),
              args = list(method = "gsynth", r = 1, tol = 1e-5,
                          max.iteration = 5000)),
    C7 = list(df = c7, args = list(method = "gsynth", r = 1, tol = 1e-10,
                                   max.iteration = 5000, normalize = TRUE))
  )
}

.cff_extract <- function(fit) {
  cf <- fit$conformal
  list(est.att = fit$est.att, est.att90 = fit$est.att90,
       att.bound = fit$att.bound, est.avg = fit$est.avg,
       est.avg.unit = fit$est.avg.unit, est.eff.calendar = fit$est.eff.calendar,
       unit = cf$unit[, c("unit", "time", "event.time", "Y", "Y.ct", "eff",
                          "CI.lower", "CI.upper", "p.value", "N", "q")],
       att = cf$att[, c("event.time", "N", "q")],
       avg = cf$avg[, c("type", "N", "q")], status = cf$status,
       n.failed.fits = cf$n.failed.fits,
       n.controls.ineligible = cf$n.controls.ineligible,
       N.range = cf$N.range, r = cf$r, lambda = cf$lambda)
}

test_that("S5-2: no-failure output equals the e99f2d2 baseline", {
  skip_on_cran()
  fx <- Sys.getenv("CFF_BASELINE",
                   testthat::test_path("fixtures", "conformal-baseline-e99f2d2.rds"))
  skip_if_not(file.exists(fx), "baseline fixture not found")
  base <- readRDS(fx)
  cfg <- .cff_configs()
  for (nm in names(cfg)) {
    fit <- NULL
    expect_no_warning(fit <- do.call(.cff_fit, c(list(cfg[[nm]]$df), cfg[[nm]]$args)))
    now <- .cff_extract(fit)
    expect_equal(now, base[[nm]], tolerance = 1e-10, info = nm)
    expect_identical(now$status, base[[nm]]$status)
    expect_true(all(fit$conformal$unit$n.failed == 0), info = nm)
    expect_true(all(fit$conformal$att$n.failed[!is.na(fit$conformal$att$N)] == 0))
    expect_true(all(fit$conformal$avg$n.failed == 0))
    if (nm == "C5") {
      u <- fit$conformal$unit
      expect_equal(u$n.excluded, fit$conformal$n.controls - (u$N - 1))
    }
  }
})

## ---- edge cases -----------------------------------------------------------
test_that("E1-E5: edge cases", {
  skip_on_cran()
  df <- .cff_sim(J = 20, Ntr = 1, seed = 1)
  ## E1: every fold fails
  ref <- .cff_fit(df, method = "gsynth", r = 1, keep.sims = TRUE)
  r <- .cff_run(df, all.fail = TRUE, method = "gsynth", r = 1, keep.sims = TRUE)
  f <- r$fit; u <- f$conformal$unit
  expect_true(all(u$CI.lower == -Inf) && all(u$CI.upper == Inf))
  expect_equal(u$n.failed, f$conformal$n.controls - u$n.excluded)
  expect_equal(u$N, ref$conformal$unit$N)
  expect_equal(.cff_nwarn(r), 1)
  expect_true(all(is.finite(u$p.value)))
  pre <- as.numeric(rownames(f$est.att)) < 1
  expect_true(all(is.na(f$est.att[pre, "CI.lower"])))
  .cff_invariants(f, 0.05)
  ## E2: q == 0 with a failure
  d8 <- .cff_sim(J = 8, seed = 9)
  ref8 <- .cff_fit(d8, method = "gsynth", r = 1)
  r8 <- .cff_run(d8, fail.calls = 3, method = "gsynth", r = 1)
  u8 <- r8$fit$conformal$unit
  expect_true(all(u8$CI.lower == -Inf) && all(u8$CI.upper == Inf))
  expect_true(all(u8$n.failed == 1) && all(u8$q == 0))
  expect_match(r8$fit$conformal$status, "infinite")
  expect_match(r8$fit$conformal$status, "failed")
  expect_equal(.cff_nwarn(r8), 1)
  .cff_log("E2 status", r8$fit$conformal$status, "| ref:", ref8$conformal$status)
  ## E3: one warning per fect() call
  r3 <- .cff_run(df, fail.calls = c(2, 7, 11), method = "gsynth", r = 1)
  expect_equal(.cff_nwarn(r3), 1)
  expect_equal(r3$fit$conformal$n.failed.fits, 3)
  expect_equal(r3$fit$conformal$n.failed.folds, 3)
  expect_match(r3$ws[grepl("failed", r3$ws)],
               "3 fold\\(s\\) failed \\(3 fit\\(s\\) returned no prediction\\)")
  ## E4: keep.sims = FALSE
  r4 <- .cff_run(df, fail.calls = 7, method = "gsynth", r = 1, keep.sims = FALSE)
  expect_null(r4$fit$conformal$folds)
  expect_true(all(r4$fit$conformal$unit$n.failed == 1))
  expect_equal(r4$fit$conformal$n.failed.fits, 1)
  expect_equal(.cff_nwarn(r4), 1)
  ## E5: ife, a control whose fold hosts no window fails (treated query call)
  d <- .cff_sim(J = 10, Ntr = 1, seed = 3)
  d <- d[!(d$id == 5 & d$time == 4), ]
  ref5 <- .cff_fit(d, method = "ife", r = 1, alpha = 0.2)
  r5 <- .cff_run(d, fail.calls = 7, method = "ife", r = 1, alpha = 0.2)
  expect_equal(r5$calls, 2 * 9 + 1)
  expect_equal(r5$fit$conformal$n.failed.fits, 1)
  expect_equal(.cff_nwarn(r5), 1)
  expect_equal(r5$fit$conformal$n.intervals.affected, 0)
  expect_true(all(r5$fit$conformal$unit$n.failed == 0))
  expect_identical(r5$fit$conformal$unit$CI.lower, ref5$conformal$unit$CI.lower)
  expect_identical(r5$fit$conformal$unit$CI.upper, ref5$conformal$unit$CI.upper)
})
