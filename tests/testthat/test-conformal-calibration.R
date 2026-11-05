## vartype = "conformal": calibration q, duality of the p-value, aggregated
## targets and the finite reference completion, tested on the internal
## functions .conformal_q() and .conformal_interval() (pure, fast, no fits).
## Spec: statsclaw-workspace/fect/runs/2026-10-28-fect-248-assembly/conformal/
## test-spec.md, sections 4.1 to 4.4.

## ---- 4.1 boundary table ----------------------------------------------------

test_that(".conformal_q gives the exact q at and next to each cutoff", {
  tab <- list(
    list(10, 0.1, 1L), list(10, 0.1 - 1e-12, 0L), list(10, 0.1 + 1e-12, 1L),
    list(10, 0.3, 2L), list(10, 0.3 - 1e-12, 1L),
    list(20, 0.05, 1L), list(20, 0.05 - 1e-12, 0L),
    list(100, 0.29, 15L), list(100, 0.29 - 1e-12, 14L),
    list(5, 0.2, 1L), list(5, 0.2 - 1e-12, 0L),
    list(4, 0.75, 2L), list(4, 0.75 - 1e-12, 1L),
    list(2, 0.5, 1L), list(2, 0.5 - 1e-12, 0L),
    list(7, 3 / 7, 2L), list(7, 3 / 7 - 1e-12, 1L),
    list(1, 0.5, 0L), list(10, 0, 0L), list(10, 0.9, 5L), list(10, 1, 5L),
    list(31, 0.1, 2L), list(61, 0.1, 3L), list(40, 0.05, 1L),
    list(0, 0.1, 0L), list(NA, 0.1, 0L), list(10, NA, 0L)
  )
  for (row in tab) {
    expect_identical(fect:::.conformal_q(row[[1]], row[[2]]), row[[3]],
                     info = paste("N =", row[[1]], "gamma =",
                                  format(row[[2]], digits = 17)))
  }
  expect_identical(fect:::.conformal_q(10L, 0.1), fect:::.conformal_q(10, 0.1))
})

## ---- 4.2 exhaustive agreement with an integer oracle ------------------------

.cc_oracle <- function(N, a, s) {
  qs <- 0:floor(N / 2)
  as.integer(max(qs[qs == 0 | (2 * qs - 1) * 10^s <= a * N]))
}

test_that(".conformal_q equals the integer oracle for two- and three-decimal gamma", {
  bad <- 0L
  for (N in 1:400) for (a in 0:100) {
    if (!identical(fect:::.conformal_q(N, a / 100), .cc_oracle(N, a, 2))) bad <- bad + 1L
  }
  expect_identical(bad, 0L)
  bad <- 0L
  for (N in 1:400) for (a in 1:99) {
    g <- as.numeric(sprintf("0.%02d", a))          # a typed literal
    if (!identical(fect:::.conformal_q(N, g), .cc_oracle(N, a, 2))) bad <- bad + 1L
  }
  expect_identical(bad, 0L)
  bad <- 0L
  for (N in 1:150) for (a in 0:1000) {
    if (!identical(fect:::.conformal_q(N, a / 1000), .cc_oracle(N, a, 3))) bad <- bad + 1L
  }
  expect_identical(bad, 0L)
})

test_that(".conformal_q is monotone, capped, and positive exactly when 1/N <= gamma", {
  g <- (0:99) / 100
  Q <- outer(1:100, g, Vectorize(function(N, gm) fect:::.conformal_q(N, gm)))
  expect_true(all(apply(Q, 1, function(x) all(diff(x) >= 0))))   # in gamma
  expect_true(all(apply(Q, 2, function(x) all(diff(x) >= 0))))   # in N
  expect_true(all(Q <= floor((1:100) / 2)))
  pos <- outer(1:100, g, function(N, gm) 1 / N <= gm)
  expect_identical(Q > 0, pos)
})

## ---- 4.3 duality at the boundary and on a grid ------------------------------

test_that("p-value and interval agree at the cutoff alpha = 0.1, N = 10", {
  ci <- function(gm) fect:::.conformal_interval(rep(0, 9), 1:9, 20, gm)
  lo <- ci(0.1 - 1e-12)
  expect_equal(unname(lo[c("q", "pred.lower", "pred.upper", "CI.lower", "CI.upper")]),
               c(0, -Inf, Inf, -Inf, Inf))
  expect_equal(unname(lo["p.value"]), 0.1)
  expect_false(lo["p.value"] <= 0.1 - 1e-12)
  at <- ci(0.1)
  expect_equal(unname(at[c("q", "pred.lower", "pred.upper", "CI.lower", "CI.upper")]),
               c(1, -9, 9, 11, 29))
  expect_equal(unname(at["p.value"]), 0.1)
  expect_true(at["p.value"] <= 0.1)
  hi <- ci(0.1 + 1e-12)
  expect_equal(unname(hi[c("q", "CI.lower", "CI.upper", "p.value")]),
               unname(at[c("q", "CI.lower", "CI.upper", "p.value")]))
})

test_that("p.value <= gamma exactly when 0 is outside the interval, over a grid", {
  viol <- 0L; formula.bad <- 0L; n <- 0L
  for (N in c(2, 4, 5, 7, 10, 20, 41, 100)) {
    J <- N - 1
    set.seed(N)
    M0 <- rnorm(J); R <- abs(rnorm(J))
    for (g in c(0.05, 0.1, 0.15, 0.2, 0.29, 0.3, 0.5, 0.75)) {
      for (d in c(-1e-12, 0, 1e-12)) {
        gm <- g + d
        for (k in 0:2) {
          failed <- if (k == 0) NULL else {
            f <- rep(FALSE, J); f[sample.int(J, min(k, J))] <- TRUE; f
          }
          fl <- if (is.null(failed)) rep(FALSE, J) else failed
          a <- (M0 - R)[!fl]; b <- (M0 + R)[!fl]
          for (y in seq(-6, 6, by = 0.05)) {
            o <- fect:::.conformal_interval(M0, R, y, gm, failed)
            n <- n + 1L
            if (!identical(unname(o["p.value"] <= gm),
                           unname(0 < o["CI.lower"] || 0 > o["CI.upper"]))) {
              viol <- viol + 1L
            }
            m <- max(sum(y < a), sum(y > b))
            if (!identical(unname(o["p.value"]), min(1, (2 * (N - m) - 1) / N))) {
              formula.bad <- formula.bad + 1L
            }
          }
        }
      }
    }
  }
  expect_gt(n, 10000L)
  expect_identical(viol, 0L)
  expect_identical(formula.bad, 0L)
})

test_that("2 * alpha is exact for two-decimal alpha", {
  for (k in 1:49) expect_identical(2 * (k / 100), (2 * k) / 100)
})

## ---- 4.4 aggregated targets and the finite reference completion -------------

test_that("averaging treated units can cover far below nominal (0.42 at 0.90)", {
  ## nine controls, 100 treated units, iid Bernoulli(0.01) outcomes, zero
  ## predictor, no effect; only whether some control has error 1 matters
  covered <- function(R, cnt) {
    o <- fect:::.conformal_interval(rep(0, 9), R, cnt / 100, 0.1)
    o[["CI.lower"]] <= 0 && 0 <= o[["CI.upper"]]
  }
  one <- c(1, rep(0, 8)); none <- rep(0, 9)
  expect_true(all(vapply(0:100, function(cnt) covered(one, cnt), logical(1))))
  cov.none <- vapply(0:100, function(cnt) covered(none, cnt), logical(1))
  expect_true(cov.none[1])
  expect_false(any(cov.none[-1]))
  coverage <- (1 - 0.99^9) + 0.99^9 * dbinom(0, 100, 0.01)
  expect_equal(coverage, 0.4208596, tolerance = 1e-7)   # 0.42085961, given to 7 digits
  expect_equal(coverage, 1 - 0.99^9 * (1 - 0.99^100), tolerance = 1e-12)
  expect_lt(coverage, 0.5)
})

test_that("failed slots never give a narrower interval than a finite completion", {
  set.seed(11)
  for (i in 1:300) {
    J <- sample(c(9, 19, 39), 1)
    gm <- sample(c(0.05, 0.1, 0.2, 0.3), 1)
    M0 <- rnorm(J); R <- abs(rnorm(J)); y <- runif(1, -3, 3)
    nf <- sample(0:min(J, 4), 1)
    failed <- rep(FALSE, J); failed[sample.int(J, nf)] <- TRUE
    M0s <- M0; Rs <- R
    for (j in which(failed)) {
      M0s[j] <- rnorm(1) * 10^sample(0:3, 1)
      Rs[j] <- abs(rnorm(1)) * 10^sample(0:3, 1)
    }
    fe <- fect:::.conformal_interval(M0, R, y, gm, failed)
    ref <- fect:::.conformal_interval(M0s, Rs, y, gm)
    expect_identical(unname(fe[c("N", "q")]), unname(ref[c("N", "q")]))
    expect_true(fe[["pred.lower"]] <= ref[["pred.lower"]])
    expect_true(fe[["pred.upper"]] >= ref[["pred.upper"]])
    expect_true(fe[["CI.lower"]] <= ref[["CI.lower"]])
    expect_true(fe[["CI.upper"]] >= ref[["CI.upper"]])
    expect_true(fe[["p.value"]] >= ref[["p.value"]])
  }
})
