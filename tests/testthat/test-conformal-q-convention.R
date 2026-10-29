## vartype = "conformal": exact q convention, alpha guard, failure bookkeeping
## (fect 2.4.9). See R/conformal.R (.conformal_q, conformal_jplus) and the
## alpha guard in fect.default(). Small panels; parallel = FALSE.

.cq_sim <- function(J = 12, TT = 14, T0 = 10, seed = 5) {
  set.seed(seed)
  N <- J + 1
  Y0 <- outer(rnorm(TT), rep(1, N)) + outer(rep(1, TT), rnorm(N)) +
    matrix(rnorm(TT * N), TT, N)
  D <- matrix(0, TT, N)
  D[seq_len(TT) > T0, 1] <- 1
  data.frame(id = rep(seq_len(N), each = TT), time = rep(seq_len(TT), N),
             Y = c(Y0 + D), D = c(D))
}

.cq_fit <- function(df, ...) {
  suppressMessages(fect(Y ~ D, data = df, index = c("id", "time"),
                        se = TRUE, vartype = "conformal", parallel = FALSE,
                        CV = FALSE, ...))
}

## integer-arithmetic oracle for a gamma = a / 10^s: (2q - 1) * 10^s <= a * N
.cq_oracle <- function(N, a, s) {
  qs <- 0:floor(N / 2)
  ok <- qs == 0 | (2 * qs - 1) * 10^s <= a * N
  max(qs[ok])
}

test_that("q uses the exact comparison fl((2q - 1) / N) <= gamma, no tolerance", {
  q <- fect:::.conformal_q
  ## gamma equal to a cutoff as a rational number admits it
  expect_equal(q(10, 0.1), 1L)
  expect_equal(q(20, 0.05), 1L)
  expect_equal(q(100, 0.29), 15L)   # 0.29 * 100 is 28.999999999999996
  expect_equal(q(2, 0.5), 1L)
  expect_equal(q(4, 0.75), 2L)
  expect_equal(q(1, 0.5), 0L)
  ## gamma below a cutoff, however slightly, does not
  expect_equal(q(10, 0.1 - 1e-12), 0L)
  expect_equal(q(10, 0.3 - 1e-12), 1L)
  expect_equal(q(20, 0.05 - 1e-12), 0L)
  expect_equal(q(100, 0.29 - 1e-12), 14L)
  expect_equal(q(2, 0.5 - 1e-12), 0L)
  expect_equal(q(4, 0.75 - 1e-12), 1L)
  ## a computed alpha is treated as the double it is
  expect_equal(q(10, 1 - 0.9), 0L)
  ## degenerate inputs give 0
  expect_equal(q(0, 0.1), 0L)
  expect_equal(q(NA_real_, 0.1), 0L)
  expect_equal(q(Inf, 0.1), 0L)
  expect_equal(q(10, NA_real_), 0L)
  expect_equal(q(10, NaN), 0L)
  ## gamma >= 1 is capped at floor(N / 2)
  expect_equal(q(11, 1), 5L)
  expect_equal(q(11, 2), 5L)
  ## integer class
  expect_type(q(10, 0.1), "integer")
})

test_that("q equals the integer-arithmetic oracle for decimal gamma", {
  q <- fect:::.conformal_q
  for (a in 0:100) {
    for (N in 1:200) {
      expect_identical(q(N, a / 100), as.integer(.cq_oracle(N, a, 2)))
    }
  }
  for (a in seq(0, 1000, by = 7)) {
    for (N in c(1, 2, 3, 10, 41, 100, 299)) {
      expect_identical(q(N, a / 1000), as.integer(.cq_oracle(N, a, 3)))
    }
  }
})

test_that("the p-value stays dual to the interval at the cutoffs", {
  set.seed(31)
  for (N in c(2, 4, 5, 10, 20, 41, 100)) {
    M0 <- rnorm(N - 1)
    R <- abs(rnorm(N - 1))
    for (g0 in c(0.05, 0.1, 0.15, 0.2, 0.29, 0.3, 0.5, 0.75)) {
      for (gam in c(g0, g0 - 1e-12, g0 + 1e-12)) {
        for (y in seq(-6, 6, by = 0.1)) {
          v <- fect:::.conformal_interval(M0, R, y = y, gamma = gam)
          out <- 0 < v["CI.lower"] || 0 > v["CI.upper"]
          expect_identical(unname(v["p.value"] <= gam), unname(out))
        }
      }
    }
  }
})

test_that("a fit uses the typed alpha: 0.1 and 1 - 0.9 differ at N = 10", {
  df <- .cq_sim(J = 9)  # N = 10
  f1 <- .cq_fit(df, method = "ife", r = 1, alpha = 0.1)
  expect_true(all(f1$conformal$unit$q == 1))
  expect_true(all(is.finite(f1$conformal$unit$CI.lower)))
  f2 <- .cq_fit(df, method = "ife", r = 1, alpha = 1 - 0.9)
  expect_true(all(f2$conformal$unit$q == 0))
  expect_true(all(is.infinite(f2$conformal$unit$CI.lower)))
  ## 2 * alpha is exact: the 90% band at alpha = 0.05 is the q of 0.1
  f3 <- .cq_fit(df, method = "ife", r = 1, alpha = 0.05)
  expect_true(all(f3$conformal$unit$q == 0))
  expect_true(all(is.finite(f3$att.bound[, "CI.lower"][!is.na(f3$att.bound[, "CI.lower"])])))
})

test_that("alpha must lie in (0, 1) for a conformal fit", {
  df <- .cq_sim()
  msg <- "needs 0 < alpha < 1"
  expect_error(.cq_fit(df, method = "ife", r = 1, alpha = 0), msg)
  expect_error(.cq_fit(df, method = "ife", r = 1, alpha = 1), msg)
  expect_error(.cq_fit(df, method = "ife", r = 1, alpha = -0.1), msg)
  expect_error(.cq_fit(df, method = "ife", r = 1, alpha = 1.5), msg)
  expect_error(.cq_fit(df, method = "ife", r = 1, alpha = NA_real_), msg)
  expect_error(.cq_fit(df, method = "ife", r = 1, alpha = c(0.05, 0.1)),
               "got alpha = 0.05, 0.10")
  expect_error(.cq_fit(df, method = "ife", r = 1, alpha = "a"), msg)
  ## alpha >= 0.5 is allowed (the 90% band is built at 2 * alpha)
  fit <- .cq_fit(df, method = "ife", r = 1, alpha = 0.6)
  expect_s3_class(fit, "fect")
  expect_equal(fit$conformal$gamma, 0.6)
})

test_that("failure bookkeeping counts controls, not fits, and matches the warning", {
  df <- .cq_sim()
  real <- fect:::.conformal_component_fit
  mk <- function(mode) {
    n <- 0L
    function(...) {
      n <<- n + 1L
      out <- real(...)
      if (n == 3L) {
        if (mode == "null") return(NULL)
        if (mode == "nan") {
          out[, ncol(out)] <- NaN
        }
      }
      out
    }
  }
  ## a fit that returns nothing
  testthat::local_mocked_bindings(.conformal_component_fit = mk("null"),
                                  .package = "fect")
  expect_warning(fit <- .cq_fit(df, method = "ife", r = 1),
                 "failed.*widened", perl = TRUE)
  cf <- fit$conformal
  expect_equal(cf$n.failed.fits, 1)
  expect_equal(cf$n.failed.folds, 1)
  expect_match(cf$status, "; 1 failed fold\\(s\\), affected intervals widened$")
  expect_true(cf$n.intervals.affected >= 1)
  expect_true(all(cf$unit$n.failed %in% 0:1))
  expect_true(all(cf$unit$N == cf$N.range[1]))
})

test_that("a non-finite fold quantity is counted as a failed fold", {
  df <- .cq_sim()
  real <- fect:::.conformal_component_fit
  n <- 0L
  testthat::local_mocked_bindings(
    .conformal_component_fit = function(...) {
      n <<- n + 1L
      out <- real(...)
      if (n == 3L) out[, ncol(out)] <- NaN
      out
    },
    .package = "fect"
  )
  expect_warning(fit <- .cq_fit(df, method = "ife", r = 1),
                 "1 fold\\(s\\) failed \\(0 fit\\(s\\) returned no prediction\\)")
  cf <- fit$conformal
  expect_equal(cf$n.failed.fits, 0)
  expect_equal(cf$n.failed.folds, 1)
  expect_match(cf$status, "; 1 failed fold\\(s\\), affected intervals widened$")
  expect_false(grepl("0 failed", cf$status))
})

test_that("a fit without failures has the plain status and no warning", {
  df <- .cq_sim()
  expect_no_warning(fit <- .cq_fit(df, method = "ife", r = 1, alpha = 0.1))
  expect_identical(fit$conformal$status, "ok")
  expect_equal(fit$conformal$n.failed.folds, 0)
  expect_equal(fit$conformal$n.failed.fits, 0)
})

test_that("the aggregate metadata makes no conservative claim", {
  df <- .cq_sim()
  fit <- .cq_fit(df, method = "ife", r = 1)
  ext <- fit$conformal$extension
  expect_false(grepl("conservativ", ext, ignore.case = TRUE))
  expect_match(ext, "not established by the unit-level theorem", fixed = TRUE)
  expect_match(ext, "separate assumptions and validation", fixed = TRUE)
})

test_that("a conformal fit prints and summarises on the merged code", {
  df <- .cq_sim()
  fit <- .cq_fit(df, method = "ife", r = 1)
  expect_null(fit$boot.diag)
  out <- capture.output(print(fit))
  expect_true(any(grepl("conformal, pooled information-matched Jackknife+", out,
                        fixed = TRUE)))
  out.s <- capture.output(summary(fit))
  expect_true(any(grepl("conformal, pooled information-matched Jackknife+",
                        out.s, fixed = TRUE)))
  expect_false(any(grepl("Bootstrap check", out.s, ignore.case = TRUE)))
})
