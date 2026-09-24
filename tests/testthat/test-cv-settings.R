## The cross-validation uses the user's CV settings in every case (fixed in
## 2.4.6). Before the fix:
##   * with se = TRUE, fect_boot() did not pass cv.rule, cv.buffer, cv.donut,
##     min.T0 or proportion on to fect_cv(), so the CV ran on fect_cv()'s
##     defaults ("1se", 1, 1, 5, 0) whatever the user set;
##   * never-treated CV (method = "gsynth", or "ife" with
##     time.component.from = "nevertreated") never received cv.rule and
##     always used "1se", with or without se.
## Found 2026-09-23. With the same seed, a run with se = TRUE (parallel
## bootstrap, which seeds with set.seed(seed) before its CV) now draws the
## same folds and reports the same CV table and r as a run with se = FALSE
## after set.seed(seed).

skip_on_cran()

## N = 30, T = 16, 10 treated units adopting at t = 11..13, two factors with
## the second one weak (loadings scaled by w = 0.3), one covariate. On these
## panels cv.rule "min" and "1se" pick different r: panel 1 for ife
## (0 vs 1), panel 2 for gsynth (1 vs 2), with set.seed(1) before the call.
.rule_panel <- function(panel) {
  withr::with_seed(panel, {
    N <- 30; TT <- 16; Ntr <- 10
    F <- matrix(rnorm(TT * 2), TT, 2)
    L <- matrix(rnorm(N * 2), N, 2); L[, 2] <- 0.3 * L[, 2]
    a <- rnorm(N); xi <- rnorm(TT); X <- matrix(rnorm(N * TT), TT, N)
    T0 <- rep(Inf, N); T0[seq_len(Ntr)] <- sample(11:13, Ntr, replace = TRUE)
    D <- sapply(seq_len(N), function(i) as.integer(seq_len(TT) >= T0[i]))
    Y <- outer(xi, rep(1, N)) + outer(rep(1, TT), a) + F %*% t(L) + 0.5 * X +
      2 * D + matrix(rnorm(N * TT), TT, N)
    data.frame(id = rep(seq_len(N), each = TT), time = rep(seq_len(TT), N),
               Y = c(Y), D = c(D), X = c(X))
  })
}

## Set the session RNG to `session`, then fit with CV over r = 0..4.
.fit_cv <- function(d, session, ...) {
  set.seed(session)
  suppressMessages(suppressWarnings(
    fect(Y ~ D + X, data = d, index = c("id", "time"), CV = TRUE, r = c(0, 4), ...)
  ))
}

## se = TRUE with the parallel bootstrap and seed = s seeds with set.seed(s)
## before its CV, the same state as .fit_cv(session = s) without a seed.
.fit_boot <- function(d, seed, ...) {
  .fit_cv(d, session = 99, seed = seed, se = TRUE, nboots = 10,
          parallel = TRUE, cores = 2, ...)
}

test_that("se = TRUE cross-validates with the user's cv.buffer, cv.donut, min.T0 and proportion", {
  d <- .rule_panel(1)
  ## cv.donut only affects block folds, so it gets its own case. (Block folds
  ## score only cells more than cv.donut steps from either end of a held-out
  ## run of cv.nobs = 3 cells, so cv.donut = 2 would leave nothing to score.)
  cases <- list(
    rolling = list(cv.buffer = 2, min.T0 = 6, proportion = 0.5),
    block   = list(cv.method = "block", cv.donut = 0)
  )
  for (nm in names(cases)) {
    set_args <- cases[[nm]]
    base_args <- if (nm == "block") list(cv.method = "block") else list()
    default <- do.call(.fit_cv, c(list(d = d, session = 3, method = "ife"), base_args))
    nose    <- do.call(.fit_cv, c(list(d = d, session = 3, method = "ife"), set_args))
    boot    <- do.call(.fit_boot, c(list(d = d, seed = 3, method = "ife"), set_args))
    ## The settings change the CV table on this panel ...
    expect_false(isTRUE(all.equal(default$CV.out.ife, nose$CV.out.ife)), info = nm)
    ## ... and se = TRUE now uses them.
    expect_equal(boot$CV.out.ife, nose$CV.out.ife, tolerance = 1e-10, info = nm)
    expect_equal(boot$r.cv, nose$r.cv, info = nm)
  }
})

test_that("se = TRUE applies cv.rule", {
  d <- .rule_panel(1)
  nose_1se <- .fit_cv(d, session = 1, method = "ife")
  nose_min <- .fit_cv(d, session = 1, method = "ife", cv.rule = "min")
  skip_if(nose_1se$r.cv == nose_min$r.cv, "this panel no longer separates the rules")
  boot_min <- .fit_boot(d, seed = 1, method = "ife", cv.rule = "min")
  expect_equal(boot_min$r.cv, nose_min$r.cv)
})

test_that("gsynth applies cv.rule, with and without se", {
  d <- .rule_panel(2)
  g_1se <- .fit_cv(d, session = 1, method = "gsynth")
  g_min <- .fit_cv(d, session = 1, method = "gsynth", cv.rule = "min")
  r_argmin <- unname(g_min$CV.out[which.min(g_min$CV.out[, "MSPE"]), "r"])
  skip_if(r_argmin == g_1se$r.cv, "this panel no longer separates the rules")
  ## "min" picks the r with the lowest CV MSPE.
  expect_equal(g_min$r.cv, r_argmin)
  ## With se = TRUE: same folds, same table, same rule.
  g_boot <- .fit_boot(d, seed = 1, method = "gsynth", cv.rule = "min",
                      cv.buffer = 2, min.T0 = 6)
  g_nose <- .fit_cv(d, session = 1, method = "gsynth", cv.rule = "min",
                    cv.buffer = 2, min.T0 = 6)
  expect_equal(g_boot$CV.out, g_nose$CV.out, tolerance = 1e-10)
  expect_equal(g_boot$r.cv, g_nose$r.cv)
})

## CFE with never-treated factors kept its own in-loop 1% rule whatever
## cv.rule was (fixed in 2.4.6). On panel 2 the lowest CV MSPE is at r = 2 and
## the 1-SE rule picks r = 1; before the fix every rule returned r = 2.
test_that("cfe with never-treated factors applies cv.rule, with and without se", {
  d <- .rule_panel(2)
  c_1se <- .fit_cv(d, session = 1, method = "cfe", time.component.from = "nevertreated")
  c_min <- .fit_cv(d, session = 1, method = "cfe", time.component.from = "nevertreated",
                   cv.rule = "min")
  r_argmin <- unname(c_min$CV.out[which.min(c_min$CV.out[, "MSPE"]), "r"])
  skip_if(r_argmin == 0, "this panel no longer has a positive lowest-MSPE r")
  ## same folds, same table; "min" picks the lowest MSPE, "1se" a smaller r
  expect_equal(c_1se$CV.out, c_min$CV.out, tolerance = 1e-10)
  expect_equal(c_min$r.cv, r_argmin)
  expect_lt(c_1se$r.cv, c_min$r.cv)
  ## with se = TRUE: same folds, same table, same rule
  c_boot <- .fit_boot(d, seed = 1, method = "cfe", time.component.from = "nevertreated",
                      cv.rule = "min")
  expect_equal(c_boot$CV.out, c_min$CV.out, tolerance = 1e-10)
  expect_equal(c_boot$r.cv, c_min$r.cv)
})
