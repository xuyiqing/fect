## fect 2.4.7 follow-ups:
##   #172  fect_cv() (method = "ife" and "mc", CV = TRUE) returns the chosen
##         model refit at the user's tol, as the gsynth path does. Before, the
##         fit made inside the CV loop at cv_tol = max(tol, 1e-3) was returned,
##         so CV = TRUE and CV = FALSE at the chosen r or lambda disagreed.
##   #174  permute = TRUE called inter_fe_mc() and inter_fe_ub() with a stale
##         argument list (no weight matrix), so every permutation errored, was
##         dropped, and the p-value was NaN ("0 permutes").
##   #175  Unscored rows of a CV table are NA, not 1e10 or 1e20, and the
##         selection rules read NA as missing. Before, a score at or above 1e9
##         (an outcome with a large scale) was read as missing too; on the
##         gsynth path the in-loop rule then never assigned r.cv and the call
##         stopped with "object 'r.cv' not found".

skip_on_cran()

.bde_data <- function(scale = 1) {
  d <- get(data("sim_gsynth", package = "fect", envir = environment()))
  d$Y <- d$Y * scale
  d
}

.bde_fit <- function(d = .bde_data(), formula = Y ~ D + X1 + X2, ...) {
  suppressWarnings(suppressMessages(
    fect::fect(formula, data = d, index = c("id", "time"),
               se = FALSE, parallel = FALSE, ...)
  ))
}

## -- #172 -------------------------------------------------------------------

test_that("#172: an mc fit with CV = TRUE equals the CV = FALSE fit at lambda.cv", {
  cv    <- .bde_fit(method = "mc", CV = TRUE, seed = 1)
  plain <- .bde_fit(method = "mc", CV = FALSE, lambda = cv$lambda.cv)
  expect_equal(cv$att.avg, plain$att.avg, tolerance = 1e-8)
  expect_equal(cv$eff, plain$eff, tolerance = 1e-8)
  expect_equal(cv$beta, plain$beta, tolerance = 1e-8)
  expect_equal(cv$att, plain$att, tolerance = 1e-8)
})

test_that("#172: an ife fit with CV = TRUE equals the CV = FALSE fit at r.cv", {
  cv    <- .bde_fit(method = "ife", CV = TRUE, r = c(0, 3), seed = 1)
  plain <- .bde_fit(method = "ife", CV = FALSE, r = cv$r.cv)
  expect_gt(cv$r.cv, 0)
  expect_equal(cv$att.avg, plain$att.avg, tolerance = 1e-8)
  expect_equal(cv$eff, plain$eff, tolerance = 1e-8)
  expect_equal(cv$beta, plain$beta, tolerance = 1e-8)
  expect_equal(cv$att, plain$att, tolerance = 1e-8)
})

test_that("#172: with tol at or above 1e-3 the CV fit is unchanged", {
  ## cv_tol = max(tol, 1e-3) equals tol here, so the loop's fit is the refit.
  cv    <- .bde_fit(method = "ife", CV = TRUE, r = c(0, 3), seed = 1, tol = 1e-3)
  plain <- .bde_fit(method = "ife", CV = FALSE, r = cv$r.cv, tol = 1e-3)
  expect_equal(cv$att.avg, plain$att.avg, tolerance = 1e-10)
})

## -- #174 -------------------------------------------------------------------

test_that("#174: permute = TRUE runs for method = 'mc', 'fe' and 'ife'", {
  for (m in c("mc", "fe", "ife")) {
    args <- list(method = m, CV = FALSE, permute = TRUE, nboots = 10, seed = 1)
    if (m == "mc")  args$lambda <- 0.05
    if (m == "ife") args$r <- 2
    perm  <- do.call(.bde_fit, args)
    plain <- do.call(.bde_fit, args[c("method", "CV", "lambda", "r")[
      c("method", "CV", "lambda", "r") %in% names(args)]])
    info <- paste("method =", m)
    expect_equal(perm$att.avg, plain$att.avg, info = info)
    expect_type(perm$permute, "list")
    expect_length(perm$permute$permute.att.avg, 10)
    expect_true(all(is.finite(perm$permute$permute.att.avg)), info = info)
    expect_true(perm$permute$p >= 0 && perm$permute$p <= 1, info = info)
  }
})

test_that("#174: the permutation fit is the estimator called with the current arguments", {
  ## one.permu() with no shuffling reproduces the plain fit's ATT.
  d <- .bde_data()
  plain <- .bde_fit(d, formula = Y ~ D, method = "mc", CV = FALSE, lambda = 0.05)
  Y <- matrix(d$Y, 30, 50); D <- matrix(d$D, 30, 50); I <- matrix(1, 30, 50)
  att <- fect:::one.permu(Y, NULL, D, I, r.cv = 0, lambda.cv = 0.05, method = "mc",
                          force = 3, tol = 1e-5, norm.para = NULL)
  expect_equal(att, abs(plain$att.avg), tolerance = 1e-8)
  att <- fect:::one.permu(Y, NULL, D, I, r.cv = 2, lambda.cv = 0.05, method = "ife",
                          force = 3, tol = 1e-5, norm.para = NULL)
  plain <- .bde_fit(d, formula = Y ~ D, method = "ife", CV = FALSE, r = 2)
  expect_equal(att, abs(plain$att.avg), tolerance = 1e-8)
})

test_that("#174: permute = TRUE stops for a method one.permu() has no estimator for", {
  ## gsynth (and cfe, ife + nevertreated) fell through the estimator branch:
  ## every permuted ATT was 0 and the p-value 0.
  expect_error(
    .bde_fit(method = "gsynth", CV = FALSE, r = 2, permute = TRUE, nboots = 5, seed = 1),
    regexp = "permute = TRUE is implemented for method = \"fe\", \"ife\" and \"mc\" \\(the fit's method is \"gsynth\"\\)")
  d <- .bde_data()
  Y <- matrix(d$Y, 30, 50); D <- matrix(d$D, 30, 50); I <- matrix(1, 30, 50)
  expect_error(
    fect:::one.permu(Y, NULL, D, I, r.cv = 2, lambda.cv = NULL, method = "gsynth",
                     force = 3, tol = 1e-5, norm.para = NULL),
    regexp = "no estimator for method = \"gsynth\"")
})

## -- #175 -------------------------------------------------------------------

test_that("#175: the gsynth CV rule chooses the same r for Y and Y * 1e5", {
  f1 <- .bde_fit(.bde_data(1),   method = "gsynth", CV = TRUE, r = c(0, 3), seed = 1)
  f2 <- .bde_fit(.bde_data(1e5), method = "gsynth", CV = TRUE, r = c(0, 3), seed = 1)
  expect_equal(unname(f1$r.cv), unname(f2$r.cv))
  expect_gt(f1$r.cv, 0)
  expect_true(all(is.finite(f1$CV.out[, "MSPE"])))
  expect_true(all(is.finite(f2$CV.out[, "MSPE"])))
  expect_gt(max(f2$CV.out[, "MSPE"]), 1e9)
  expect_equal(f2$CV.out[, "MSPE"], 1e10 * f1$CV.out[, "MSPE"], tolerance = 1e-6)
  expect_equal(f2$att.avg, 1e5 * f1$att.avg, tolerance = 1e-6)
})

test_that("#175: an ife fit with CV = TRUE chooses the same r for Y and Y * 1e5", {
  f1 <- .bde_fit(.bde_data(1),   method = "ife", CV = TRUE, r = c(0, 3), seed = 1)
  f2 <- .bde_fit(.bde_data(1e5), method = "ife", CV = TRUE, r = c(0, 3), seed = 1)
  expect_equal(unname(f1$r.cv), unname(f2$r.cv))
  expect_gt(max(f2$CV.out.ife[, "MSPE"]), 1e9)
  expect_equal(f2$CV.out.ife[, "MSPE"], 1e10 * f1$CV.out.ife[, "MSPE"], tolerance = 1e-6)
})

test_that("#175: rows the mc loop never scored are NA and lambda.cv is a scored row", {
  f <- .bde_fit(method = "mc", CV = TRUE, seed = 1, nlambda = 12)
  tab <- f$CV.out.mc
  expect_false(any(is.finite(tab) & abs(tab) >= 1e9))
  scored <- is.finite(tab[, "MSPE"])
  expect_true(any(scored))
  expect_true(all(is.na(tab[!scored, "MSPE"])))
  expect_true(scored[which(f$lambda.seq == f$lambda.cv)])
})

test_that("#175: the in-loop rule and the selection rules read NA as unscored", {
  improves <- fect:::.fect_cv_improves
  best     <- fect:::.fect_cv_best
  expect_equal(best(c(NA_real_, NA_real_)), Inf)
  expect_equal(best(c(NA_real_, 4, 5)), 4)
  expect_true(improves(c(NA_real_, NA_real_), 5e12))    # first scored row wins
  expect_true(improves(c(NA_real_, 4), 3.9))           # more than 1% better
  expect_false(improves(c(NA_real_, 4), 3.97))         # within 1%
  expect_false(improves(c(NA_real_, 4), NaN))          # a failed row never wins
  expect_false(improves(c(NA_real_, 4), Inf))
  ## a failed row (NA) in the table: the rules select among the finite rows
  means <- c(NA_real_, 2, 1.5, 1.6)
  ses   <- c(NA_real_, 0.1, 0.2, 0.1)
  expect_equal(fect:::.fect_apply_cv_rule(means, ses, rule = "min"), 3L)
  expect_equal(fect:::.fect_apply_cv_rule(means, ses, rule = "1se"), 3L)
  expect_equal(fect:::.fect_apply_cv_rule(means, ses, rule = "1pct"), 3L)
  expect_equal(fect:::.fect_apply_cv_rule(c(NA_real_, 1e12, 1e12 * 0.9), NULL, rule = "min"), 3L)
})
