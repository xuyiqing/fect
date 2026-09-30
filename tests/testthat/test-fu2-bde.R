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

