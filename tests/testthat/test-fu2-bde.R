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

