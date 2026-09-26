## ---------------------------------------------------------------
## fect 2.4.7: final fixes on PR #151 (b1dded6), group Y: fect_mspe()
## refits, causal moderation (cm = TRUE), interFE()'s covariate checks.
## Every test_that block below fails on b1dded6 and passes after the fix.
## Self-contained: helpers prefixed .pfy_.
## ---------------------------------------------------------------

.pfy_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## One of fect's datasets, loaded into a local environment.
.pfy_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "fect", envir = e)
  e[[name]]
}

## A wrapper in the style of gsynth::gsynth(): it fits method = "gsynth"
## with its own defaults (force = "unit", tol = 0.001) and stores its own
## call. `inference` is an argument fect() does not have.
.pfy_wrapper <- function(formula, data, index, r = 2, force = "unit",
                         CV = FALSE, se = FALSE, inference = NULL,
                         tol = 0.001) {
  out <- fect::fect(formula, data = data, index = index, method = "gsynth",
                    r = r, CV = CV, se = se, force = force, tol = tol,
                    min.T0 = 5, parallel = FALSE)
  out$call <- match.call()
  out
}


## -- S2  fect_mspe() refits the fit's own model ----------------------------

test_that("S2: fect_mspe() refits a wrapper's fit with the wrapper, not with fect's defaults", {
  skip_on_cran()
  sg <- .pfy_data("simgsynth")
  g  <- .pfy_quiet(.pfy_wrapper(Y ~ D + X1 + X2, data = sg,
                                index = c("id", "time")))
  gi <- .pfy_quiet(.pfy_wrapper(Y ~ D + X1 + X2, data = sg,
                                index = c("id", "time"),
                                inference = "parametric"))
  ## the same model as a fect() fit
  f  <- .pfy_quiet(fect::fect(Y ~ D + X1 + X2, data = sg,
                              index = c("id", "time"), method = "gsynth",
                              r = 2, CV = FALSE, force = "unit", tol = 0.001,
                              min.T0 = 5, se = FALSE, parallel = FALSE))
  expect_identical(g$Y.ct, f$Y.ct)
  ref <- .pfy_quiet(fect::fect_mspe(f, seed = 1, k = 3))$summary
  ## b1dded6 replayed the wrapper's call into fect(): method = "fe" and
  ## force = "two-way", and stopped on `inference` ("unused argument")
  expect_identical(.pfy_quiet(fect::fect_mspe(g, seed = 1, k = 3))$summary$MSPE,
                   ref$MSPE)
  expect_identical(.pfy_quiet(fect::fect_mspe(gi, seed = 1, k = 3))$summary$MSPE,
                   ref$MSPE)
})


## -- S4  fect_mspe() evaluates the fit's call where the caller can see it -----

test_that("S4: fect_mspe() finds data and variables of the caller for fits made without a formula", {
  skip_on_cran()
  run <- function() {
    dloc <- .pfy_data("simgsynth")
    k <- 2
    fs <- .pfy_quiet(fect::fect(data = dloc, Y = "Y", D = "D",
                                X = c("X1", "X2"), index = c("id", "time"),
                                method = "ife", r = 2, CV = FALSE,
                                se = FALSE, parallel = FALSE,
                                min.T0 = k + 3))
    ff <- .pfy_quiet(fect::fect(Y ~ D + X1 + X2, data = dloc,
                                index = c("id", "time"), method = "ife",
                                r = 2, CV = FALSE, se = FALSE,
                                parallel = FALSE, min.T0 = 5))
    ## b1dded6: "object 'dloc' not found"; with the data found, min.T0 =
    ## k + 3 took fect_mspe()'s own k (3), so the refits used min.T0 = 6
    c(.pfy_quiet(fect::fect_mspe(fs, seed = 1, k = 3))$summary$MSPE,
      .pfy_quiet(fect::fect_mspe(ff, seed = 1, k = 3))$summary$MSPE)
  }
  m <- run()
  expect_identical(m[1], m[2])
  ## a formula fit made inside a function, scored outside it (b1dded6:
  ## "object 'dloc2' not found"): the formula's environment is searched
  mk <- function() {
    dloc2 <- .pfy_data("simgsynth")
    .pfy_quiet(fect::fect(Y ~ D + X1 + X2, data = dloc2,
                          index = c("id", "time"), method = "ife", r = 2,
                          CV = FALSE, se = FALSE, parallel = FALSE))
  }
  fo <- mk()
  expect_identical(.pfy_quiet(fect::fect_mspe(fo, seed = 1, k = 3))$summary$MSPE,
                   m[2])
})
