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


## -- S3  causal moderation: est.cm ------------------------------------------

test_that("S3: cm = TRUE stores est.cm without SEs and with cross-validation", {
  skip_on_cran()
  d <- .pfy_data("simdata")
  cmfit <- function(...) {
    .pfy_quiet(fect::fect(Y ~ D + X1 + X2, data = d, index = c("id", "time"),
                          force = "two-way", parallel = FALSE, cm = TRUE, ...))
  }
  f0 <- cmfit(method = "fe", se = FALSE)
  f1 <- cmfit(method = "fe", se = TRUE, nboots = 5, seed = 1)
  ## b1dded6: no est.cm without SEs
  expect_false(is.null(f0$est.cm))
  expect_identical(f0$est.cm, f1$est.cm)
  i0 <- cmfit(method = "ife", r = 2, CV = FALSE, se = FALSE)
  i1 <- cmfit(method = "ife", r = 2, CV = FALSE, se = TRUE, nboots = 5,
              seed = 1)
  expect_identical(i0$est.cm, i1$est.cm)
  ## b1dded6: no est.cm after cross-validation
  cv <- cmfit(method = "ife", r = c(0, 3), CV = TRUE, se = FALSE, seed = 1)
  ir <- cmfit(method = "ife", r = as.numeric(cv$r.cv), CV = FALSE, se = FALSE)
  expect_false(is.null(cv$est.cm))
  expect_identical(cv$est.cm, ir$est.cm)
  ## the over-identification test runs without SEs, as with them
  t0 <- .pfy_quiet(fect::fect_iden(f0, moderator = "X1"))
  t1 <- .pfy_quiet(fect::fect_iden(f1, moderator = "X1"))
  expect_identical(t0$e1$stat, t1$e1$stat)
})

test_that("S3: cm = TRUE with never-treated controls stops", {
  sg <- .pfy_data("simgsynth")
  ## b1dded6: ran, without est.cm and without a message
  expect_error(
    .pfy_quiet(fect::fect(Y ~ D, data = sg, index = c("id", "time"),
                          method = "ife", r = 1, CV = FALSE,
                          time.component.from = "nevertreated", cm = TRUE,
                          se = FALSE, parallel = FALSE)),
    "time.component.from = \"notyettreated\"", fixed = TRUE)
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


## -- S5  interFE(): covariates absorbed by the model, for every force --------

test_that("S5: interFE() stops with the reason for every force value", {
  d <- withr::with_seed(7, {
    N <- 30
    TT <- 10
    d <- expand.grid(time = seq_len(TT), id = seq_len(N))
    a <- stats::rnorm(N)
    g <- stats::rnorm(TT)
    d$X1 <- stats::rnorm(nrow(d))
    d$Zu <- a[d$id]              # constant over time within units
    d$Zc <- 1.5                  # constant everywhere
    d$Zut <- a[d$id] + g[d$time] # unit-level plus period-level
    d$Y <- 1 + 0.5 * d$X1 + a[d$id] + g[d$time] + stats::rnorm(nrow(d))
    d
  })
  fit <- function(f, force) {
    fect::interFE(f, data = d, index = c("id", "time"), force = force, r = 0)
  }
  ## b1dded6: an NaN coefficient, without a message
  expect_error(fit(Y ~ X1 + Zc, "none"),
               "Variable \"Zc\" does not vary (it is absorbed by the intercept). Remove it.",
               fixed = TRUE)
  expect_error(fit(Y ~ X1 + Zut, "two-way"),
               "Variable \"Zut\" is the sum of a unit-level and a period-level variable (it is absorbed by the unit and time fixed effects). Remove it.",
               fixed = TRUE)
  ## estimable under the other fixed effects: unchanged
  expect_true(all(is.finite(fit(Y ~ X1 + Zut, "unit")$beta)))
  expect_true(all(is.finite(fit(Y ~ X1 + Zu, "none")$beta)))
  expect_error(fit(Y ~ X1 + Zu, "two-way"),
               "Variable \"Zu\" does not vary over time within units (it is absorbed by the unit fixed effects). Remove it.",
               fixed = TRUE)
})
