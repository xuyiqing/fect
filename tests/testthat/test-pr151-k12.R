## ---------------------------------------------------------------
## fect 2.4.7: fect_mspe() and r.cv.rolling() score the model's full
## prediction at the held-out cells (fect #163). Every test_that block
## below fails on c4331a4 and passes after the fix.
## Self-contained: helpers prefixed .k12_.
##
## sim_gsynth is made with two-way fixed effects, two factors, X1
## (coefficient 1) and X2 (coefficient 3); var(sim_gsynth$error) = 0.991.
## ---------------------------------------------------------------

.k12_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## One of fect's datasets, loaded into a local environment.
.k12_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "fect", envir = e)
  e[[name]]
}

## The fits below are made in the test body with `data = d`, so that
## fect_mspe() finds `d` where it is called (a fit made inside a helper
## with `data = data` would make fect_mspe() look up the name `data`).


## -- a correctly specified model scores near the error variance ------------

test_that("K12: a correctly specified model scores near the error variance", {
  skip_on_cran()
  d <- .k12_data("sim_gsynth")
  s2 <- stats::var(d$error)
  ife <- .k12_quiet(fect::fect(Y ~ D + X1 + X2, data = d,
    index = c("id", "time"), CV = FALSE, se = FALSE, parallel = FALSE,
    method = "ife", force = "two-way", r = 2))
  gsy <- .k12_quiet(fect::fect(Y ~ D + X1 + X2, data = d,
    index = c("id", "time"), CV = FALSE, se = FALSE, parallel = FALSE,
    method = "gsynth", force = "two-way", r = 2))
  score <- function(fits) {
    .k12_quiet(fect::fect_mspe(fits, seed = 1234, cv.method = "block",
                               k = 10))$summary
  }
  ## not-yet-treated fit with covariates (c4331a4: 44.3, X times beta left out)
  s.ife <- score(ife)
  expect_lt(abs(s.ife$MSPE - s2), 0.5)
  ## never-treated fit, held-out control cells (c4331a4: about 124, only
  ## the factor part was scored)
  s.gsy <- score(gsy)
  expect_lt(abs(s.gsy$MSPE - s2), 0.5)
  ## never-treated fit scored on the not-yet-treated fit's cells, which
  ## include treated units' pre-treatment cells (c4331a4: about 124)
  s.both <- score(list(ife = ife, gsynth = gsy))
  expect_lt(abs(s.both$MSPE[s.both$Model == "gsynth"] - s2), 0.5)
  expect_identical(s.both$Hidden_N[s.both$Model == "gsynth"],
                   s.both$Hidden_N[s.both$Model == "ife"])
  ## the same model with normalize = TRUE: covariates enter on the
  ## outcome's scale
  ife.n <- .k12_quiet(fect::fect(Y ~ D + X1 + X2, data = d,
    index = c("id", "time"), CV = FALSE, se = FALSE, parallel = FALSE,
    method = "ife", force = "two-way", r = 2, normalize = TRUE))
  s.n <- score(ife.n)
  expect_lt(abs(s.n$MSPE - s.ife$MSPE), 0.05 * s.ife$MSPE)
})


## -- book chapter 06 section 6.5.2: never-treated CFE models no longer tie --

test_that("K12: never-treated CFE models with different specifications score differently", {
  skip_on_cran()
  d <- .k12_data("sim_gsynth")
  m1 <- .k12_quiet(fect::fect(Y ~ D + X1 + X2, data = d,
    index = c("id", "time"), CV = FALSE, se = FALSE, parallel = FALSE,
    method = "gsynth", force = "two-way", r = 2, max.iteration = 20000))
  m2 <- .k12_quiet(fect::fect(Y ~ D + X1 + X2, data = d,
    index = c("id", "time"), CV = FALSE, se = FALSE, parallel = FALSE,
    method = "cfe", force = "two-way", time.component.from = "nevertreated",
    r = 2, max.iteration = 20000))
  m3 <- .k12_quiet(fect::fect(Y ~ D + X1 + X2, data = d,
    index = c("id", "time"), CV = FALSE, se = FALSE, parallel = FALSE,
    method = "cfe", force = "two-way", time.component.from = "nevertreated",
    Q.type = "linear", r = 2, max.iteration = 20000))
  m4 <- .k12_quiet(fect::fect(Y ~ D, data = d,
    index = c("id", "time"), CV = FALSE, se = FALSE, parallel = FALSE,
    method = "cfe", force = "two-way", time.component.from = "nevertreated",
    r = 0, max.iteration = 20000))
  s <- .k12_quiet(fect::fect_mspe(list(gsynth_r2 = m1, CFE_r2 = m2,
                                       CFE_linear_r2 = m3, CFE_r0_noX = m4),
                                  seed = 1234, k = 10))$summary
  mspe <- stats::setNames(s$MSPE, s$Model)
  ## c4331a4: the three CFE models all scored 110.33 (a prediction of 0)
  expect_lt(mspe[["CFE_r2"]], 3)
  expect_lt(abs(mspe[["CFE_r2"]] - mspe[["gsynth_r2"]]),
            0.05 * mspe[["gsynth_r2"]])
  expect_lt(mspe[["CFE_r2"]], mspe[["CFE_linear_r2"]])
  expect_lt(mspe[["CFE_linear_r2"]], mspe[["CFE_r0_noX"]])
  expect_gt(mspe[["CFE_r0_noX"]], 5 * mspe[["CFE_r2"]])
})


## -- r.cv.rolling(method = "gsynth") recovers the two factors ---------------

test_that("K12: r.cv.rolling(method = 'gsynth') picks r = 2 on data made with two factors", {
  skip_on_cran()
  d <- .k12_data("sim_gsynth")
  cv <- .k12_quiet(fect::r.cv.rolling(Y ~ D + X1 + X2, data = d,
                                      index = c("id", "time"),
                                      method = "gsynth", r.max = 3,
                                      min.T0 = 10, force = "two-way",
                                      seed = 1, verbose = FALSE))
  ## c4331a4: r.cv = 0, MSPE about 130 for every r
  expect_identical(cv$r.cv, 2L)
  expect_lt(cv$mspe$mspe[cv$mspe$r == 2], 2)
})


## -- covariates lower the score of a model that needs them -----------------

test_that("K12: adding the true covariates lowers the score", {
  skip_on_cran()
  d <- .k12_data("sim_gsynth")
  with.x <- .k12_quiet(fect::fect(Y ~ D + X1 + X2, data = d,
    index = c("id", "time"), CV = FALSE, se = FALSE, parallel = FALSE,
    method = "ife", force = "two-way", r = 2))
  no.x <- .k12_quiet(fect::fect(Y ~ D, data = d,
    index = c("id", "time"), CV = FALSE, se = FALSE, parallel = FALSE,
    method = "ife", force = "two-way", r = 2))
  s <- .k12_quiet(fect::fect_mspe(list(with_X = with.x, no_X = no.x),
                                  seed = 1234))$summary
  mspe <- stats::setNames(s$MSPE, s$Model)
  ## c4331a4: 45.10 with the covariates, 56.48 without
  expect_lt(mspe[["with_X"]], 0.2 * mspe[["no_X"]])
})
