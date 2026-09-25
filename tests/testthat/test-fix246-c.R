## ---------------------------------------------------------------
## fix246 section C: input handling (items C1-C6).
##
## Regression tests for the 2.4.6 correctness run. Each block fails on
## fect dev @ 412d7ae and passes after the fix, except the blocks marked
## "guard", which pin behaviour that must not change. Fits are small and
## serial (parallel = FALSE, se = FALSE unless the item is about SEs).
## ---------------------------------------------------------------

.fc_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## gsynth's simulated panel: 50 units x 30 periods, units 101-105 treated
## from period 21, 45 never-treated units, covariates X1 and X2.
.fc_sim <- function() {
  e <- new.env()
  utils::data("simgsynth", package = "fect", envir = e)
  e$simgsynth
}

.fc_fit <- function(formula, data = .fc_sim(), method = "gsynth", r = 2, ...) {
  fect::fect(formula, data = data, index = c("id", "time"), method = method,
             force = "two-way", r = r, CV = FALSE, se = FALSE,
             parallel = FALSE, ...)
}


## -- C1  formula terms must be bare column names ------------------------

test_that("C1: a transformed outcome stops instead of fitting the raw column", {
  d <- .fc_sim()
  d$Ypos <- d$Y + 20
  ## 412d7ae fitted log(Ypos) ~ D as Ypos ~ D
  expect_error(
    .fc_quiet(.fc_fit(log(Ypos) ~ D, data = d)),
    "fect() formulas take bare column names only; not a column name: `log(Ypos)`",
    fixed = TRUE
  )
})

test_that("C1: factor(), interactions, `.` and other expressions stop", {
  d <- .fc_sim()
  set.seed(1)
  d$g3 <- sample(1:3, nrow(d), replace = TRUE)
  ## 412d7ae fitted factor(g3) as one linear slope and X1 * X2 as X1 + X2
  expect_error(.fc_quiet(.fc_fit(Y ~ D + X1 + factor(g3), data = d)),
               "not a column name: `factor(g3)`", fixed = TRUE)
  expect_error(.fc_quiet(.fc_fit(Y ~ D + X1 * X2, data = d)),
               "not a column name: `X1 * X2`", fixed = TRUE)
  expect_error(.fc_quiet(.fc_fit(Y ~ D + X1:X2, data = d)),
               "not a column name: `X1:X2`", fixed = TRUE)
  expect_error(.fc_quiet(.fc_fit(Y ~ ., data = d)),
               "not a column name: `.`", fixed = TRUE)
  expect_error(.fc_quiet(.fc_fit(Y ~ D + X1 - X2, data = d)),
               "bare column names only", fixed = TRUE)
  expect_error(.fc_quiet(.fc_fit(Y ~ D + X1 + 2, data = d)),
               "not a column name: `2`", fixed = TRUE)
  ## every bad term is named in one message, which says what to do
  err <- tryCatch(.fc_fit(log(Y + 20) ~ D + I(X1^2) + X2, data = d),
                  error = conditionMessage)
  expect_match(err, "not column names: `log(Y + 20)`, `I(X1^2)`", fixed = TRUE)
  expect_match(err, "Create the variable first", fixed = TRUE)
  expect_match(err, "model.matrix()", fixed = TRUE)
})

test_that("C1: the outcome on the right-hand side, or no treatment, stops", {
  ## 412d7ae dropped the repeated outcome silently (all.vars() deduplicates)
  expect_error(.fc_quiet(.fc_fit(Y ~ D + X1 + Y)),
               "The outcome \"Y\" also appears on the right-hand side",
               fixed = TRUE)
  expect_error(.fc_quiet(.fc_fit(Y ~ 1)),
               "needs a treatment variable on the right-hand side", fixed = TRUE)
})

test_that("C1 guard: intercept specifiers are accepted and ignored", {
  d <- .fc_sim()
  ref <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = d))
  forms <- list(Y ~ D + X1 + X2 - 1, Y ~ 0 + D + X1 + X2, Y ~ D + X1 + X2 + 0,
                Y ~ 1 + D + X1 + X2, Y ~ -1 + D + X1 + X2,
                Y ~ D + X1 + X2 - 0, Y ~ D + X1 + X2 + X1)
  for (f in forms) {
    o <- .fc_quiet(.fc_fit(f, data = d))
    expect_identical(o$att.avg, ref$att.avg)
    expect_identical(o$eff, ref$eff)
    expect_identical(o$beta, ref$beta)
    expect_identical(o$X, c("X1", "X2"))
  }
})

test_that("C1: interFE() stops on D:time instead of splitting it into D and time", {
  d <- .fc_sim()
  ## 412d7ae read D:time as the two regressors D and time
  expect_error(
    .fc_quiet(fect::interFE(Y ~ X1 + X2 + D:time, data = d,
                            index = c("id", "time"), r = 2,
                            force = "two-way")),
    "interFE() formulas take bare column names only; not a column name: `D:time`",
    fixed = TRUE
  )
  ## the documented workaround (build the column first) runs, and an
  ## intercept specifier changes nothing
  d$Dtime <- d$D * d$time
  a <- .fc_quiet(fect::interFE(Y ~ D + Dtime + X1 + X2, data = d,
                               index = c("id", "time"), r = 2,
                               force = "two-way"))
  b <- .fc_quiet(fect::interFE(Y ~ 0 + D + Dtime + X1 + X2, data = d,
                               index = c("id", "time"), r = 2,
                               force = "two-way"))
  expect_identical(a$beta, b$beta)
  expect_length(a$beta, 4)
})
