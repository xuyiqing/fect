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


## -- C2  X together with a formula; duplicated covariate names -----------

test_that("C2: X given together with a formula stops", {
  d <- .fc_sim()
  ## 412d7ae ignored X here and fitted Y ~ D without covariates
  expect_error(
    .fc_quiet(.fc_fit(Y ~ D, data = d, X = c("X1", "X2"))),
    "Covariates were given both in the formula and in `X`", fixed = TRUE
  )
  ## an unquoted name is never looked up as a column: it stops too
  expect_error(
    .fc_quiet(.fc_fit(Y ~ D, data = d, X = X1)),
    "Covariates were given both in the formula and in `X`", fixed = TRUE
  )
  expect_error(
    .fc_quiet(fect::interFE(Y ~ D + X1, data = d, X = "X2",
                            index = c("id", "time"), r = 2)),
    "Covariates were given both in the formula and in `X`", fixed = TRUE
  )
})

test_that("C2 guard: X = NULL with a formula (gsynth's call) still runs", {
  d <- .fc_sim()
  ref <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = d))
  a <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = d, X = NULL))
  ## gsynth forwards its own X argument (default NULL) and missing Y and D
  wrapper <- function(formula, data, Y, D, X = NULL) {
    fect::fect(formula = formula, data = data, Y = Y, D = D, X = X,
               index = c("id", "time"), method = "gsynth",
               force = "two-way", r = 2, CV = FALSE, se = FALSE,
               parallel = FALSE)
  }
  b <- .fc_quiet(wrapper(Y ~ D + X1 + X2, data = d))
  expect_identical(a$att.avg, ref$att.avg)
  expect_identical(b$att.avg, ref$att.avg)
  expect_identical(b$beta, ref$beta)
})

test_that("C2: duplicated covariate names, or X naming Y or D, stop", {
  d <- .fc_sim()
  fit_x <- function(X) {
    fect::fect(data = d, Y = "Y", D = "D", X = X, index = c("id", "time"),
               method = "gsynth", force = "two-way", r = 2, CV = FALSE,
               se = FALSE, parallel = FALSE)
  }
  ## 412d7ae ran with beta = (-2.03, 2.03) for X = c("X1", "X1")
  expect_error(.fc_quiet(fit_x(c("X1", "X1"))),
               "`X` contains duplicated covariate names: \"X1\"", fixed = TRUE)
  expect_error(.fc_quiet(fit_x(c("X1", "Y"))),
               "`X` contains the outcome or the treatment variable: \"Y\"",
               fixed = TRUE)
  expect_error(.fc_quiet(fit_x(c("D", "X2"))),
               "`X` contains the outcome or the treatment variable: \"D\"",
               fixed = TRUE)
})


## -- C3  non-numeric covariates -------------------------------------------

test_that("C3: factor, character and Date covariates stop with a clear message", {
  d <- .fc_sim()
  set.seed(2)
  d$Xf <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
  d$Xc <- as.character(d$Xf)
  d$Xd <- as.Date("2000-01-01") + seq_len(nrow(d))
  ## 412d7ae: "Calling var(x) on a factor x is defunct." for the factor and
  ## a false "Variable \"Xc\" is unit-invariant." for the character column
  expect_error(.fc_quiet(.fc_fit(Y ~ D + X1 + Xf, data = d)),
               "Covariate \"Xf\" is a factor; fect() needs numeric covariates",
               fixed = TRUE)
  expect_error(.fc_quiet(.fc_fit(Y ~ D + X1 + Xc, data = d)),
               "Covariate \"Xc\" is character; fect() needs numeric covariates",
               fixed = TRUE)
  expect_error(.fc_quiet(.fc_fit(Y ~ D + X1 + Xd, data = d)),
               "Covariate \"Xd\" is of class \"Date\"", fixed = TRUE)
  expect_error(.fc_quiet(.fc_fit(Y ~ D + X1 + Xf, data = d, method = "fe",
                                 r = 0)),
               "model.matrix(~ Xf, data)", fixed = TRUE)
})

test_that("C3 guard: a logical covariate is used as 0/1", {
  d <- .fc_sim()
  d$Xl <- d$X1 > 0
  d$Xn <- as.numeric(d$Xl)
  a <- .fc_quiet(.fc_fit(Y ~ D + X2 + Xl, data = d))
  b <- .fc_quiet(.fc_fit(Y ~ D + X2 + Xn, data = d))
  expect_identical(a$att.avg, b$att.avg)
  expect_identical(unname(a$beta), unname(b$beta))
})


## -- C4  time index given as factor or character --------------------------

test_that("C4: a factor or character time index in numeric order fits like the numbers", {
  d <- .fc_sim()
  ref <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = d))
  ## 412d7ae sorted both as text ("1", "10", ...) and gsynth stopped with
  ## "Gsynth can't be used when treatments have reversals."
  df <- d
  df$time <- factor(df$time)
  a <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = df))
  dc <- d
  dc$time <- as.character(dc$time)
  b <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = dc))
  for (o in list(a, b)) {
    expect_identical(o$att.avg, ref$att.avg)
    expect_identical(unname(o$eff), unname(ref$eff))
    expect_identical(o$rawtime, ref$rawtime)
  }
  ## levels that are increasing numbers become those numbers
  dy <- d
  dy$time <- factor(dy$time + 2000)
  y <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = dy))
  expect_identical(y$att.avg, ref$att.avg)
  expect_identical(y$rawtime, ref$rawtime + 2000)
})

test_that("C4: a factor with non-numeric levels is read in level order", {
  d <- .fc_sim()
  ref <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = d))
  lev <- paste0("2018-", 1:30)               # not in time order as text
  dl <- d
  dl$time <- factor(paste0("2018-", dl$time), levels = lev)
  o <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = dl))
  expect_identical(o$att.avg, ref$att.avg)
  expect_identical(unname(o$eff), unname(ref$eff))
  expect_identical(o$rawtime, lev)
  expect_identical(rownames(o$eff), lev)
  expect_identical(levels(o$data.long$time), lev)
  expect_identical(as.character(o$data.long$time), as.character(dl$time))
})

test_that("C4: ife with a factor time index keeps the event-time axis", {
  d <- .fc_sim()
  ref <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = d, method = "ife"))
  df <- d
  df$time <- factor(df$time)
  o <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = df, method = "ife"))
  ## 412d7ae: 22 event times instead of 30
  expect_identical(o$time, ref$time)
  expect_identical(o$att, ref$att)
  expect_identical(unname(o$eff), unname(ref$eff))
})

test_that("C4: a character time index with non-numbers stops", {
  d <- .fc_sim()
  d$time <- sprintf("t%02d", d$time)
  ## 412d7ae sorted these as text (here in time order) and ran
  expect_error(
    .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = d, method = "fe", r = 0)),
    "The time index \"time\" is character and some values are not numbers (for example \"t01\")",
    fixed = TRUE
  )
})

test_that("C4 guard: Date and unbalanced factor indices with unused levels", {
  d <- .fc_sim()
  ref <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = d))
  dd <- d
  dd$time <- as.Date("2000-01-01") + 7 * dd$time
  o <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = dd))
  expect_identical(o$att.avg, ref$att.avg)
  expect_s3_class(o$rawtime, "Date")
  ## the original gsynth #13 setting: unbalanced panel, factor indices with
  ## unused levels
  ub <- d[!(d$id >= 140 & d$time <= 3), ]
  ub.ref <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = ub))
  uf <- ub
  uf$time <- factor(uf$time + 2000, levels = 1990:2040)
  uf$id <- factor(uf$id, levels = 90:160)
  o2 <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = uf))
  expect_identical(o2$att.avg, ub.ref$att.avg)
})

test_that("C4 guard: the dropped-period message shows the level label", {
  d <- .fc_sim()
  ## no control unit is observed in period 25, so fect drops that period
  d <- d[!(d$time == 25 & d$id > 105), ]
  d$time <- factor(paste0("w", d$time), levels = paste0("w", 1:30))
  expect_message(
    .fc_fit(Y ~ D + X1 + X2, data = d, method = "fe", r = 0),
    "There are not any observations under control at w25", fixed = TRUE
  )
})
