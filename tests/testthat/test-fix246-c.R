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


## -- C5 + C6  collinear and FE-absorbed covariates ------------------------

## capture the warnings of a call, silencing messages
.fc_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(
    expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")
  )
  list(value = val, warnings = w)
}

test_that("C5: an exactly collinear covariate is dropped with a warning (V10)", {
  d <- .fc_sim()
  d$X3 <- 2 * d$X1
  d$X4 <- d$X1 + d$X2
  for (r in c(0, 2)) {
    ref <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = d, r = r))
    ## 412d7ae: ATT 5.0841 -> 7.1057 (r = 0) and 5.5435 -> 6.5098 (r = 2)
    ## with X3 = 2 * X1; validX stayed 1 and nothing was reported
    for (z in c("X3", "X4")) {
      f <- stats::reformulate(c("D", "X1", "X2", z), response = "Y")
      o <- .fc_warnings(.fc_fit(f, data = d, r = r))
      expect_length(o$warnings, 1)
      expect_match(o$warnings,
                   paste0("Dropped 1 covariate that cannot be estimated on the cells used to fit the model: \"",
                          z, "\" is a linear combination of other covariates"),
                   fixed = TRUE)
      fit <- o$value
      expect_identical(fit$att.avg, ref$att.avg)
      expect_identical(unname(fit$eff), unname(ref$eff))
      expect_identical(fit$X, c("X1", "X2", z))
      expect_identical(rownames(fit$beta), c("X1", "X2", z))
      expect_identical(unname(fit$beta[1:2, 1]), unname(ref$beta[, 1]))
      expect_true(is.na(fit$beta[3, 1]))
      expect_equal(fit$validX, 1)
    }
  }
})

test_that("C5: collinear covariates no longer crash ife, fe, mc and cfe", {
  d <- .fc_sim()
  d$X3 <- 2 * d$X1
  fits <- list(
    ife = function(f) .fc_fit(f, data = d, method = "ife"),
    fe  = function(f) .fc_fit(f, data = d, method = "fe", r = 0),
    mc  = function(f) fect::fect(f, data = d, index = c("id", "time"),
                                 method = "mc", lambda = 0.1, CV = FALSE,
                                 se = FALSE, parallel = FALSE),
    cfe = function(f) .fc_fit(f, data = d, method = "cfe", r = 0)
  )
  for (m in names(fits)) {
    ref <- .fc_quiet(fits[[m]](Y ~ D + X1 + X2))
    ## 412d7ae: "inv(): matrix is singular" for ife, fe and mc
    o <- .fc_warnings(fits[[m]](Y ~ D + X1 + X2 + X3))
    expect_match(o$warnings, "\"X3\" is a linear combination", fixed = TRUE)
    expect_identical(o$value$att.avg, ref$att.avg)
    expect_identical(unname(o$value$eff), unname(ref$eff))
  }
  ## like lm(), the earlier column of a collinear set is kept
  ref <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = d, method = "fe", r = 0))
  o <- .fc_warnings(.fc_fit(Y ~ D + X3 + X1 + X2, data = d, method = "fe",
                            r = 0))
  expect_match(o$warnings, "\"X1\" is a linear combination", fixed = TRUE)
  expect_identical(o$value$att.avg, ref$att.avg)
  expect_identical(rownames(o$value$beta), c("X3", "X1", "X2"))
  expect_true(is.na(o$value$beta["X1", 1]))
  expect_equal(o$value$beta["X3", 1], ref$beta["X1", 1] / 2)
})

test_that("C6: covariates absorbed by the fixed effects are dropped with correct labels", {
  d <- .fc_sim()
  d$Zu <- stats::ave(d$X1, d$id)     # constant over time within each unit
  d$Zt <- stats::ave(d$X1, d$time)   # constant across units in each period
  ref <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = d, method = "fe", r = 0))
  ## 412d7ae stopped with swapped labels: Zu "is unit-invariant", Zt "is
  ## time-invariant", while ?fect said such covariates are dropped
  o <- .fc_warnings(.fc_fit(Y ~ D + X1 + X2 + Zu, data = d, method = "fe",
                            r = 0))
  expect_match(o$warnings,
               "\"Zu\" does not vary over time within units, so it is absorbed by the unit fixed effects",
               fixed = TRUE)
  expect_identical(o$value$att.avg, ref$att.avg)
  o <- .fc_warnings(.fc_fit(Y ~ D + X1 + X2 + Zt, data = d, method = "fe",
                            r = 0))
  expect_match(o$warnings,
               "\"Zt\" does not vary across units within periods, so it is absorbed by the time fixed effects",
               fixed = TRUE)
  expect_identical(o$value$att.avg, ref$att.avg)
  ## the same check under never-treated fitting (gsynth)
  refg <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = d))
  o <- .fc_warnings(.fc_fit(Y ~ D + X1 + Zu + X2, data = d))
  expect_match(o$warnings, "\"Zu\" does not vary over time within units",
               fixed = TRUE)
  expect_identical(o$value$att.avg, refg$att.avg)
  expect_identical(rownames(o$value$beta), c("X1", "Zu", "X2"))
})

test_that("C6: when every covariate is dropped the fit has no covariates", {
  d <- .fc_sim()
  d$Zu <- stats::ave(d$X1, d$id)
  d$Zt <- stats::ave(d$X1, d$time)
  ref <- .fc_quiet(.fc_fit(Y ~ D, data = d, method = "fe", r = 0))
  o <- .fc_warnings(.fc_fit(Y ~ D + Zu + Zt, data = d, method = "fe", r = 0))
  expect_match(o$warnings, "^Dropped 2 covariates that cannot be estimated")
  expect_match(o$warnings, "Their coefficients are reported as NA.",
               fixed = TRUE)
  fit <- o$value
  expect_identical(fit$att.avg, ref$att.avg)
  expect_identical(unname(fit$eff), unname(ref$eff))
  expect_equal(fit$validX, 0)
  expect_identical(fit$X, c("Zu", "Zt"))
  expect_identical(dim(fit$beta), c(2L, 1L))
  expect_true(all(is.na(fit$beta)))
})

test_that("C5: a covariate with no variation on the estimation cells is dropped with a warning", {
  d <- .fc_sim()
  set.seed(3)
  ## zero for every never-treated unit: gsynth estimates beta on those only
  d$Ztr <- ifelse(d$D == 1, stats::rnorm(nrow(d)), 0)
  ref <- .fc_quiet(.fc_fit(Y ~ D + X1 + X2, data = d))
  ## 412d7ae shed it silently (same numbers, beta NA, no message)
  o <- .fc_warnings(.fc_fit(Y ~ D + X1 + X2 + Ztr, data = d))
  expect_match(o$warnings,
               "\"Ztr\" has no variation on the cells used to estimate the covariate coefficients",
               fixed = TRUE)
  expect_identical(o$value$att.avg, ref$att.avg)
  expect_identical(unname(o$value$eff), unname(ref$eff))
})

test_that("C5: a covariate absorbed by an extra fixed effect (cfe) is dropped", {
  d <- .fc_sim()
  d$rt <- paste((d$id %% 5) + 1, d$time)      # region x period effects
  d$Zrt <- stats::ave(d$X1, d$rt)             # constant within region-period
  ref <- .fc_quiet(fect::fect(Y ~ D + X1, data = d,
                              index = c("id", "time", "rt"), method = "cfe",
                              r = 0, CV = FALSE, se = FALSE, parallel = FALSE))
  ## 412d7ae fitted it anyway (beta -0.467, a different ATT)
  o <- .fc_warnings(fect::fect(Y ~ D + X1 + Zrt, data = d,
                               index = c("id", "time", "rt"), method = "cfe",
                               r = 0, CV = FALSE, se = FALSE,
                               parallel = FALSE))
  expect_match(o$warnings, "\"Zrt\" is absorbed by the fixed effects",
               fixed = TRUE)
  expect_identical(o$value$att.avg, ref$att.avg)
})

test_that("C5: standard errors with a dropped covariate: NA rows, same inference", {
  skip_on_cran()
  d <- .fc_sim()
  d$X3 <- 2 * d$X1
  fit_se <- function(f) {
    fect::fect(f, data = d, index = c("id", "time"), method = "fe",
               force = "two-way", se = TRUE, nboots = 20, seed = 11,
               parallel = FALSE)
  }
  ref <- .fc_quiet(fit_se(Y ~ D + X1 + X2))
  ## 412d7ae: "inv(): matrix is singular"
  fit <- .fc_quiet(fit_se(Y ~ D + X1 + X3 + X2))
  expect_identical(fit$est.avg, ref$est.avg)
  expect_identical(fit$att.avg.boot, ref$att.avg.boot)
  expect_identical(rownames(fit$est.beta), c("X1", "X3", "X2"))
  expect_true(all(is.na(fit$est.beta["X3", ])))
  expect_identical(unname(fit$est.beta[c("X1", "X2"), ]),
                   unname(ref$est.beta))
  expect_identical(dim(fit$beta.boot), c(3L, ncol(ref$beta.boot)))
  expect_true(all(is.na(fit$beta.boot[2, ])))
  expect_identical(fit$beta.boot[c(1, 3), ], ref$beta.boot)
})

test_that("C5: the C++ (X'X)^-1 falls back to a pseudo-inverse when singular", {
  set.seed(4)
  x <- array(stats::rnorm(60), c(10, 3, 2))
  x[, , 2] <- 2 * x[, , 1]
  ## 412d7ae: "inv(): matrix is singular"
  a <- fect:::XXinv(x)
  b <- fect:::wXXinv(x, matrix(1, 10, 3))
  xx <- matrix(c(sum(x[, , 1]^2), rep(sum(x[, , 1] * x[, , 2]), 2),
                 sum(x[, , 2]^2)), 2)
  expect_true(all(is.finite(a)))
  expect_equal(a, MASS::ginv(xx), tolerance = 1e-8)
  expect_equal(b, MASS::ginv(xx), tolerance = 1e-8)
  ## guard: a well-conditioned X'X keeps the exact inverse
  y <- array(stats::rnorm(60), c(10, 3, 2))
  yy <- matrix(c(sum(y[, , 1]^2), rep(sum(y[, , 1] * y[, , 2]), 2),
                 sum(y[, , 2]^2)), 2)
  expect_equal(fect:::XXinv(y), solve(yy), tolerance = 1e-12)
})

test_that("C5 guard: full-rank covariates give no new warning", {
  d <- .fc_sim()
  for (m in c("fe", "gsynth", "ife")) {
    expect_no_warning(
      suppressMessages(.fc_fit(Y ~ D + X1 + X2, data = d, method = m,
                               r = if (m == "fe") 0 else 2))
    )
  }
})
