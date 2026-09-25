## ---------------------------------------------------------------
## fect 2.4.6 correctness fixes, group C (input handling) and the
## loo / routing items B8b-B8e. Run 2026-09-24-fix246-correctness; one
## block (or more) per item. Every item has at least one block that fails
## on dev @ 412d7ae and passes after the fix.
## ---------------------------------------------------------------

.fix246c_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## Collect the messages of `expr` (warnings are muffled); errors propagate.
.fix246c_messages <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(
    suppressWarnings(expr),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(value = val, messages = msgs)
}

## Collect the warnings of `expr` (messages are muffled); errors propagate.
.fix246c_warnings <- function(expr) {
  warns <- character(0)
  val <- withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = val, warnings = warns)
}

## Balanced panel: `Ntr` units treated from period T0 + 1 on, `Nco`
## never-treated units, two factors, two covariates.
.fix246c_panel <- function(Nco = 20, Ntr = 5, TT = 12, T0 = 8, seed = 7) {
  set.seed(seed)
  N <- Ntr + Nco
  d <- expand.grid(time = seq_len(TT), id = seq_len(N))[, c("id", "time")]
  F <- matrix(stats::rnorm(TT * 2), TT, 2)
  L <- matrix(stats::rnorm(N * 2), N, 2)
  d$X1 <- stats::rnorm(nrow(d))
  d$X2 <- stats::rnorm(nrow(d))
  d$D <- as.numeric(d$id <= Ntr & d$time > T0)
  d$Y <- 1 + stats::rnorm(N)[d$id] + stats::rnorm(TT)[d$time] +
    rowSums(F[d$time, ] * L[d$id, ]) + 0.5 * d$X1 + 0.3 * d$X2 + 2 * d$D +
    stats::rnorm(nrow(d), sd = 0.3)
  d
}


## -- B8b  loo refits keep loading.bound ------------------------------------

test_that("B8b: loo pre-trend refits use the simplex bound of the fit", {
  skip_on_cran()
  d <- .fix246c_panel()
  fit_loo <- function(...) .fix246c_quiet(fect::fect(
    Y ~ D + X1, data = d, index = c("id", "time"), method = "gsynth",
    r = 2, CV = FALSE, loo = TRUE, nboots = 5, parallel = FALSE, seed = 1, ...
  ))
  f.s <- fit_loo(loading.bound = "simplex", gamma.loading = 1)
  f.n <- fit_loo()
  expect_identical(f.s$loading.bound, "simplex")
  ## Before 2.4.6 every refit ran unbounded, so the two sets of loo
  ## pre-trend estimates were identical.
  expect_false(isTRUE(all.equal(f.s$pre.est.att[, "ATT"], f.n$pre.est.att[, "ATT"])))
  ## Each loo estimate is the bounded model refit with that period held out.
  for (kk in c(0, -3)) {
    pl <- .fix246c_quiet(fect::fect(
      Y ~ D + X1, data = d, index = c("id", "time"), method = "gsynth",
      r = 2, CV = FALSE, loading.bound = "simplex", gamma.loading = 1,
      placeboTest = TRUE, placebo.period = kk, se = FALSE, parallel = FALSE
    ))
    expect_equal(unname(f.s$pre.est.att[as.character(kk), "ATT"]),
                 pl$att[which(pl$time == kk)], tolerance = 1e-10, info = paste("kk =", kk))
  }
})


## -- B8c  loo refits keep time.component.from ------------------------------

test_that("B8c: loo refits of a never-treated cfe fit use the never-treated estimator", {
  skip_on_cran()
  d <- .fix246c_panel()
  f <- .fix246c_quiet(fect::fect(
    Y ~ D + X1, data = d, index = c("id", "time"), method = "cfe",
    time.component.from = "nevertreated", r = 2, CV = FALSE, loo = TRUE,
    nboots = 5, parallel = FALSE, seed = 1
  ))
  ## Each loo estimate equals the same model refit with that period held
  ## out. Before 2.4.6 the refits ran the not-yet-treated cfe estimator
  ## (differences up to 0.69 on this panel).
  for (kk in c(0, -1, -5)) {
    pl <- .fix246c_quiet(fect::fect(
      Y ~ D + X1, data = d, index = c("id", "time"), method = "cfe",
      time.component.from = "nevertreated", r = 2, CV = FALSE,
      placeboTest = TRUE, placebo.period = kk, se = FALSE, parallel = FALSE
    ))
    expect_equal(unname(f$pre.est.att[as.character(kk), "ATT"]),
                 pl$att[which(pl$time == kk)], tolerance = 1e-10, info = paste("kk =", kk))
  }
})


## -- B8d  loo refits keep para.error -----------------------------------------

test_that("B8d: parametric loo refits use the fit's para.error", {
  skip_on_cran()
  d <- .fix246c_panel()
  res <- .fix246c_messages(fect::fect(
    Y ~ D + X1, data = d, index = c("id", "time"), method = "gsynth",
    r = 2, CV = FALSE, loo = TRUE, se = TRUE, vartype = "parametric",
    para.error = "wild", nboots = 5, parallel = FALSE, seed = 1
  ))
  used <- regmatches(res$messages,
                     regexpr("para.error = \"[a-z]+\"", res$messages))
  ## one message for the main fit and one per refit (8 pre-periods);
  ## before 2.4.6 the refits reported "empirical" (the "auto" choice)
  expect_length(used, 9L)
  expect_true(all(used == "para.error = \"wild\""))
})


## -- B8e  ife + never-treated with se = TRUE, CV = FALSE -------------------

test_that("B8e: ife + nevertreated gives the same point estimate with and without SEs", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  fit <- function(...) .fix246c_quiet(fect::fect(
    Y ~ D + X1 + X2, data = simgsynth, index = c("id", "time"),
    force = "two-way", r = 2, CV = FALSE, parallel = FALSE, ...
  ))
  a <- fit(method = "ife", time.component.from = "nevertreated", se = FALSE)
  g <- fit(method = "gsynth", se = FALSE)
  expect_equal(a$att.avg, g$att.avg, tolerance = 1e-12)
  ## Before 2.4.6 the se = TRUE fits ran the not-yet-treated ife estimator
  ## (att.avg 5.57884 instead of 5.54329).
  for (vt in c("bootstrap", "parametric", "jackknife")) {
    b <- fit(method = "ife", time.component.from = "nevertreated", se = TRUE,
             vartype = vt, nboots = 10, seed = 1)
    expect_equal(b$att.avg, a$att.avg, tolerance = 1e-12, info = vt)
    expect_equal(b$eff, a$eff, tolerance = 1e-10, info = vt)
    expect_identical(b$method, "gsynth", info = vt)
  }
  ## and the bootstrap is the gsynth bootstrap (same draws, same SE)
  b <- fit(method = "ife", time.component.from = "nevertreated", se = TRUE,
           nboots = 10, seed = 1)
  bg <- fit(method = "gsynth", se = TRUE, nboots = 10, seed = 1)
  expect_equal(b$est.avg, bg$est.avg, tolerance = 1e-12)
})


## -- C1  formula terms must be bare column names ----------------------------

.fix246c_simdata <- function() {
  if (!exists("simdata", inherits = TRUE)) utils::data("simdata", package = "fect")
  d <- get("simdata", inherits = TRUE)
  d$g3 <- d$id %% 3
  d
}

test_that("C1: fect() stops on formula terms that are not bare column names", {
  d <- .fix246c_simdata()
  fit <- function(f) fect::fect(f, data = d, index = c("id", "time"),
                                method = "fe", se = FALSE, parallel = FALSE)
  ## before 2.4.6 these ran silently: log(Y + 20) ~ D fitted Y ~ D,
  ## factor(g3) became one linear slope, X1 * X2 became X1 + X2
  expect_error(fit(log(Y + 20) ~ D), "not a column name: `log\\(Y \\+ 20\\)`")
  expect_error(fit(Y ~ D + factor(g3)), "not a column name: `factor\\(g3\\)`")
  expect_error(fit(Y ~ D + X1 * X2), "`X1 \\* X2`")
  expect_error(fit(Y ~ D + X1:X2), "`X1:X2`")
  expect_error(fit(Y ~ D + I(X1^2) + log(X2 + 10)),
               "not column names: `I\\(X1\\^2\\)`, `log\\(X2 \\+ 10\\)`")
  expect_error(fit(Y ~ D + X1 - X2), "`-X2`")
  expect_error(fit(Y ~ D + 2), "`2`")
  expect_error(fit(Y ~ .), "not a column name: `\\.`")
  expect_error(fit(Y ~ D + X1 + Y), "outcome \"Y\" also appears on the right-hand side")
  expect_error(fit(Y ~ 1), "needs a treatment variable")
  ## the message says what to do
  expect_error(fit(Y ~ D + factor(g3)), "model.matrix")
})

test_that("C1: bare names and intercept specifiers fit exactly as before", {
  skip_on_cran()
  d <- .fix246c_simdata()
  fit <- function(f) .fix246c_quiet(fect::fect(
    f, data = d, index = c("id", "time"), method = "fe", se = FALSE,
    parallel = FALSE
  ))
  base <- fit(Y ~ D + X1 + X2)
  forms <- list(Y ~ D + X1 + X2 + 0, Y ~ 0 + D + X1 + X2, Y ~ D + X1 + X2 - 1,
                Y ~ -1 + D + X1 + X2, Y ~ 1 + D + X1 + X2, Y ~ D + X1 + X2 + 1,
                Y ~ D + X1 + X2 - 0, Y ~ D + X1 + X2 + X1)
  for (f in forms) {
    g <- fit(f)
    expect_identical(g$att.avg, base$att.avg, info = deparse1(f))
    expect_identical(g$eff, base$eff, info = deparse1(f))
    expect_identical(g$beta, base$beta, info = deparse1(f))
    expect_identical(g$X, c("X1", "X2"), info = deparse1(f))
  }
  ## the formula and the column-name interfaces agree
  h <- .fix246c_quiet(fect::fect(
    data = d, Y = "Y", D = "D", X = c("X1", "X2"), index = c("id", "time"),
    method = "fe", se = FALSE, parallel = FALSE
  ))
  expect_identical(h$att.avg, base$att.avg)
})

test_that("C1: interFE() formulas take bare column names only", {
  d <- .fix246c_simdata()
  ## #62: D:time was read as D and time
  expect_error(fect::interFE(Y ~ X1 + D:time, data = d, index = c("id", "time")),
               "interFE\\(\\) formulas take bare column names only; not a column name: `D:time`")
  expect_error(fect::interFE(log(Y + 20) ~ X1, data = d, index = c("id", "time")),
               "`log\\(Y \\+ 20\\)`")
  expect_error(fect::interFE(Y ~ X1 + nosuchcol, data = d, index = c("id", "time")),
               "variable \"nosuchcol\" is not in the data set")
  ## bare names keep working and match the column-name interface (tol is
  ## given because the two methods have different defaults)
  a <- fect::interFE(Y ~ X1 + X2, data = d, index = c("id", "time"), r = 1,
                     tol = 1e-5)
  b <- fect::interFE(data = d, Y = "Y", X = c("X1", "X2"),
                     index = c("id", "time"), r = 1, tol = 1e-5)
  expect_identical(a$beta, b$beta)
  expect_identical(a$X, c("X1", "X2"))
})


## -- C2  X with a formula; duplicated covariate names ------------------------

test_that("C2: X given together with a formula stops (fect and interFE)", {
  d <- .fix246c_simdata()
  fit <- function(...) fect::fect(data = d, index = c("id", "time"),
                                  method = "fe", se = FALSE, parallel = FALSE, ...)
  ## before 2.4.6 X was ignored and the fit had no covariates
  expect_error(fit(Y ~ D, X = c("X1", "X2")),
               "Covariates were given in `X` together with a formula")
  expect_error(fit(Y ~ D + X1, X = "X2"), "not both")
  ## an unquoted name was never even evaluated
  expect_error(fit(Y ~ D, X = X1), "together with a formula")
  ## X = NULL with a formula is fine (gsynth's wrapper passes it)
  f <- .fix246c_quiet(fit(Y ~ D + X1, X = NULL))
  expect_identical(f$X, "X1")
  expect_error(fect::interFE(Y ~ X1, data = d, index = c("id", "time"), X = "X2"),
               "Covariates were given in `X` together with a formula")
})

test_that("C2: duplicated covariate names, or X naming Y or D, stop", {
  d <- .fix246c_simdata()
  fit <- function(X) fect::fect(data = d, Y = "Y", D = "D", X = X,
                                index = c("id", "time"), method = "fe",
                                se = FALSE, parallel = FALSE)
  ## before 2.4.6 a duplicated name entered the model twice (exactly
  ## collinear): garbage coefficients, or "inv(): matrix is singular" as here
  expect_error(fit(c("X1", "X2", "X1")), "duplicated covariate names: \"X1\"")
  expect_error(fit(c("X1", "D")), "outcome or the treatment variable \\(\"D\"\\)")
  expect_error(fit(c("Y", "X1")), "outcome or the treatment variable \\(\"Y\"\\)")
  expect_error(fit(c("X1", "nosuchcol")), "variable \"nosuchcol\" is not in the data set")
})


## -- C3  non-numeric covariates ---------------------------------------------

test_that("C3: non-numeric covariates stop with a clear message", {
  d <- .fix246c_simdata()
  d$Xf <- factor(ifelse(d$X1 > 0, "hi", "lo"))
  d$Xc <- ifelse(d$X1 > 0, "hi", "lo")
  d$Xd <- as.Date("2000-01-01") + seq_len(nrow(d))
  fit <- function(f) fect::fect(f, data = d, index = c("id", "time"),
                                method = "fe", se = FALSE, parallel = FALSE)
  ## before 2.4.6: "Calling var(x) on a factor x is defunct" (factor) and a
  ## false "unit-invariant" stop (character)
  expect_error(fit(Y ~ D + Xf),
               "Covariate \"Xf\" is a factor; fect\\(\\) needs numeric covariates")
  expect_error(fit(Y ~ D + Xf), "model.matrix\\(~ Xf, data\\)")
  expect_error(fit(Y ~ D + X1 + Xc), "Covariate \"Xc\" is character")
  expect_error(fit(Y ~ D + X1 + Xd), "Covariate \"Xd\" is of class \"Date\"")
})

test_that("C3: logical covariates are used as 0/1", {
  skip_on_cran()
  d <- .fix246c_simdata()
  d$Xl <- d$X1 > 0
  d$Xn <- as.numeric(d$Xl)
  fit <- function(f) .fix246c_quiet(fect::fect(
    f, data = d, index = c("id", "time"), method = "fe", se = FALSE,
    parallel = FALSE
  ))
  a <- fit(Y ~ D + X2 + Xl)
  b <- fit(Y ~ D + X2 + Xn)
  expect_identical(a$att.avg, b$att.avg)
  expect_identical(unname(a$beta), unname(b$beta))
})


## -- C4  time index given as a factor or as character ------------------------

test_that("C4: a factor time index is used in its level order", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  fit <- function(d, m) .fix246c_quiet(fect::fect(
    Y ~ D + X1 + X2, data = d, index = c("id", "time"), method = m, r = 2,
    CV = FALSE, se = FALSE, parallel = FALSE
  ))
  base.g <- fit(simgsynth, "gsynth")
  base.i <- fit(simgsynth, "ife")
  ## levels that are increasing numbers: the same fit as the numbers.
  ## Before 2.4.6 the levels were sorted as text ("1", "10", "11", ...):
  ## gsynth stopped with a false "reversals" error and ife ran with 22
  ## event times instead of 30.
  d1 <- simgsynth
  d1$time <- factor(d1$time)
  g1 <- fit(d1, "gsynth")
  i1 <- fit(d1, "ife")
  expect_identical(g1$att.avg, base.g$att.avg)
  expect_identical(g1$eff, base.g$eff)
  expect_identical(g1$rawtime, base.g$rawtime)
  expect_identical(i1$att.avg, base.i$att.avg)
  expect_identical(i1$time, base.i$time)
  ## other levels: level order, with the labels kept for the output
  lev <- paste0("p", 1:30)
  d2 <- simgsynth
  d2$time <- factor(paste0("p", d2$time), levels = lev)
  g2 <- fit(d2, "gsynth")
  expect_identical(g2$att.avg, base.g$att.avg)
  expect_identical(unname(g2$eff), unname(base.g$eff))
  expect_identical(g2$rawtime, lev)
  expect_identical(rownames(g2$eff), lev)
  expect_true(is.factor(g2$data.long$time))
  expect_identical(levels(g2$data.long$time), lev)
  expect_setequal(as.character(g2$data.long$time), lev)
})

test_that("C4: a period without controls is reported by its label", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  lev <- paste0("p", 1:30)
  d <- simgsynth
  d$time <- factor(paste0("p", d$time), levels = lev)
  tr <- unique(d$id[d$D == 1])
  ## no control unit is observed in the last period
  d <- d[!(d$time == "p30" & !d$id %in% tr), ]
  res <- .fix246c_messages(fect::fect(
    Y ~ D + X1 + X2, data = d, index = c("id", "time"), method = "gsynth",
    r = 2, CV = FALSE, se = FALSE, parallel = FALSE
  ))
  expect_true(any(grepl("under control at p30, drop that period", res$messages)))
  expect_identical(res$value$rawtime, lev[1:29])
})

test_that("C4: a character time index must hold numbers", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  fit <- function(d) .fix246c_quiet(fect::fect(
    Y ~ D + X1 + X2, data = d, index = c("id", "time"), method = "gsynth",
    r = 2, CV = FALSE, se = FALSE, parallel = FALSE
  ))
  base <- fit(simgsynth)
  d3 <- simgsynth
  d3$time <- as.character(d3$time)
  g3 <- fit(d3)
  ## before 2.4.6: ordered as text, and gsynth stopped ("reversals")
  expect_identical(g3$att.avg, base$att.avg)
  expect_identical(g3$eff, base$eff)
  expect_identical(g3$rawtime, base$rawtime)
  d4 <- simgsynth
  d4$time <- sprintf("t%02d", d4$time)
  expect_error(fit(d4), paste0("The time index \"time\" is character and some ",
                               "values are not numbers \\(for example \"t01\"\\)"))
  ## Date indices are used as before
  d5 <- simgsynth
  d5$time <- as.Date("2000-01-01") + 31 * d5$time
  g5 <- fit(d5)
  expect_identical(g5$att.avg, base$att.avg)
  expect_s3_class(g5$rawtime, "Date")
})

test_that("C4b: numbers as factor levels out of numeric order give a warning; level order is kept", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  fit <- function(d) .fix246c_warnings(fect::fect(
    Y ~ D + X1 + X2, data = d, index = c("id", "time"), method = "ife",
    r = 2, CV = FALSE, se = FALSE, parallel = FALSE
  ))
  c4b <- "is a factor whose levels are numbers but not in increasing order"
  ## factor(as.character(time)) has the levels "1", "10", "11", ..., "9"
  d1 <- simgsynth
  d1$time <- factor(as.character(d1$time))
  w1 <- fit(d1)
  msg <- paste0(
    "The time index \"time\" is a factor whose levels are numbers but not in ",
    "increasing order (for example \"1\", \"10\", \"11\"); fect orders the ",
    "periods by the factor levels. If that is not the time order, pass a ",
    "numeric index (e.g. as.numeric(as.character(time))) or reorder the levels."
  )
  expect_identical(sum(w1$warnings == msg), 1L)
  expect_identical(sum(grepl(c4b, w1$warnings, fixed = TRUE)), 1L)
  ## the levels are still used in their order: the same fit as labels that
  ## are not numbers, in the same order
  d2 <- simgsynth
  d2$time <- factor(paste0("p", d2$time), levels = paste0("p", levels(d1$time)))
  w2 <- fit(d2)
  expect_false(any(grepl(c4b, w2$warnings, fixed = TRUE)))
  expect_identical(w1$value$att.avg, w2$value$att.avg)
  expect_identical(unname(w1$value$eff), unname(w2$value$eff))
  expect_identical(w1$value$rawtime, levels(d1$time))
  ## levels that are increasing numbers: no warning
  d3 <- simgsynth
  d3$time <- factor(d3$time)
  expect_false(any(grepl(c4b, fit(d3)$warnings, fixed = TRUE)))
})


## -- C5 + C6  collinear and FE-absorbed covariates --------------------------

## simgsynth plus covariates that cannot be estimated: X3 = 2 X1 and
## X4 = X1 + X2 (collinear), U (unit-level), Tm (period-level), UT = U + Tm,
## Z0 (zero for every never-treated unit).
.fix246c_cov_data <- function() {
  data(simgsynth, package = "fect")
  d <- simgsynth
  d$X3 <- 2 * d$X1
  d$X4 <- d$X1 + d$X2
  ids <- sort(unique(d$id))
  set.seed(11)
  d$U <- stats::rnorm(length(ids))[match(d$id, ids)]
  d$Tm <- stats::rnorm(max(d$time))[d$time]
  d$UT <- d$U + d$Tm
  tr <- unique(d$id[d$D == 1])
  d$Z0 <- ifelse(d$id %in% tr, stats::rnorm(nrow(d)), 0)
  d
}

test_that("C5: an exactly collinear covariate is dropped with a warning; estimates equal those without it", {
  skip_on_cran()
  d <- .fix246c_cov_data()
  for (m in c("gsynth", "ife")) {
    fit <- function(f) fect::fect(f, data = d, index = c("id", "time"),
                                  method = m, r = 2, CV = FALSE, se = FALSE,
                                  parallel = FALSE)
    base <- .fix246c_quiet(fit(Y ~ D + X1 + X2))
    ## before 2.4.6: gsynth att.avg 6.5098 instead of 5.5433 with no warning;
    ## ife stopped with "inv(): matrix is singular"
    for (f in list(Y ~ D + X1 + X2 + X3, Y ~ D + X1 + X2 + X4)) {
      w <- .fix246c_warnings(fit(f))
      xn <- all.vars(f)[5]
      expect_length(w$warnings, 1L)
      expect_match(w$warnings, paste0(
        "Dropped 1 covariate that cannot be estimated on the cells used to fit ",
        "the model: \"", xn, "\" is a linear combination of other covariates ",
        "\\(after removing the fixed effects\\)\\. Its coefficient is reported as NA\\."))
      g <- w$value
      expect_identical(g$att.avg, base$att.avg, info = paste(m, xn))
      expect_identical(g$eff, base$eff, info = paste(m, xn))
      expect_identical(rownames(g$beta), c("X1", "X2", xn))
      expect_identical(unname(g$beta[1:2, 1]), unname(base$beta[, 1]))
      expect_true(is.na(g$beta[3, 1]))
      expect_identical(g$X, c("X1", "X2", xn))
      expect_equal(g$validX, 1)
    }
  }
  ## lm()'s rule: the earlier column of a collinear set is kept
  w <- .fix246c_warnings(fect::fect(Y ~ D + X3 + X1 + X2, data = d,
                                    index = c("id", "time"), method = "fe",
                                    se = FALSE, parallel = FALSE))
  expect_match(w$warnings, "\"X1\" is a linear combination", all = FALSE)
  expect_true(is.na(w$value$beta["X1", 1]))
  expect_false(is.na(w$value$beta["X3", 1]))
})

test_that("C5: a covariate with no variation on the estimation cells is dropped (never-treated fit)", {
  skip_on_cran()
  d <- .fix246c_cov_data()
  fit <- function(f) fect::fect(f, data = d, index = c("id", "time"),
                                method = "gsynth", r = 2, CV = FALSE,
                                se = FALSE, parallel = FALSE)
  base <- .fix246c_quiet(fit(Y ~ D + X1 + X2))
  ## Z0 is 0 for every never-treated unit, which is where gsynth fits the
  ## covariates. Before 2.4.6 it was dropped silently.
  w <- .fix246c_warnings(fit(Y ~ D + X1 + X2 + Z0))
  expect_match(w$warnings, paste0("\"Z0\" has no variation on the cells used to ",
                                  "estimate the covariate coefficients"))
  expect_identical(w$value$att.avg, base$att.avg)
  expect_true(is.na(w$value$beta["Z0", 1]))
})

test_that("C5: bootstrap outputs keep one row per requested covariate", {
  skip_on_cran()
  d <- .fix246c_cov_data()
  fit <- function(f) fect::fect(f, data = d, index = c("id", "time"),
                                method = "gsynth", r = 2, CV = FALSE,
                                se = TRUE, nboots = 10, seed = 1,
                                parallel = FALSE)
  base <- .fix246c_quiet(fit(Y ~ D + X1 + X2))
  g <- .fix246c_quiet(fit(Y ~ D + X1 + X3 + X2))
  expect_identical(g$est.avg, base$est.avg)
  expect_identical(rownames(g$est.beta), c("X1", "X3", "X2"))
  expect_identical(colnames(g$est.beta), colnames(base$est.beta))
  expect_identical(unname(g$est.beta[c(1, 3), ]), unname(base$est.beta))
  expect_true(all(is.na(g$est.beta["X3", ])))
  expect_identical(dim(g$beta.boot), c(3L, ncol(base$beta.boot)))
  expect_identical(unname(g$beta.boot[c(1, 3), ]), unname(base$beta.boot))
  expect_true(all(is.na(g$beta.boot[2, ])))
})

test_that("C5: the C++ inverse of X'X falls back to a generalized inverse when singular", {
  set.seed(3)
  x <- array(stats::rnorm(24), dim = c(4, 3, 2))
  w1 <- matrix(1, 4, 3)
  xx <- matrix(c(sum(x[, , 1]^2), sum(x[, , 1] * x[, , 2]),
                 sum(x[, , 1] * x[, , 2]), sum(x[, , 2]^2)), 2, 2)
  ## well-conditioned: the ordinary inverse, as before
  expect_equal(fect:::XXinv(x), solve(xx), tolerance = 1e-10)
  expect_equal(fect:::wXXinv(x, w1), solve(xx), tolerance = 1e-10)
  ## singular: before 2.4.6 "inv(): matrix is singular". The fallback is a
  ## symmetric generalized inverse G of X'X (X'X G X'X = X'X, G X'X G = G),
  ## also when the collinear covariates are on very different scales.
  for (k in c(2, 2e8)) {
    x[, , 2] <- k * x[, , 1]
    xs <- matrix(c(1, k, k, k^2) * sum(x[, , 1]^2), 2, 2)
    for (G in list(fect:::XXinv(x), fect:::wXXinv(x, w1))) {
      expect_true(all(is.finite(G)))
      expect_equal(G, t(G), tolerance = 1e-10)
      expect_equal(xs %*% G %*% xs, xs, tolerance = 1e-10)
      expect_equal(G %*% xs %*% G, G, tolerance = 1e-10)
    }
  }
  ## an all-zero covariate: its row and column of G are 0
  x[, , 2] <- 0
  for (G in list(fect:::XXinv(x), fect:::wXXinv(x, w1))) {
    expect_equal(G, matrix(c(1 / sum(x[, , 1]^2), 0, 0, 0), 2, 2),
                 tolerance = 1e-12)
  }
})

test_that("C5b: the C++ inverse of X'X does not depend on the covariates' units", {
  set.seed(3)
  x <- array(stats::rnorm(60), dim = c(5, 4, 3))
  ## full rank, but columns in very different units: rcond of the raw X'X is
  ## about 1e-22, rcond after scaling X'X to unit diagonal is not small
  x[, , 1] <- x[, , 1] * 1e8
  x[, , 3] <- x[, , 3] * 1e-3
  w <- matrix(stats::runif(20, 0.5, 2), 5, 4)
  gram <- function(x, w) {
    p <- dim(x)[3]
    outer(seq_len(p), seq_len(p),
          Vectorize(function(k, m) sum(w * x[, , k] * x[, , m])))
  }
  for (wt in list(NULL, w)) {
    xx <- gram(x, if (is.null(wt)) 1 else wt)
    G <- if (is.null(wt)) fect:::XXinv(x) else fect:::wXXinv(x, wt)
    ## the ordinary inverse; compared on the unit-diagonal scale so that every
    ## entry counts. The first version of the guard returned pinv(X'X) here,
    ## which zeroes the directions of the small-scale covariates.
    sc <- sqrt(diag(xx))
    expect_equal(G * outer(sc, sc), solve(xx / outer(sc, sc)), tolerance = 1e-8)
  }
})

test_that("C5b: covariates in very different units give the fit in the original units", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  d0 <- simgsynth
  set.seed(4)
  d0$wt <- stats::runif(nrow(d0), 0.5, 2)
  fit <- function(d, f, args) .fix246c_quiet(do.call(fect::fect, c(list(
    f, data = d, index = c("id", "time"), CV = FALSE, se = FALSE,
    parallel = FALSE), args)))
  specs <- list(
    fe = list(method = "fe"),
    ife = list(method = "ife", r = 2),
    mc = list(method = "mc", lambda = 0.01),
    gsynth = list(method = "gsynth", r = 2),
    ife.weighted = list(method = "ife", r = 2, W = "wt")
  )
  for (nm in names(specs)) {
    base <- fit(d0, Y ~ D + X1 + X2, specs[[nm]])
    for (s in c(1e8, 1e9)) {
      d <- d0
      d$GDP <- d0$X1 * s
      g <- fit(d, Y ~ D + GDP + X2, specs[[nm]])
      ## first version of the C5 guard: fe att.avg 4.0686 instead of 5.0852,
      ## gsynth 4.9664 instead of 5.5433, X2 coefficient about 1e-16
      info <- paste(nm, "x", s)
      expect_equal(g$att.avg, base$att.avg, tolerance = 1e-10, info = info)
      expect_equal(g$eff, base$eff, tolerance = 1e-8, info = info)
      expect_equal(unname(g$beta[, 1] * c(s, 1)), unname(base$beta[, 1]),
                   tolerance = 1e-8, info = info)
    }
  }
})

test_that("C5b: parametric SEs do not depend on the covariates' units", {
  skip_on_cran()
  data(simgsynth, package = "fect")
  fit <- function(d, f) .fix246c_quiet(fect::fect(
    f, data = d, index = c("id", "time"), method = "gsynth", r = 2,
    CV = FALSE, se = TRUE, vartype = "parametric", nboots = 20, seed = 1,
    parallel = FALSE
  ))
  base <- fit(simgsynth, Y ~ D + X1 + X2)
  d <- simgsynth
  d$GDP <- d$X1 * 1e8
  g <- fit(d, Y ~ D + GDP + X2)
  ## first version of the C5 guard: S.E. 0.9256 instead of 0.2886
  expect_equal(g$att.avg, base$att.avg, tolerance = 1e-10)
  expect_equal(g$est.avg[, "S.E."], base$est.avg[, "S.E."], tolerance = 1e-4)
})

test_that("C6: covariates absorbed by the fixed effects are dropped with a warning naming the fixed effect", {
  skip_on_cran()
  d <- .fix246c_cov_data()
  fit <- function(f, m = "gsynth") fect::fect(
    f, data = d, index = c("id", "time"), method = m, r = 2, CV = FALSE,
    se = FALSE, parallel = FALSE, force = "two-way"
  )
  base <- .fix246c_quiet(fit(Y ~ D + X1 + X2))
  ## before 2.4.6: stop with swapped labels (U was called "unit-invariant")
  reasons <- c(
    U = "does not vary over time within units, so it is absorbed by the unit fixed effects",
    Tm = "does not vary across units within periods, so it is absorbed by the time fixed effects",
    UT = "is the sum of a unit-level and a period-level variable, so it is absorbed by the unit and time fixed effects"
  )
  for (v in names(reasons)) {
    w <- .fix246c_warnings(fit(stats::reformulate(c("D", "X1", "X2", v), "Y")))
    expect_match(w$warnings, paste0("\"", v, "\" ", reasons[[v]]), fixed = TRUE)
    expect_identical(w$value$att.avg, base$att.avg, info = v)
    expect_identical(w$value$eff, base$eff, info = v)
    expect_true(is.na(w$value$beta[v, 1]))
  }
})

test_that("C6: when every covariate is dropped the fit has no covariates", {
  skip_on_cran()
  d <- .fix246c_cov_data()
  fit <- function(f) fect::fect(f, data = d, index = c("id", "time"),
                                method = "fe", se = FALSE, parallel = FALSE)
  base <- .fix246c_quiet(fit(Y ~ D))
  w <- .fix246c_warnings(fit(Y ~ D + U + Tm))
  ## one warning (no extra "Multi-colinearity" warning)
  expect_length(w$warnings, 1L)
  expect_match(w$warnings, "^Dropped 2 covariates .*Their coefficients are reported as NA\\.$")
  expect_identical(w$value$att.avg, base$att.avg)
  expect_identical(w$value$eff, base$eff)
  expect_identical(dim(w$value$beta), c(2L, 1L))
  expect_identical(rownames(w$value$beta), c("U", "Tm"))
  expect_true(all(is.na(w$value$beta)))
  expect_equal(w$value$validX, 0)
})

test_that("C6: a covariate absorbed by an extra cfe fixed effect is dropped", {
  skip_on_cran()
  d <- .fix246c_cov_data()
  d$region <- d$id %% 3
  d$rt <- d$region * 100 + d$time
  set.seed(5)
  d$RT <- stats::rnorm(1000)[d$rt] # constant within region x period
  fit <- function(f) fect::fect(f, data = d, index = c("id", "time", "rt"),
                                method = "cfe", r = 0, se = FALSE,
                                parallel = FALSE)
  base <- .fix246c_quiet(fit(Y ~ D + X1))
  w <- .fix246c_warnings(fit(Y ~ D + X1 + RT))
  expect_match(w$warnings, "\"RT\" is absorbed by the fixed effects", fixed = TRUE)
  expect_identical(w$value$att.avg, base$att.avg)
})

test_that("C6: covariates are checked only against the fixed effects in the model", {
  skip_on_cran()
  d <- .fix246c_cov_data()
  ## with unit effects only, a period-level covariate is estimable
  w <- .fix246c_warnings(fect::fect(Y ~ D + X1 + Tm, data = d,
                                    index = c("id", "time"), method = "fe",
                                    force = "unit", se = FALSE,
                                    parallel = FALSE))
  expect_length(w$warnings, 0L)
  expect_true(all(is.finite(w$value$beta)))
})

test_that("C6: interFE() names the fixed effect that absorbs a covariate", {
  d <- .fix246c_cov_data()
  ## before 2.4.6 the labels were swapped ("unit-invariant" for U)
  expect_error(fect::interFE(Y ~ X1 + U, data = d, index = c("id", "time")),
               "Variable \"U\" does not vary over time within units \\(it is absorbed by the unit fixed effects\\)")
  expect_error(fect::interFE(Y ~ X1 + Tm, data = d, index = c("id", "time")),
               "Variable \"Tm\" does not vary across units within periods \\(it is absorbed by the time fixed effects\\)")
  ## a covariate is checked only against the fixed effects in the model
  f <- fect::interFE(Y ~ X1 + Tm, data = d, index = c("id", "time"), force = "unit")
  expect_true(all(is.finite(f$beta)))
})

test_that("C5/C6: dropped covariates get NA rows in every beta output", {
  ## the restore helper on a hand-built fit (a binary model with SEs)
  out <- list(beta = matrix(c(1, 2), 2, 1), marginal = matrix(c(0.1, 0.2), 2, 1),
              est.beta = matrix(1:10, 2, 5, dimnames = list(NULL, letters[1:5])),
              est.marginal = matrix(11:20, 2, 5), beta.boot = matrix(1:6, 2, 3),
              att.avg.boot = matrix(0, 1, 3), validX = 1)
  r <- fect:::.fect_restore_dropped_covariates(out, keep = c(1L, 3L), p.all = 3L,
                                               se = TRUE, binary = TRUE)
  expect_identical(r$beta[, 1], c(1, NA, 2))
  expect_identical(r$marginal[, 1], c(0.1, NA, 0.2))
  expect_identical(r$est.beta[2, ], setNames(rep(NA_real_, 5), letters[1:5]))
  expect_identical(unname(r$est.beta[3, 1]), 2)
  expect_identical(dim(r$est.marginal), c(3L, 5L))
  expect_identical(r$beta.boot[, 1], c(1, NA, 2))
  expect_equal(r$validX, 1)
  ## nothing kept: NA matrices of the right shape, validX 0
  r0 <- fect:::.fect_restore_dropped_covariates(
    list(beta = NA, att.avg.boot = matrix(0, 1, 4), validX = 0),
    keep = integer(0), p.all = 2L, se = TRUE, binary = TRUE)
  expect_identical(dim(r0$beta), c(2L, 1L))
  expect_identical(dim(r0$marginal), c(2L, 1L))
  expect_identical(dim(r0$est.beta), c(2L, 5L))
  expect_identical(dim(r0$est.marginal), c(2L, 5L))
  expect_identical(dim(r0$beta.boot), c(2L, 4L))
  expect_equal(r0$validX, 0)
})
