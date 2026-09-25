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
