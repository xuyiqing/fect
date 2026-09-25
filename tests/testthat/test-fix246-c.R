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
                 pl$att[which(pl$time == kk)], tolerance = 1e-10, info = kk)
  }
})
