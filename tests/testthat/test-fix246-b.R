## ---------------------------------------------------------------
## fect 2.4.6 correctness fixes, group B (CV, weights in CV, implied
## weights, plots, guards). Run 2026-09-24-fix246-correctness; one block
## (or more) per item (B1-B9). Every test here fails on dev @ 412d7ae
## and passes after the fix.
## ---------------------------------------------------------------

.fix246b_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## Collect the messages of `expr` (warnings are muffled); errors propagate.
.fix246b_messages <- function(expr) {
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

## Balanced panel with `Nco` never-treated units and `Ntr` units treated
## from period T0 + 1 on; two factors, two covariates.
.fix246b_panel <- function(Nco, Ntr = 3, TT = 30, T0 = 20, seed = 32) {
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


## -- B1  criterion = "pc" -------------------------------------------------

test_that("B1: criterion = 'pc' picks the lowest-PC r for gsynth and never-treated CFE", {
  skip_on_cran()
  data(turnout, package = "fect")
  ## On turnout the PC is lowest at r = 4 while the MSPE rule picks r = 2.
  for (m in c("gsynth", "cfe")) {
    fit <- .fix246b_quiet(fect::fect(
      turnout ~ policy_edr + policy_mail_in + policy_motor, data = turnout,
      index = c("abb", "year"), method = m, time.component.from = "nevertreated",
      force = "two-way", CV = TRUE, r = c(0, 5), k = 5, criterion = "pc",
      se = FALSE, parallel = FALSE, seed = 1
    ))
    r.pc <- unname(fit$CV.out[which.min(fit$CV.out[, "PC"]), "r"])
    expect_equal(r.pc, 4, info = m)
    expect_equal(unname(fit$r.cv), r.pc, info = m)
    ## the final fit is estimated at the selected r
    expect_equal(ncol(fit$factor), r.pc, info = m)
  }
})

test_that("B1: fect's IFE CV table under criterion = 'pc' has the right columns", {
  skip_on_cran()
  data(turnout, package = "fect")
  fit <- .fix246b_quiet(fect::fect(
    turnout ~ policy_edr + policy_mail_in + policy_motor, data = turnout,
    index = c("abb", "year"), method = "ife", force = "two-way", CV = TRUE,
    r = c(0, 5), k = 5, criterion = "pc", se = FALSE, parallel = FALSE, seed = 1
  ))
  tab <- fit$CV.out.ife
  expect_identical(colnames(tab), c("r", "sigma2", "IC", "PC", "MSPTATT", "MSE"))
  ## every row holds its own statistics (no sentinel or NA left behind)
  expect_true(all(is.finite(tab)))
  expect_true(all(tab[, "PC"] < 1e19))
  expect_equal(unname(fit$r.cv), unname(tab[which.min(tab[, "PC"]), "r"]))
})
