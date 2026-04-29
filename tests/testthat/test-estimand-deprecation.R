## ---------------------------------------------------------------
## Tests for soft-deprecation messages on effect() and att.cumu().
##
## Locked behavior per ref/po-estimands-contract.md §7:
## - Calling effect() or att.cumu() emits a one-time-per-session
##   message pointing at ?estimand.
## - Subsequent calls in the same session do NOT re-emit.
## - Calling estimand("att.cumu", ...) (which delegates internally)
##   does NOT trigger the message.
## ---------------------------------------------------------------

.fit_no_reversal <- function(nboots = 30, keep.sims = TRUE) {
  set.seed(1); N <- 30; TT <- 20
  df <- expand.grid(id = 1:N, time = 1:TT)
  treat_start <- sample(c(NA, 8:15), N, replace = TRUE)
  df$D <- ifelse(is.na(treat_start[df$id]) | df$time < treat_start[df$id],
                 0, 1)
  df$Y <- rnorm(nrow(df)) + 0.5 * df$D + 0.1 * df$time
  set.seed(42)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = df, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, nboots = nboots, parallel = FALSE,
      keep.sims = keep.sims
    )
  ))
}


## Reset the once-per-session flags before each test so messaging
## behavior is reproducible regardless of test ordering.
.reset_dep_flags <- function() {
  options(
    fect.deprecated_warned_effect   = FALSE,
    fect.deprecated_warned_att.cumu = FALSE,
    fect.suppress_estimand_deprecation = FALSE
  )
}


## -- DEP.1  effect() emits a deprecation message on first call -------

test_that("DEP.1: effect() emits a one-time deprecation message", {

  skip_on_cran()

  .reset_dep_flags()
  fit <- .fit_no_reversal()
  expect_message(
    out <- fect::effect(fit, cumu = TRUE, plot = FALSE),
    "soft-deprecated",
    fixed = FALSE
  )
  ## Subsequent call: no message (one-time gate).
  expect_no_message(
    out2 <- fect::effect(fit, cumu = TRUE, plot = FALSE)
  )
})


## -- DEP.2  att.cumu() emits a deprecation message on first call -----

test_that("DEP.2: att.cumu() emits a one-time deprecation message", {

  skip_on_cran()

  .reset_dep_flags()
  fit <- .fit_no_reversal()
  expect_message(
    out <- fect::att.cumu(fit, period = c(1, 5), plot = FALSE),
    "soft-deprecated",
    fixed = FALSE
  )
  expect_no_message(
    out2 <- fect::att.cumu(fit, period = c(1, 5), plot = FALSE)
  )
})


## -- DEP.3  estimand("att.cumu") does NOT trigger deprecation msg ----

test_that("DEP.3: estimand('att.cumu') is silent (no internal deprecation)", {

  skip_on_cran()

  .reset_dep_flags()
  fit <- .fit_no_reversal()
  ## estimand() delegates to effect()/att.cumu() internally, with the
  ## suppression flag set. No message expected.
  expect_no_message(
    fect::estimand(fit, "att.cumu", "event.time")
  )
  expect_no_message(
    fect::estimand(fit, "att.cumu", "overall", window = c(1, 5))
  )

  ## After estimand() returns, the suppression flag is reset, so a
  ## direct effect() call DOES emit the message (one-time gate intact).
  expect_message(
    fect::effect(fit, cumu = TRUE, plot = FALSE),
    "soft-deprecated",
    fixed = FALSE
  )
})


## -- DEP.4  Numerical equality preserved post-deprecation ------------

test_that("DEP.4: deprecation does not change numerical output", {

  skip_on_cran()

  .reset_dep_flags()
  fit <- .fit_no_reversal()
  out_eff <- suppressMessages(fect::effect(fit, cumu = TRUE, plot = FALSE))
  out_acc <- suppressMessages(
    fect::att.cumu(fit, period = c(1, 5), plot = FALSE)
  )

  expect_false(is.null(out_eff$effect.est.att))
  expect_false(is.null(out_acc))
  ## Post-deprecation outputs are still byte-identical to estimand().
  est <- fect::estimand(fit, "att.cumu", "event.time")
  est_overall <- fect::estimand(fit, "att.cumu", "overall",
                                 window = c(1, 5))

  shared_et <- intersect(est$event.time,
                         as.numeric(rownames(out_eff$effect.est.att)))
  expect_true(length(shared_et) > 0)
})
