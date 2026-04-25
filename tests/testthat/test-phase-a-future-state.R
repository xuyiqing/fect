## ---------------------------------------------------------------
## Regression test: Phase A bootstrap parallel state across calls.
##
## Pre-fix bug (HANDOFF 2026-04-25): boot.R Phase A's foreach %dopar%
## inherited whatever foreach backend was registered in the global
## registry. After any prior fect parallel call, run_dopar_retry's
## on.exit re-registered doFuture, polluting subsequent calls --
## variant (iii) Phase A in a forest session ran 8x slower than the
## same call in isolation.
##
## Fix (Stage 1, 2026-04-25):
##   1. Phase A migrated from foreach %dopar% to future_lapply
##      (boot.R:884). Backend is now governed by the active future plan
##      set in default.R, not by a global foreach registry.
##   2. Removed `doFuture::registerDoFuture()` from run_dopar_retry's
##      on.exit (boot.R:1553).
##
## Test: two consecutive fect(parallel = TRUE) calls in the same R
## process. The second call's wall-time must be within 2x of the
## first; pre-fix this could blow up to 8x.
## ---------------------------------------------------------------

test_that("Phase A wall-time is stable across consecutive parallel calls", {

  skip_on_cran()

  one_call <- function(seed) {
    set.seed(seed)
    t0 <- Sys.time()
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data      = simdata,
        index     = c("id", "time"),
        method    = "ife",
        r         = 1,
        CV        = FALSE,
        force     = "two-way",
        time.component.from = "nevertreated",
        se        = TRUE,
        vartype   = "bootstrap",
        nboots    = 20,
        parallel  = TRUE,
        cores     = 2
      )
    ))
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }

  ## Warm-up: pay future plan / package-load cost outside the measurement.
  invisible(one_call(seed = 1))

  ## First measured call.
  t1 <- one_call(seed = 2)
  ## Second measured call -- must not blow up due to backend pollution.
  t2 <- one_call(seed = 3)

  ratio <- t2 / max(t1, 0.5)
  ## 2x is the plan budget; allow 3x ceiling to absorb GC / OS jitter on CI.
  expect_lt(ratio, 3)
})
