## ---------------------------------------------------------------
## Regression test: r.cv.rolling supports both method = "ife"
## (IFE-EM, time.component.from = "notyettreated") and
## method = "gsynth" (GSC, time.component.from = "nevertreated").
##
## Pre-fold-in (PR #119 as merged): r.cv.rolling guarded on
## method == "ife" and hard-coded the inner fect call to
## time.component.from = "notyettreated", because the GSC path did not
## populate Y.ct.full at masked control positions.
##
## Stage 2 (committed separately) fixed Y.ct.full on the GSC path.
## Stage 3 fold-in (this branch):
##   - Drops the method == "ife" guard in r.cv.rolling.
##   - Accepts method = "gsynth" and routes the inner fect call to
##     time.component.from = "nevertreated".
##
## The bundled fect demo datasets (simdata, simgsynth) are balanced
## panels: masking the last cv.nobs observations of every control
## unit drops those time points entirely, so cv.sample.rolling cannot
## exercise method = "gsynth" cleanly on them. We synthesize a small
## unbalanced panel inline for the gsynth test, and use simdata for
## the IFE-EM test (works because simdata has many ineligible /
## treated units that retain data at the masked time points).
## ---------------------------------------------------------------

.make_unbalanced_panel <- function(seed = 11) {
  set.seed(seed)
  N  <- 30L; T <- 25L
  ids   <- rep(seq_len(N), each = T)
  times <- rep(seq_len(T), N)
  F1 <- rnorm(T); L1 <- rnorm(N)
  err <- rnorm(N * T, sd = 0.3)
  mu  <- rep(L1, each = T) * rep(F1, N)
  Y   <- mu + err
  treated      <- 26L:30L
  treat_starts <- c(18L, 19L, 20L, 21L, 22L)
  D <- rep(0L, N * T)
  for (k in seq_along(treated)) {
    u <- treated[k]
    D[ids == u & times >= treat_starts[k]] <- 1L
  }
  Y <- Y + D * 1.5
  df <- data.frame(id = ids, time = times, Y = Y, D = D)
  ## Drop trailing observations from a few units to make panel unbalanced
  ## so the rolling-CV mask doesn't drop entire time points.
  drop_keys <- c(
    paste(27L, c(24L, 25L), sep = "_"),
    paste(28L, c(22L, 23L, 24L, 25L), sep = "_"),
    paste(1L,  c(23L, 24L, 25L), sep = "_"),
    paste(2L,  c(24L, 25L), sep = "_")
  )
  df_keys <- paste(df$id, df$time, sep = "_")
  df[!df_keys %in% drop_keys, ]
}

test_that("r.cv.rolling works for method = 'ife' (IFE-EM)", {

  skip_on_cran()

  if (!exists("simdata", inherits = TRUE)) {
    suppressWarnings(try(utils::data("simdata", package = "fect"),
                         silent = TRUE))
  }
  skip_if_not(exists("simdata", inherits = TRUE),
              "simdata dataset not available")

  set.seed(42)
  res <- suppressWarnings(suppressMessages(
    fect::r.cv.rolling(
      Y ~ D, data = simdata, index = c("id", "time"),
      method  = "ife",
      r.max   = 3,
      cv.nobs = 3,
      cv.rule = "1se",
      verbose = FALSE
    )
  ))

  expect_true(is.list(res))
  expect_true(is.integer(res$r.cv))
  expect_true(res$r.cv >= 0L && res$r.cv <= 3L)
  expect_true(is.data.frame(res$mspe))
  expect_equal(nrow(res$mspe), 4L)
  expect_true(any(is.finite(res$mspe$mspe)))
})

test_that("r.cv.rolling works for method = 'gsynth' (GSC)", {

  skip_on_cran()

  df <- .make_unbalanced_panel(seed = 11)

  set.seed(42)
  res <- suppressWarnings(suppressMessages(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method  = "gsynth",
      r.max   = 2,
      cv.nobs = 2,
      cv.rule = "1se",
      force   = "two-way",
      verbose = FALSE
    )
  ))

  expect_true(is.list(res))
  expect_true(is.integer(res$r.cv))
  expect_true(res$r.cv >= 0L && res$r.cv <= 2L)
  expect_true(is.data.frame(res$mspe))
  expect_equal(nrow(res$mspe), 3L)
  ## At least one r must have produced a finite MSPE (otherwise the GSC
  ## inner fits silently fell through and the dispatch is broken).
  expect_true(any(is.finite(res$mspe$mspe)))
})

test_that("r.cv.rolling rejects unknown method", {

  skip_on_cran()

  df <- .make_unbalanced_panel(seed = 11)

  expect_error(
    fect::r.cv.rolling(
      Y ~ D, data = df, index = c("id", "time"),
      method  = "mc",
      r.max   = 1,
      verbose = FALSE
    ),
    regexp = "should be one of"
  )
})
