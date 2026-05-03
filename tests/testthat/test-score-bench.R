## ---------------------------------------------------------------
## Runtime benchmarks for cv.method timing comparison.
## Heavy: each test fits the same DGP under {loo, all_units, treated_units}
## and prints elapsed time per call. Skip on CRAN.
##
## Originally Section I of test-score-unify.R; split out 2026-05-03.
## Shared fixtures live in helper-score-unify.R.
## ---------------------------------------------------------------

## =================================================================
## Section I: Runtime Benchmarks (cv.method timing comparison)
## Informational only — skip on CRAN
## =================================================================

test_that("BENCH1: IFE timing comparison (loo vs all_units vs treated_units)", {

  skip_on_cran()

  cat("\n=== Runtime Benchmark: IFE nevertreated cv.method timing ===\n")

  set.seed(42)
  t_loo <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        time.component.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = "loo",
        se              = FALSE,
        parallel        = FALSE
      )
    ))
  )
  cat(sprintf("  IFE loo:            %6.2f sec (elapsed)\n", t_loo["elapsed"]))

  set.seed(42)
  t_au <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        time.component.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = "all_units",
        se              = FALSE,
        parallel        = FALSE
      )
    ))
  )
  cat(sprintf("  IFE all_units:      %6.2f sec (elapsed)\n", t_au["elapsed"]))

  set.seed(42)
  t_tu <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        time.component.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = "treated_units",
        se              = FALSE,
        parallel        = FALSE
      )
    ))
  )
  cat(sprintf("  IFE treated_units:  %6.2f sec (elapsed)\n", t_tu["elapsed"]))

  cat(sprintf("  Speedup (loo/all_units):      %.2fx\n",
              t_loo["elapsed"] / max(t_au["elapsed"], 0.001)))
  cat(sprintf("  Speedup (loo/treated_units):  %.2fx\n",
              t_loo["elapsed"] / max(t_tu["elapsed"], 0.001)))
  cat("=== End IFE Benchmark ===\n")

  # All three must complete — that's the real test
  expect_true(TRUE)
})

test_that("BENCH2: CFE timing comparison (loo vs all_units vs treated_units)", {

  skip_on_cran()

  cat("\n=== Runtime Benchmark: CFE nevertreated cv.method timing ===\n")

  set.seed(42)
  t_loo <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "cfe",
        force           = "two-way",
        time.component.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = "loo",
        se              = FALSE,
        parallel        = FALSE
      )
    ))
  )
  cat(sprintf("  CFE loo:            %6.2f sec (elapsed)\n", t_loo["elapsed"]))

  set.seed(42)
  t_au <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "cfe",
        force           = "two-way",
        time.component.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = "all_units",
        se              = FALSE,
        parallel        = FALSE
      )
    ))
  )
  cat(sprintf("  CFE all_units:      %6.2f sec (elapsed)\n", t_au["elapsed"]))

  set.seed(42)
  t_tu <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "cfe",
        force           = "two-way",
        time.component.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = "treated_units",
        se              = FALSE,
        parallel        = FALSE
      )
    ))
  )
  cat(sprintf("  CFE treated_units:  %6.2f sec (elapsed)\n", t_tu["elapsed"]))

  cat(sprintf("  Speedup (loo/all_units):      %.2fx\n",
              t_loo["elapsed"] / max(t_au["elapsed"], 0.001)))
  cat(sprintf("  Speedup (loo/treated_units):  %.2fx\n",
              t_loo["elapsed"] / max(t_tu["elapsed"], 0.001)))
  cat("=== End CFE Benchmark ===\n")

  expect_true(TRUE)
})

