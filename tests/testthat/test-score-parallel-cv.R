## ---------------------------------------------------------------
## Tests for parallel CV folds in fect_nevertreated.
## Section originally also labeled "G" (now disambiguated via filename
## as test-score-parallel-cv.R; the integration "G" is in
## test-score-nevertreated.R).
##
## Originally Section G "Parallel CV Folds" of test-score-unify.R;
## split out 2026-05-03. Shared fixtures live in helper-score-unify.R.
## ---------------------------------------------------------------

## ---------------------------------------------------------------
## Section G: Parallel CV Folds in fect_nevertreated
##
## Tests for REQ-parallel-cv: verifies that parallel=TRUE produces
## identical results to parallel=FALSE (sequential), that
## reproducibility holds under parallelism with fixed seeds, that
## default behavior is unchanged, that the LOO path is unaffected,
## and that edge cases work correctly.
##
## Follows test-spec.md for REQ-parallel-cv.
## Tolerances: 1e-10 for CV score differences (per test-spec.md).
## ---------------------------------------------------------------

## -- G.1  Sequential-Parallel Equivalence: IFE, all_units ----------

test_that("G.1: parallel CV matches sequential — IFE, all_units", {

  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  set.seed(123)
  result_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2
    )
  ))

  ## r.cv must be identical
  expect_identical(result_seq$r.cv, result_par$r.cv)

  ## CV.out matrix must match within tolerance 1e-10
  cv_diff <- max(abs(result_seq$CV.out - result_par$CV.out))
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out max diff = %.2e (tolerance = 1e-10)", cv_diff))
})

## -- G.2  Sequential-Parallel Equivalence: IFE, treated_units ------

test_that("G.2: parallel CV matches sequential — IFE, treated_units", {

  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "treated_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  set.seed(123)
  result_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "treated_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2
    )
  ))

  expect_identical(result_seq$r.cv, result_par$r.cv)

  cv_diff <- max(abs(result_seq$CV.out - result_par$CV.out))
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out max diff = %.2e (tolerance = 1e-10)", cv_diff))
})

## -- G.3  Sequential-Parallel Equivalence: CFE, all_units ----------

test_that("G.3: parallel CV matches sequential — CFE, all_units", {

  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "cfe",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  set.seed(123)
  result_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "cfe",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2
    )
  ))

  expect_identical(result_seq$r.cv, result_par$r.cv)

  cv_diff <- max(abs(result_seq$CV.out - result_par$CV.out))
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out max diff = %.2e (tolerance = 1e-10)", cv_diff))
})

## -- G.4  Sequential-Parallel Equivalence: CFE, treated_units ------

test_that("G.4: parallel CV matches sequential — CFE, treated_units", {

  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "cfe",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "treated_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  set.seed(123)
  result_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "cfe",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "treated_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2
    )
  ))

  expect_identical(result_seq$r.cv, result_par$r.cv)

  cv_diff <- max(abs(result_seq$CV.out - result_par$CV.out))
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out max diff = %.2e (tolerance = 1e-10)", cv_diff))
})

## -- G.5  Reproducibility Under Parallelism ------------------------

test_that("G.5: parallel CV is reproducible with same seed", {

  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_par1 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2,
      seed         = 12345
    )
  ))

  set.seed(123)
  result_par2 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2,
      seed         = 12345
    )
  ))

  expect_identical(result_par1$r.cv, result_par2$r.cv)
  expect_identical(result_par1$CV.out, result_par2$CV.out)
})

## -- G.6  Default Behavior Unchanged -------------------------------

test_that("G.6: default (no parallel arg) behaves as sequential", {

  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  ## Call without specifying parallel or cores — should default to sequential
  set.seed(123)
  result_default <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE
    )
  ))

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  expect_identical(result_default$r.cv, result_seq$r.cv)
  expect_identical(result_default$CV.out, result_seq$CV.out)
})

## -- G.7  LOO Path Unaffected --------------------------------------

test_that("G.7: LOO path is unaffected by parallel flag", {

  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_loo_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      cv.method    = "loo",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  set.seed(123)
  result_loo_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      cv.method    = "loo",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2
    )
  ))

  ## LOO should produce identical results regardless of parallel flag
  expect_identical(result_loo_seq$r.cv, result_loo_par$r.cv)
  expect_identical(result_loo_seq$CV.out, result_loo_par$CV.out)
})

## -- G.8  Edge Case: k = 1 ----------------------------------------

test_that("G.8: edge case — k = 1 with parallel=TRUE runs without error", {

  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 1,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2
    )
  ))

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 1,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  expect_identical(result_par$r.cv, result_seq$r.cv)
})

## -- G.9  Edge Case: cores = 1 ------------------------------------

test_that("G.9: edge case — cores = 1 behaves as sequential", {

  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_c1 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 1
    )
  ))

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  cv_diff <- max(abs(result_c1$CV.out - result_seq$CV.out))
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out max diff = %.2e (tolerance = 1e-10)", cv_diff))
})

## -- G.10  Edge Case: cores = NULL (auto-detect) ------------------

test_that("G.10: edge case — cores = NULL auto-detects and runs", {

  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_auto <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = NULL
    )
  ))

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  cv_diff <- max(abs(result_auto$CV.out - result_seq$CV.out))
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out max diff = %.2e (tolerance = 1e-10)", cv_diff))
})

## -- G.11  Edge Case: parallel=FALSE with cores specified ----------

test_that("G.11: edge case — parallel=FALSE ignores cores", {

  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  set.seed(123)
  result_no_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE,
      cores        = 4
    )
  ))

  set.seed(123)
  result_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = FALSE
    )
  ))

  expect_identical(result_no_par$r.cv, result_seq$r.cv)
  expect_identical(result_no_par$CV.out, result_seq$CV.out)
})

## -- G.12  Timing Benchmark: Parallel vs Sequential ----------------
## This is an informational benchmark; we report speedup but do not
## hard-fail on specific speedup thresholds (per test-spec.md:
## "Not a pass/fail test, but auditor should measure and report").

test_that("G.12: timing benchmark — parallel vs sequential with 10 cores", {

  skip_on_cran()

  ## Larger dataset for meaningful timing differences
  dat <- make_factor_data(N = 100, TT = 30, Ntr = 25, r = 2, seed = 99)

  cat("\n=== Parallel CV Timing Benchmark ===\n")

  ## --- IFE, all_units ---
  set.seed(123)
  t_seq_ife_au <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "ife",
        force        = "two-way",
        time.component.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "all_units",
        se           = FALSE,
        parallel     = FALSE
      )
    ))
  )

  set.seed(123)
  t_par_ife_au <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "ife",
        force        = "two-way",
        time.component.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "all_units",
        se           = FALSE,
        parallel     = TRUE,
        cores        = 10
      )
    ))
  )

  cat(sprintf("  IFE all_units  sequential: %6.2f sec\n", t_seq_ife_au["elapsed"]))
  cat(sprintf("  IFE all_units  parallel:   %6.2f sec (10 cores)\n", t_par_ife_au["elapsed"]))
  speedup_ife_au <- t_seq_ife_au["elapsed"] / max(t_par_ife_au["elapsed"], 0.001)
  cat(sprintf("  Speedup: %.2fx\n\n", speedup_ife_au))

  ## --- IFE, treated_units ---
  set.seed(123)
  t_seq_ife_tu <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "ife",
        force        = "two-way",
        time.component.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "treated_units",
        se           = FALSE,
        parallel     = FALSE
      )
    ))
  )

  set.seed(123)
  t_par_ife_tu <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "ife",
        force        = "two-way",
        time.component.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "treated_units",
        se           = FALSE,
        parallel     = TRUE,
        cores        = 10
      )
    ))
  )

  cat(sprintf("  IFE treated_units  sequential: %6.2f sec\n", t_seq_ife_tu["elapsed"]))
  cat(sprintf("  IFE treated_units  parallel:   %6.2f sec (10 cores)\n", t_par_ife_tu["elapsed"]))
  speedup_ife_tu <- t_seq_ife_tu["elapsed"] / max(t_par_ife_tu["elapsed"], 0.001)
  cat(sprintf("  Speedup: %.2fx\n\n", speedup_ife_tu))

  ## --- CFE, all_units ---
  set.seed(123)
  t_seq_cfe_au <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "cfe",
        force        = "two-way",
        time.component.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "all_units",
        se           = FALSE,
        parallel     = FALSE
      )
    ))
  )

  set.seed(123)
  t_par_cfe_au <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "cfe",
        force        = "two-way",
        time.component.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "all_units",
        se           = FALSE,
        parallel     = TRUE,
        cores        = 10
      )
    ))
  )

  cat(sprintf("  CFE all_units  sequential: %6.2f sec\n", t_seq_cfe_au["elapsed"]))
  cat(sprintf("  CFE all_units  parallel:   %6.2f sec (10 cores)\n", t_par_cfe_au["elapsed"]))
  speedup_cfe_au <- t_seq_cfe_au["elapsed"] / max(t_par_cfe_au["elapsed"], 0.001)
  cat(sprintf("  Speedup: %.2fx\n\n", speedup_cfe_au))

  ## --- CFE, treated_units ---
  set.seed(123)
  t_seq_cfe_tu <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "cfe",
        force        = "two-way",
        time.component.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "treated_units",
        se           = FALSE,
        parallel     = FALSE
      )
    ))
  )

  set.seed(123)
  t_par_cfe_tu <- system.time(
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data         = dat,
        index        = c("id", "time"),
        method       = "cfe",
        force        = "two-way",
        time.component.from = "nevertreated",
        CV           = TRUE,
        r            = c(0, 5),
        k            = 10,
        cv.method    = "treated_units",
        se           = FALSE,
        parallel     = TRUE,
        cores        = 10
      )
    ))
  )

  cat(sprintf("  CFE treated_units  sequential: %6.2f sec\n", t_seq_cfe_tu["elapsed"]))
  cat(sprintf("  CFE treated_units  parallel:   %6.2f sec (10 cores)\n", t_par_cfe_tu["elapsed"]))
  speedup_cfe_tu <- t_seq_cfe_tu["elapsed"] / max(t_par_cfe_tu["elapsed"], 0.001)
  cat(sprintf("  Speedup: %.2fx\n\n", speedup_cfe_tu))

  cat("=== End Parallel CV Timing Benchmark ===\n")

  ## Informational — always passes; speedup is reported in test output
  expect_true(TRUE)
})

## -- G.13  Property: Backend Cleanup After Parallel CV -------------

test_that("G.13: parallel backend is restored after fect() returns", {

  skip_on_cran()

  dat <- make_factor_data(N = 50, TT = 20, Ntr = 15, r = 2, seed = 42)

  old_plan <- future::plan()

  set.seed(123)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data         = dat,
      index        = c("id", "time"),
      method       = "ife",
      force        = "two-way",
      time.component.from = "nevertreated",
      CV           = TRUE,
      r            = c(0, 3),
      k            = 5,
      cv.method    = "all_units",
      se           = FALSE,
      parallel     = TRUE,
      cores        = 2
    )
  ))

  new_plan <- future::plan()

  ## The future plan class should be restored
  expect_identical(class(old_plan), class(new_plan))
})
