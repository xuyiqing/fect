## ---------------------------------------------------------------
## Tests for Phase 1: parallel CV in cv.R IFE notyettreated path
##
## Covers test-spec.md sections P, A, S, T, R, D, L, E, and G
## for REQ-parallel-cv-phase1.
##
## RNG discipline: set.seed() immediately before EVERY fect() call.
## cv.sample() consumes the global RNG for fold masks; any intervening
## RNG call between two fect() calls produces different fold splits.
##
## See: statsclaw-workspace/fect/runs/REQ-parallel-cv-phase1/test-spec.md
## ---------------------------------------------------------------

suppressWarnings(data("simdata", package = "fect"))

## =================================================================
## Section P: Numerical Identity (serial == parallel, IFE notyettreated)
## =================================================================

## -- P.1  IFE x notyettreated x all_units -------------------------

test_that("P.1: IFE notyettreated all_units serial == parallel", {

  skip_on_cran()

  set.seed(42)
  fit_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      cv.method = "all_units",
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  set.seed(42)
  fit_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      cv.method = "all_units",
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  ## r.cv: named integer — must be identical (including the "r" name attribute)
  expect_identical(fit_seq$r.cv, fit_par$r.cv)

  ## CV.out.ife scores: within 1e-10
  cv_diff <- max(abs(fit_seq$CV.out.ife - fit_par$CV.out.ife), na.rm = TRUE)
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out.ife max diff = %.2e (tolerance 1e-10)", cv_diff))

  ## att.avg: same r.cv => same final model => identical
  expect_identical(fit_seq$att.avg, fit_par$att.avg)
})

## -- P.2  IFE x notyettreated x treated_units ---------------------

test_that("P.2: IFE notyettreated treated_units serial == parallel", {

  skip_on_cran()

  set.seed(42)
  fit_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      cv.method = "treated_units",
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  set.seed(42)
  fit_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      cv.method = "treated_units",
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  expect_identical(fit_seq$r.cv, fit_par$r.cv)

  cv_diff <- max(abs(fit_seq$CV.out.ife - fit_par$CV.out.ife), na.rm = TRUE)
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out.ife max diff = %.2e (tolerance 1e-10)", cv_diff))

  expect_identical(fit_seq$att.avg, fit_par$att.avg)
})

## -- P.3  IFE x nevertreated x all_units (helper-migration fixture) --
## Note: simdata triggers the nevertreated early-exit path
## (insufficient pre-treatment records -> r.cv=0, CV.out=NULL).
## The G.1-G.4 tests in test-score-unify.R exercise the substantive
## path with synthetic data; this test pins the migration regression.

test_that("P.3: IFE nevertreated all_units — helper migration preserves output", {

  skip_on_cran()

  fixture <- readRDS(test_path("fixtures", "nevertreated_ife_all.rds"))

  set.seed(42)
  fit_new <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      cv.method = "all_units",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  expect_identical(fit_new$r.cv, fixture$r.cv)

  ## CV.out may be NULL when early-exit fires; handle that gracefully
  if (is.null(fixture$CV.out) && is.null(fit_new$CV.out)) {
    expect_true(TRUE)  ## both NULL: identical
  } else {
    cv_diff <- max(abs(fit_new$CV.out - fixture$CV.out), na.rm = TRUE)
    expect_true(cv_diff < 1e-10,
      info = sprintf("CV.out max diff = %.2e (tolerance 1e-10)", cv_diff))
  }

  expect_identical(fit_new$att.avg, fixture$att.avg)
})

## -- P.4  IFE x nevertreated x treated_units (helper-migration fixture) --

test_that("P.4: IFE nevertreated treated_units — helper migration preserves output", {

  skip_on_cran()

  fixture <- readRDS(test_path("fixtures", "nevertreated_ife_tr.rds"))

  set.seed(42)
  fit_new <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      cv.method = "treated_units",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  expect_identical(fit_new$r.cv, fixture$r.cv)

  if (is.null(fixture$CV.out) && is.null(fit_new$CV.out)) {
    expect_true(TRUE)
  } else {
    cv_diff <- max(abs(fit_new$CV.out - fixture$CV.out), na.rm = TRUE)
    expect_true(cv_diff < 1e-10,
      info = sprintf("CV.out max diff = %.2e (tolerance 1e-10)", cv_diff))
  }

  expect_identical(fit_new$att.avg, fixture$att.avg)
})

## -- P.5  CFE x nevertreated x all_units (helper-migration fixture) --

test_that("P.5: CFE nevertreated all_units — helper migration preserves output", {

  skip_on_cran()

  fixture <- readRDS(test_path("fixtures", "nevertreated_cfe_all.rds"))

  set.seed(42)
  fit_new <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "cfe",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      cv.method = "all_units",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  expect_identical(fit_new$r.cv, fixture$r.cv)

  if (is.null(fixture$CV.out) && is.null(fit_new$CV.out)) {
    expect_true(TRUE)
  } else {
    cv_diff <- max(abs(fit_new$CV.out - fixture$CV.out), na.rm = TRUE)
    expect_true(cv_diff < 1e-10,
      info = sprintf("CV.out max diff = %.2e (tolerance 1e-10)", cv_diff))
  }

  expect_identical(fit_new$att.avg, fixture$att.avg)
})

## -- P.6  CFE x nevertreated x treated_units (helper-migration fixture) --

test_that("P.6: CFE nevertreated treated_units — helper migration preserves output", {

  skip_on_cran()

  fixture <- readRDS(test_path("fixtures", "nevertreated_cfe_tr.rds"))

  set.seed(42)
  fit_new <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "cfe",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      cv.method = "treated_units",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  expect_identical(fit_new$r.cv, fixture$r.cv)

  if (is.null(fixture$CV.out) && is.null(fit_new$CV.out)) {
    expect_true(TRUE)
  } else {
    cv_diff <- max(abs(fit_new$CV.out - fixture$CV.out), na.rm = TRUE)
    expect_true(cv_diff < 1e-10,
      info = sprintf("CV.out max diff = %.2e (tolerance 1e-10)", cv_diff))
  }

  expect_identical(fit_new$att.avg, fixture$att.avg)
})

## =================================================================
## Section A: API Back-Compatibility
## =================================================================

## -- A.1  parallel = TRUE unchanged (threshold gates to serial on simdata) --

test_that("A.1: parallel=TRUE on small panel is backward-compatible (threshold gates to serial)", {

  skip_on_cran()

  set.seed(42)
  fit_ref <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      cv.method = "all_units",
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  ## simdata Nco * TT < 20000 => parallel=TRUE gates to serial => results identical
  set.seed(42)
  fit_true <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      cv.method = "all_units",
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = TRUE
    )
  ))

  expect_identical(fit_ref$r.cv, fit_true$r.cv)

  cv_diff <- max(abs(fit_ref$CV.out.ife - fit_true$CV.out.ife), na.rm = TRUE)
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out.ife max diff = %.2e (tolerance 1e-10)", cv_diff))
})

## -- A.2  parallel = FALSE deterministic ----------------------------

test_that("A.2: parallel=FALSE two identical calls produce identical results", {

  skip_on_cran()

  set.seed(42)
  fit1 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  set.seed(42)
  fit2 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  expect_identical(fit1$r.cv, fit2$r.cv)

  cv_diff <- max(abs(fit1$CV.out.ife - fit2$CV.out.ife), na.rm = TRUE)
  expect_true(cv_diff < 1e-12,
    info = sprintf("CV.out.ife max diff = %.2e (tolerance 1e-12)", cv_diff))
})

## -- A.3  Invalid parallel value raises informative error -----------

test_that("A.3: invalid parallel value raises an informative error", {

  expect_error(
    fect::fect(
      Y ~ D,
      data     = simdata,
      index    = c("id", "time"),
      method   = "ife",
      CV       = TRUE,
      parallel = "all"
    ),
    regexp = "parallel.*must be"
  )

  expect_error(
    fect::fect(
      Y ~ D,
      data     = simdata,
      index    = c("id", "time"),
      method   = "ife",
      CV       = TRUE,
      parallel = 2L
    ),
    regexp = "parallel.*must be"
  )
})

## =================================================================
## Section S: Selective Switches
## =================================================================

## -- S.1  parallel = "cv" engages CV, leaves bootstrap serial ------

test_that("S.1: parallel='cv' emits CV parallel banner but not bootstrap banner", {

  skip_on_cran()

  msgs <- capture.output({
    set.seed(42)
    suppressWarnings(
      fect::fect(
        Y ~ D,
        data      = simdata,
        index     = c("id", "time"),
        method    = "ife",
        r         = 0:2,
        CV        = TRUE,
        k         = 3,
        cv.method = "all_units",
        time.component.from = "notyettreated",
        se        = TRUE,
        nboots    = 5,
        parallel  = "cv",
        cores     = 2
      )
    )
  }, type = "message")

  ## CV banner present
  expect_true(any(grepl("Parallel CV", msgs)),
    info = "Expected 'Parallel CV' banner when parallel='cv'")

  ## Bootstrap parallel banner absent
  expect_false(any(grepl("Parallel computing", msgs)),
    info = "Did not expect 'Parallel computing' banner when parallel='cv' (boot should be serial)")
})

## -- S.2  parallel = "boot" engages bootstrap, leaves CV serial ----

test_that("S.2: parallel='boot' emits bootstrap banner but not CV parallel banner", {

  skip_on_cran()

  msgs <- capture.output({
    set.seed(42)
    suppressWarnings(
      fect::fect(
        Y ~ D,
        data      = simdata,
        index     = c("id", "time"),
        method    = "ife",
        r         = 0:2,
        CV        = TRUE,
        k         = 3,
        cv.method = "all_units",
        time.component.from = "notyettreated",
        se        = TRUE,
        nboots    = 5,
        parallel  = "boot",
        cores     = 2
      )
    )
  }, type = "message")

  ## Bootstrap banner present
  expect_true(any(grepl("Parallel computing", msgs)),
    info = "Expected 'Parallel computing' banner when parallel='boot'")

  ## CV banner absent
  expect_false(any(grepl("Parallel CV", msgs)),
    info = "Did not expect 'Parallel CV' banner when parallel='boot' (CV should be serial)")
})

## =================================================================
## Section T: Threshold Gating and Override
## =================================================================

## -- T.1  parallel = TRUE on small panel -> serial (no CV banner) --

test_that("T.1: parallel=TRUE on small panel runs CV serially (threshold gate)", {

  skip_on_cran()

  ## simdata: Nco * TT well below 20,000 threshold => parallel auto-gates to serial
  msgs <- capture.output({
    set.seed(42)
    suppressWarnings(
      fect::fect(
        Y ~ D,
        data      = simdata,
        index     = c("id", "time"),
        method    = "ife",
        r         = 0:3,
        CV        = TRUE,
        k         = 5,
        cv.method = "all_units",
        time.component.from = "notyettreated",
        se        = FALSE,
        parallel  = TRUE
      )
    )
  }, type = "message")

  expect_false(any(grepl("Parallel CV", msgs)),
    info = "Threshold gate should suppress 'Parallel CV' banner for small simdata panel")
})

## -- T.2  parallel = "cv" on small panel -> parallel engages (override) --

test_that("T.2: parallel='cv' on small panel forces parallel (override)", {

  skip_on_cran()

  msgs <- capture.output({
    set.seed(42)
    suppressWarnings(
      fect::fect(
        Y ~ D,
        data      = simdata,
        index     = c("id", "time"),
        method    = "ife",
        r         = 0:3,
        CV        = TRUE,
        k         = 5,
        cv.method = "all_units",
        time.component.from = "notyettreated",
        se        = FALSE,
        parallel  = "cv",
        cores     = 2
      )
    )
  }, type = "message")

  expect_true(any(grepl("Parallel CV", msgs)),
    info = "Explicit parallel='cv' should override threshold and emit 'Parallel CV' banner")
})

## -- T.3  parallel = "cv" on small panel -> results correct --------

test_that("T.3: parallel='cv' on small panel produces numerically correct output", {

  skip_on_cran()

  set.seed(42)
  fit_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      cv.method = "all_units",
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  set.seed(42)
  fit_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      cv.method = "all_units",
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  expect_identical(fit_seq$r.cv, fit_par$r.cv)

  cv_diff <- max(abs(fit_seq$CV.out.ife - fit_par$CV.out.ife), na.rm = TRUE)
  expect_true(cv_diff < 1e-10,
    info = sprintf("CV.out.ife max diff = %.2e (tolerance 1e-10)", cv_diff))
})

## =================================================================
## Section R: Plan Restoration
## =================================================================

## -- R.1  Plan restored after normal return ------------------------

test_that("R.1: future plan restored after fect() with parallel='cv'", {

  skip_on_cran()

  future::plan(future::sequential)
  plan_before_class <- class(future::plan())

  set.seed(42)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  plan_after_class <- class(future::plan())
  expect_identical(plan_before_class, plan_after_class)

  future::plan(future::sequential)  ## cleanup
})

## -- R.2  Plan restored even when fect errors internally -----------

test_that("R.2: future plan restored even when fect returns an error", {

  skip_on_cran()

  future::plan(future::sequential)
  plan_before_class <- class(future::plan())

  ## min.T0 = 9999L forces "All treated units have been removed" error
  tryCatch(
    suppressMessages(suppressWarnings(
      fect::fect(
        Y ~ D,
        data      = simdata,
        index     = c("id", "time"),
        method    = "ife",
        r         = 0:3,
        CV        = TRUE,
        k         = 5,
        se        = FALSE,
        parallel  = "cv",
        cores     = 2,
        min.T0    = 9999L
      )
    )),
    error = function(e) NULL
  )

  plan_after_class <- class(future::plan())
  expect_identical(plan_before_class, plan_after_class)

  future::plan(future::sequential)  ## cleanup
})

## =================================================================
## Section D: Reproducibility
## =================================================================

## -- D.1  Same seed -> bit-identical parallel results --------------

test_that("D.1: two parallel runs with same seed produce identical CV.out.ife and r.cv", {

  skip_on_cran()

  set.seed(42)
  fit_par1 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      cv.method = "all_units",
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  set.seed(42)
  fit_par2 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      k         = 5,
      cv.method = "all_units",
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  expect_identical(fit_par1$r.cv, fit_par2$r.cv)
  expect_identical(fit_par1$CV.out.ife, fit_par2$CV.out.ife)
})

## =================================================================
## Section L: LOO Path Unaffected
## =================================================================

## -- L.1  LOO path identical regardless of parallel ---------------
## Note: test-spec.md L.1 listed time.component.from = "notyettreated"
## with cv.method = "loo", which is an invalid combination (LOO is only
## supported by fect_nevertreated.R). Per tester audit, L.1 is a
## test-spec error. This test uses time.component.from = "nevertreated"
## (matching G.7 in test-score-unify.R which already covers this path
## substantively with ntdata). With simdata the nevertreated path fires
## an early-exit (insufficient pre-treatment records), so r.cv = 0 and
## CV.out = NULL in both modes — the test confirms no error and identity.

test_that("L.1: LOO path produces identical results for any parallel value (nevertreated)", {

  skip_on_cran()

  set.seed(42)
  fit_loo_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      cv.method = "loo",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  set.seed(42)
  fit_loo_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:3,
      CV        = TRUE,
      cv.method = "loo",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  expect_identical(fit_loo_seq$r.cv, fit_loo_par$r.cv)

  ## CV.out.ife may be NULL (early-exit); handle both cases
  if (is.null(fit_loo_seq$CV.out.ife) && is.null(fit_loo_par$CV.out.ife)) {
    expect_true(TRUE)
  } else {
    cv_diff <- max(abs(fit_loo_seq$CV.out.ife - fit_loo_par$CV.out.ife), na.rm = TRUE)
    expect_true(cv_diff < 1e-12,
      info = sprintf("CV.out.ife max diff = %.2e (tolerance 1e-12)", cv_diff))
  }
})

## =================================================================
## Section E: Edge Cases
## =================================================================

## -- E.1  k = 1 with parallel = "cv" runs serial silently ----------

test_that("E.1: k=1 with parallel='cv' runs serial silently (no error)", {

  skip_on_cran()

  expect_no_error({
    set.seed(42)
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data      = simdata,
        index     = c("id", "time"),
        method    = "ife",
        r         = 0:2,
        CV        = TRUE,
        k         = 1,
        se        = FALSE,
        parallel  = "cv",
        cores     = 2
      )
    ))
  })
})

## -- E.2  cores = 1 with parallel = "cv" runs without error --------

test_that("E.2: cores=1 with parallel='cv' runs without error", {

  skip_on_cran()

  expect_no_error({
    set.seed(42)
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data      = simdata,
        index     = c("id", "time"),
        method    = "ife",
        r         = 0:2,
        CV        = TRUE,
        k         = 3,
        se        = FALSE,
        parallel  = "cv",
        cores     = 1
      )
    ))
  })
})

## -- E.3  parallel = FALSE with cores = 4 ignores cores (no CV banner) --

test_that("E.3: parallel=FALSE with cores=4 runs sequentially (no CV banner)", {

  skip_on_cran()

  msgs <- capture.output({
    set.seed(42)
    suppressWarnings(
      fect::fect(
        Y ~ D,
        data      = simdata,
        index     = c("id", "time"),
        method    = "ife",
        r         = 0:2,
        CV        = TRUE,
        k         = 3,
        se        = FALSE,
        parallel  = FALSE,
        cores     = 4
      )
    )
  }, type = "message")

  expect_false(any(grepl("Parallel CV", msgs)),
    info = "parallel=FALSE should suppress CV parallelism even when cores=4")
})

## -- E.4  c("cv","boot") form is accepted without error ------------

test_that("E.4: parallel=c('cv','boot') is accepted without error", {

  skip_on_cran()

  expect_no_error({
    set.seed(42)
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data      = simdata,
        index     = c("id", "time"),
        method    = "ife",
        r         = 0:2,
        CV        = TRUE,
        k         = 3,
        se        = FALSE,
        parallel  = c("cv", "boot"),
        cores     = 2
      )
    ))
  })
})

## -- E.5  criterion = "pc" with parallel = "cv" — Bug #3 regression guard --
## Bug #3: cv.R parallel branch guard used `criterion != "PC"` (uppercase)
## while the canonical form is `"pc"` (lowercase). Fixed in 379751d.

test_that("E.5: criterion='pc' with parallel='cv' runs without error (Bug #3 guard)", {

  skip_on_cran()

  expect_no_error({
    set.seed(42)
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data      = simdata,
        index     = c("id", "time"),
        method    = "ife",
        r         = 0:3,
        CV        = TRUE,
        k         = 5,
        criterion = "pc",
        se        = FALSE,
        parallel  = "cv",
        cores     = 2
      )
    ))
  })
})

## -- E.6  Bootstrap + CV interaction: se=TRUE, nboots=5, parallel=TRUE — no deadlock --

test_that("E.6: bootstrap + CV interaction (se=TRUE, nboots=5, parallel=TRUE) runs to completion", {

  skip_on_cran()

  ## Defensive: clear any leaked future plan / parallel backend from
  ## upstream tests in the same suite run. Without this, an inherited
  ## doFuture / doParallel registration with stale cluster handles can
  ## cause "incorrect number of dimensions" errors in the bootstrap
  ## phase. Test passes 73/73 in isolation; this guard makes it pass
  ## inside the full devtools::test() ordering as well.
  suppressWarnings({
    try(future::plan(future::sequential), silent = TRUE)
    try(foreach::registerDoSEQ(),         silent = TRUE)
  })

  expect_no_error({
    set.seed(42)
    suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data      = simdata,
        index     = c("id", "time"),
        method    = "ife",
        r         = 0:2,
        CV        = TRUE,
        k         = 3,
        se        = TRUE,
        nboots    = 5,
        parallel  = TRUE
      )
    ))
  })
})

## =================================================================
## Section G: Regression guard — nevertreated parallel paths
## =================================================================
## G.1-G.4 in test-score-unify.R already cover the nevertreated IFE/CFE
## fold-parallel paths with synthetic data (ntdata, N=50 TT=20 Ntr=15)
## that exercises substantive CV scoring (r.cv > 0, CV.out non-NULL).
## No duplication needed here; the P.3-P.6 fixture tests above guard
## against the simdata early-exit case.

## =================================================================
## Section M: MC Parallel CV (Phase 2)
## =================================================================
## Tests for the MC notyettreated parallel branch added in Phase 2.
## Mirror the P-section style for numerical identity and the T/R-section
## style for threshold gating and plan restoration.
##
## RNG discipline: set.seed() immediately before EVERY fect() call.
## simdata is used throughout (small panel; Nco * TT < 20000 threshold).
## parallel = "cv" overrides the threshold for MC explicitly.
## =================================================================

## -- M.1  MC notyettreated serial == parallel (numerical identity) --

test_that("M.1: MC notyettreated serial == parallel (lambda.cv and CV.out.mc within 1e-10)", {

  skip_on_cran()

  set.seed(42)
  fit_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "mc",
      CV        = TRUE,
      k         = 5,
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  set.seed(42)
  fit_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "mc",
      CV        = TRUE,
      k         = 5,
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  ## lambda.cv: plain numeric (no name attribute needed for MC)
  expect_equal(fit_seq$lambda.cv, fit_par$lambda.cv, tolerance = 1e-10)

  ## CV.out.mc: parallel computes all lambdas; serial may break_check early, leaving some
  ## rows at the 1e20 sentinel. Compare only the rows both paths computed (MSPE < 1e19).
  ## The optimal lambda (lambda.cv) is selected from the rows both paths evaluated,
  ## so numerical identity of that selection is the key correctness property.
  both_computed <- fit_seq$CV.out.mc[, "MSPE"] < 1e19 & fit_par$CV.out.mc[, "MSPE"] < 1e19
  if (any(both_computed)) {
    cv_diff <- max(abs(
      fit_seq$CV.out.mc[both_computed, ] - fit_par$CV.out.mc[both_computed, ]
    ), na.rm = TRUE)
    expect_true(cv_diff < 1e-10,
      info = sprintf("CV.out.mc max diff (computed rows only) = %.2e (tolerance 1e-10)", cv_diff))
  }

  ## att.avg: same lambda.cv => same final model => identical
  expect_equal(fit_seq$att.avg, fit_par$att.avg, tolerance = 1e-10)
})

## -- M.2  MC parallel: banner contains "Parallel CV (MC)" string ----

test_that("M.2: parallel='cv' with method='mc' emits MC-specific banner", {

  skip_on_cran()

  msgs <- capture.output({
    set.seed(42)
    suppressWarnings(
      fect::fect(
        Y ~ D,
        data      = simdata,
        index     = c("id", "time"),
        method    = "mc",
        CV        = TRUE,
        k         = 3,
        time.component.from = "notyettreated",
        se        = FALSE,
        parallel  = "cv",
        cores     = 2
      )
    )
  }, type = "message")

  ## Banner contains "Parallel CV (MC):" to distinguish from IFE banner
  expect_true(any(grepl("Parallel CV \\(MC\\)", msgs)),
    info = "Expected 'Parallel CV (MC)' banner when parallel='cv', method='mc'")
})

## -- M.3  MC parallel: threshold gate on small panel ---------------

test_that("M.3: parallel=TRUE on small MC panel runs serially (no MC banner)", {

  skip_on_cran()

  ## simdata: Nco * TT well below 20,000 threshold => parallel=TRUE gates to serial
  msgs <- capture.output({
    set.seed(42)
    suppressWarnings(
      fect::fect(
        Y ~ D,
        data      = simdata,
        index     = c("id", "time"),
        method    = "mc",
        CV        = TRUE,
        k         = 5,
        time.component.from = "notyettreated",
        se        = FALSE,
        parallel  = TRUE
      )
    )
  }, type = "message")

  expect_false(any(grepl("Parallel CV \\(MC\\)", msgs)),
    info = "Threshold gate should suppress MC banner for small simdata panel")
})

## -- M.4  MC parallel: plan restored after normal return ------------

test_that("M.4: future plan restored after MC parallel CV call", {

  skip_on_cran()

  future::plan(future::sequential)
  plan_before_class <- class(future::plan())

  set.seed(42)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "mc",
      CV        = TRUE,
      k         = 5,
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  plan_after_class <- class(future::plan())
  expect_identical(plan_before_class, plan_after_class)

  future::plan(future::sequential)  ## cleanup
})

## -- M.5  method = "both" with parallel = "cv" — LIFO plan restore --

test_that("M.5: method='both' parallel='cv' restores plan correctly (LIFO)", {

  skip_on_cran()

  future::plan(future::sequential)
  plan_before_class <- class(future::plan())

  set.seed(42)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "both",
      r         = 0:2,
      CV        = TRUE,
      k         = 3,
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  plan_after_class <- class(future::plan())
  expect_identical(plan_before_class, plan_after_class)

  future::plan(future::sequential)  ## cleanup
})

## -- M.6  MC reproducibility: two parallel runs bit-identical -------

test_that("M.6: two MC parallel runs with same seed produce identical results", {

  skip_on_cran()

  set.seed(42)
  fit_par1 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "mc",
      CV        = TRUE,
      k         = 5,
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  set.seed(42)
  fit_par2 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "mc",
      CV        = TRUE,
      k         = 5,
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  expect_identical(fit_par1$lambda.cv, fit_par2$lambda.cv)
  expect_identical(fit_par1$CV.out.mc, fit_par2$CV.out.mc)
})

## =================================================================
## Section C: CFE Nevertreated Parallel CV (Phase 3)
## =================================================================
## Tests for the CFE nevertreated flat r×k parallel dispatch added in
## Phase 3. Mirror the M-section structure.
##
## RNG discipline: set.seed() immediately before EVERY fect() call.
## Synthetic data (ntdata_cfe) is used because simdata triggers an
## early-exit (insufficient pre-treatment records). ntdata_cfe is
## generated with a fixed seed inside each test so the test is
## self-contained and reproducible.
##
## See: statsclaw-workspace/fect/runs/REQ-parallel-cv-phase3/spec.md
## =================================================================

## Helper: generate reproducible synthetic panel for CFE tests
## N = 38 (Nco=30, Ntr=8), TT=30, treatment from period 20.
## Nco * TT = 900 < 60000 threshold: use parallel = "cv" to force.
.make_cfe_panel <- function(seed = 1234) {
  set.seed(seed)
  TT <- 30; Nco <- 30; Ntr <- 8; N <- Nco + Ntr
  panel <- expand.grid(id = 1:N, time = 1:TT)
  panel <- panel[order(panel$id, panel$time), ]
  fac1 <- sin(seq(0, pi, length.out = TT))
  fac2 <- cos(seq(0, pi, length.out = TT))
  for (i in 1:N) {
    lambda1_i <- rnorm(1, 0, 0.5); lambda2_i <- rnorm(1, 0, 0.5)
    panel$Y[panel$id == i] <- 2*lambda1_i*fac1 + lambda2_i*fac2 + rnorm(TT, 0, 0.3)
  }
  panel$D <- ifelse(panel$id > Nco & panel$time >= 20, 1, 0)
  panel
}

## -- C.1  CFE nevertreated all_units serial == parallel -------------

test_that("C.1: CFE nevertreated all_units serial == parallel (within 1e-10)", {

  skip_on_cran()

  ntdata_cfe <- .make_cfe_panel()

  set.seed(42)
  fit_seq <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = ntdata_cfe,
      index     = c("id", "time"),
      method    = "cfe",
      r         = 0:2,
      CV        = TRUE,
      k         = 3,
      cv.method = "all_units",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  set.seed(42)
  fit_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = ntdata_cfe,
      index     = c("id", "time"),
      method    = "cfe",
      r         = 0:2,
      CV        = TRUE,
      k         = 3,
      cv.method = "all_units",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  ## r.cv must match
  expect_identical(fit_seq$r.cv, fit_par$r.cv)

  ## CV.out scores within 1e-10
  if (!is.null(fit_seq$CV.out) && !is.null(fit_par$CV.out)) {
    cv_diff <- max(abs(fit_seq$CV.out - fit_par$CV.out), na.rm = TRUE)
    expect_true(cv_diff < 1e-10,
      info = sprintf("CV.out max diff = %.2e (tolerance 1e-10)", cv_diff))
  }

  ## att.avg: same r.cv => same final model
  expect_equal(fit_seq$att.avg, fit_par$att.avg, tolerance = 1e-10)
})

## -- C.2  CFE parallel banner contains "Parallel CV (CFE)" ----------

test_that("C.2: parallel='cv' with method='cfe' emits CFE-specific banner", {

  skip_on_cran()

  ntdata_cfe <- .make_cfe_panel()

  msgs <- capture.output({
    set.seed(42)
    suppressWarnings(
      fect::fect(
        Y ~ D,
        data      = ntdata_cfe,
        index     = c("id", "time"),
        method    = "cfe",
        r         = 0:2,
        CV        = TRUE,
        k         = 3,
        cv.method = "all_units",
        time.component.from = "nevertreated",
        se        = FALSE,
        parallel  = "cv",
        cores     = 2
      )
    )
  }, type = "message")

  ## Banner contains "Parallel CV (CFE):" to distinguish from IFE banner
  expect_true(any(grepl("Parallel CV \\(CFE\\)", msgs)),
    info = "Expected 'Parallel CV (CFE)' banner when parallel='cv', method='cfe'")
})

## -- C.3  CFE threshold gate: parallel=TRUE on small panel -> serial --

test_that("C.3: parallel=TRUE on small CFE panel runs serially (CFE threshold gate)", {

  skip_on_cran()

  ntdata_cfe <- .make_cfe_panel()
  ## ntdata_cfe: Nco * TT = 30 * 30 = 900 < 60000 => parallel=TRUE gates to serial

  msgs <- capture.output({
    set.seed(42)
    suppressWarnings(
      fect::fect(
        Y ~ D,
        data      = ntdata_cfe,
        index     = c("id", "time"),
        method    = "cfe",
        r         = 0:2,
        CV        = TRUE,
        k         = 3,
        cv.method = "all_units",
        time.component.from = "nevertreated",
        se        = FALSE,
        parallel  = TRUE
      )
    )
  }, type = "message")

  expect_false(any(grepl("Parallel CV \\(CFE\\)", msgs)),
    info = "CFE threshold (60k) should suppress CFE parallel banner for small panel")
})

## -- C.4  CFE parallel: plan restored after normal return -----------

test_that("C.4: future plan restored after CFE parallel CV call", {

  skip_on_cran()

  ntdata_cfe <- .make_cfe_panel()

  future::plan(future::sequential)
  plan_before_class <- class(future::plan())

  set.seed(42)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = ntdata_cfe,
      index     = c("id", "time"),
      method    = "cfe",
      r         = 0:2,
      CV        = TRUE,
      k         = 3,
      cv.method = "all_units",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  plan_after_class <- class(future::plan())
  expect_identical(plan_before_class, plan_after_class)

  future::plan(future::sequential)  ## cleanup
})

## -- C.5  CFE reproducibility: two parallel runs bit-identical ------

test_that("C.5: two CFE parallel runs with same seed produce identical results", {

  skip_on_cran()

  ntdata_cfe <- .make_cfe_panel()

  set.seed(42)
  fit_par1 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = ntdata_cfe,
      index     = c("id", "time"),
      method    = "cfe",
      r         = 0:2,
      CV        = TRUE,
      k         = 3,
      cv.method = "all_units",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  set.seed(42)
  fit_par2 <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = ntdata_cfe,
      index     = c("id", "time"),
      method    = "cfe",
      r         = 0:2,
      CV        = TRUE,
      k         = 3,
      cv.method = "all_units",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  expect_identical(fit_par1$r.cv, fit_par2$r.cv)
  if (!is.null(fit_par1$CV.out) && !is.null(fit_par2$CV.out)) {
    expect_identical(fit_par1$CV.out, fit_par2$CV.out)
  }
})

## -- C.6  CFE notyettreated + CV = TRUE => CV forced to FALSE -------
## Architectural confirmation: planner's finding that method="cfe" with
## time.component.from="notyettreated" sets CV=FALSE in default.R.

test_that("C.6: CFE notyettreated+CV=TRUE is silently forced to CV=FALSE (no CV output)", {

  skip_on_cran()

  set.seed(42)
  fit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = simdata,
      index     = c("id", "time"),
      method    = "cfe",
      CV        = TRUE,
      time.component.from = "notyettreated",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  ## CV forced to FALSE => no CV output, r.cv = 0 (single r candidate)
  expect_null(fit$CV.out.ife,
    info = "CV.out.ife should be NULL when CV is forced FALSE for cfe+notyettreated")
  expect_null(fit$CV.out,
    info = "CV.out should be NULL when CV is forced FALSE for cfe+notyettreated")
})

## =================================================================
## Section N: Nevertreated Flat-Dispatch Migration Regression (Phase 3)
## =================================================================
## Pin the post-migration serial == parallel identity for both
## IFE and CFE nevertreated paths using the Phase 3 fixture data.
##
## Fixtures generated from HEAD 5531166 (pre-refactor serial path)
## and stored in tests/testthat/fixtures/phase3-baseline-*.rds.
## Serial path remains identical; parallel path (future_lapply) must
## match within 1e-10.
##
## RNG discipline: set.seed(42) immediately before each fect() call.
## =================================================================

## -- N.1  IFE nevertreated all_units: parallel == serial (fixture) --

test_that("N.1: IFE nevertreated all_units parallel matches phase3 fixture (within 1e-10)", {

  skip_on_cran()

  fixture <- readRDS(test_path("fixtures", "phase3-baseline-ife-all-serial.rds"))
  ntdata <- fixture$panel  ## NOTE: if fixture was saved without $panel, use inline
  ## The fixture was saved with the panel as part of the environment — regenerate inline
  set.seed(if (!is.null(fixture$panel_seed)) fixture$panel_seed else 1234)
  TT <- 30; Nco <- 30; Ntr <- 8; N <- Nco + Ntr
  ntdata <- expand.grid(id = 1:N, time = 1:TT)
  ntdata <- ntdata[order(ntdata$id, ntdata$time), ]
  fac1 <- sin(seq(0, pi, length.out = TT)); fac2 <- cos(seq(0, pi, length.out = TT))
  for (i in 1:N) {
    lambda1_i <- rnorm(1, 0, 0.5); lambda2_i <- rnorm(1, 0, 0.5)
    ntdata$Y[ntdata$id == i] <- 2*lambda1_i*fac1 + lambda2_i*fac2 + rnorm(TT, 0, 0.3)
  }
  ntdata$D <- ifelse(ntdata$id > Nco & ntdata$time >= 20, 1, 0)

  set.seed(42)
  fit_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = ntdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:2,
      CV        = TRUE,
      cv.rule   = "1pct",   # phase3 fixture pinned under legacy 1% rule
      cv.prop   = 0.1,      # phase3 fixture pinned under pre-0.2 default
      cv.donut  = 0,        # phase3 fixture pinned under pre-1 default
      k         = 3,
      cv.method = "all_units",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  expect_identical(fit_par$r.cv, fixture$fit$r.cv)
  att_diff <- max(abs(fit_par$att.avg - fixture$fit$att.avg), na.rm = TRUE)
  expect_true(att_diff < 1e-10,
    info = sprintf("IFE nevertreated all_units parallel vs fixture: att.avg diff = %.2e", att_diff))
})

## -- N.2  IFE nevertreated treated_units: serial still matches fixture --

test_that("N.2: IFE nevertreated treated_units serial matches phase3 fixture (within 1e-10)", {

  skip_on_cran()

  fixture <- readRDS(test_path("fixtures", "phase3-baseline-ife-tr-serial.rds"))
  set.seed(if (!is.null(fixture$panel_seed)) fixture$panel_seed else 1234)
  TT <- 30; Nco <- 30; Ntr <- 8; N <- Nco + Ntr
  ntdata <- expand.grid(id = 1:N, time = 1:TT)
  ntdata <- ntdata[order(ntdata$id, ntdata$time), ]
  fac1 <- sin(seq(0, pi, length.out = TT)); fac2 <- cos(seq(0, pi, length.out = TT))
  for (i in 1:N) {
    lambda1_i <- rnorm(1, 0, 0.5); lambda2_i <- rnorm(1, 0, 0.5)
    ntdata$Y[ntdata$id == i] <- 2*lambda1_i*fac1 + lambda2_i*fac2 + rnorm(TT, 0, 0.3)
  }
  ntdata$D <- ifelse(ntdata$id > Nco & ntdata$time >= 20, 1, 0)

  set.seed(42)
  fit_ser <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = ntdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:2,
      CV        = TRUE,
      cv.rule   = "1pct",   # phase3 fixture pinned under legacy 1% rule
      cv.prop   = 0.1,      # phase3 fixture pinned under pre-0.2 default
      cv.donut  = 0,        # phase3 fixture pinned under pre-1 default
      k         = 3,
      cv.method = "treated_units",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  expect_identical(fit_ser$r.cv, fixture$fit$r.cv)
  att_diff <- max(abs(fit_ser$att.avg - fixture$fit$att.avg), na.rm = TRUE)
  expect_true(att_diff < 1e-10,
    info = sprintf("IFE nevertreated treated_units serial vs fixture: att.avg diff = %.2e", att_diff))
})

## -- N.3  CFE nevertreated all_units: parallel == serial (fixture) --

test_that("N.3: CFE nevertreated all_units parallel matches phase3 fixture (within 1e-10)", {

  skip_on_cran()

  fixture <- readRDS(test_path("fixtures", "phase3-baseline-cfe-all-serial.rds"))
  set.seed(if (!is.null(fixture$panel_seed)) fixture$panel_seed else 1234)
  TT <- 30; Nco <- 30; Ntr <- 8; N <- Nco + Ntr
  ntdata <- expand.grid(id = 1:N, time = 1:TT)
  ntdata <- ntdata[order(ntdata$id, ntdata$time), ]
  fac1 <- sin(seq(0, pi, length.out = TT)); fac2 <- cos(seq(0, pi, length.out = TT))
  for (i in 1:N) {
    lambda1_i <- rnorm(1, 0, 0.5); lambda2_i <- rnorm(1, 0, 0.5)
    ntdata$Y[ntdata$id == i] <- 2*lambda1_i*fac1 + lambda2_i*fac2 + rnorm(TT, 0, 0.3)
  }
  ntdata$D <- ifelse(ntdata$id > Nco & ntdata$time >= 20, 1, 0)

  set.seed(42)
  fit_par <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = ntdata,
      index     = c("id", "time"),
      method    = "cfe",
      r         = 0:2,
      CV        = TRUE,
      cv.rule   = "1pct",   # phase3 fixture pinned under legacy 1% rule
      cv.prop   = 0.1,      # phase3 fixture pinned under pre-0.2 default
      cv.donut  = 0,        # phase3 fixture pinned under pre-1 default
      k         = 3,
      cv.method = "all_units",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  expect_identical(fit_par$r.cv, fixture$fit$r.cv)
  att_diff <- max(abs(fit_par$att.avg - fixture$fit$att.avg), na.rm = TRUE)
  expect_true(att_diff < 1e-10,
    info = sprintf("CFE nevertreated all_units parallel vs fixture: att.avg diff = %.2e", att_diff))
})

## -- N.4  CFE nevertreated treated_units: serial still matches fixture --

test_that("N.4: CFE nevertreated treated_units serial matches phase3 fixture (within 1e-10)", {

  skip_on_cran()

  fixture <- readRDS(test_path("fixtures", "phase3-baseline-cfe-tr-serial.rds"))
  set.seed(if (!is.null(fixture$panel_seed)) fixture$panel_seed else 1234)
  TT <- 30; Nco <- 30; Ntr <- 8; N <- Nco + Ntr
  ntdata <- expand.grid(id = 1:N, time = 1:TT)
  ntdata <- ntdata[order(ntdata$id, ntdata$time), ]
  fac1 <- sin(seq(0, pi, length.out = TT)); fac2 <- cos(seq(0, pi, length.out = TT))
  for (i in 1:N) {
    lambda1_i <- rnorm(1, 0, 0.5); lambda2_i <- rnorm(1, 0, 0.5)
    ntdata$Y[ntdata$id == i] <- 2*lambda1_i*fac1 + lambda2_i*fac2 + rnorm(TT, 0, 0.3)
  }
  ntdata$D <- ifelse(ntdata$id > Nco & ntdata$time >= 20, 1, 0)

  set.seed(42)
  fit_ser <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = ntdata,
      index     = c("id", "time"),
      method    = "cfe",
      r         = 0:2,
      CV        = TRUE,
      k         = 3,
      cv.method = "treated_units",
      cv.rule   = "1pct",   # phase3 fixture pinned under legacy 1% rule
      cv.prop   = 0.1,      # phase3 fixture pinned under pre-0.2 default
      cv.donut  = 0,        # phase3 fixture pinned under pre-1 default
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = FALSE
    )
  ))

  expect_identical(fit_ser$r.cv, fixture$fit$r.cv)
  att_diff <- max(abs(fit_ser$att.avg - fixture$fit$att.avg), na.rm = TRUE)
  expect_true(att_diff < 1e-10,
    info = sprintf("CFE nevertreated treated_units serial vs fixture: att.avg diff = %.2e", att_diff))
})

## -- N.6  method = "both" (sequential IFE + CFE) plan restoration ----
## Formalized from Phase 3 audit spot-check. Runs IFE nevertreated parallel
## CV followed immediately by CFE nevertreated parallel CV. Each call
## registers an on.exit() that restores the caller's sequential plan.
## After both calls complete, the outer plan must be unchanged.

test_that("N.6: plan restored after sequential IFE and CFE nevertreated parallel CV", {

  skip_on_cran()

  ## Nevertreated panel: Nco=30, Ntr=8, TT=30.  Nco*TT=900 < both thresholds;
  ## use parallel="cv" to force parallel dispatch regardless of size.
  set.seed(1234)
  TT <- 30; Nco <- 30; Ntr <- 8; N <- Nco + Ntr
  ntdata <- expand.grid(id = 1:N, time = 1:TT)
  ntdata <- ntdata[order(ntdata$id, ntdata$time), ]
  fac1 <- sin(seq(0, pi, length.out = TT))
  fac2 <- cos(seq(0, pi, length.out = TT))
  for (i in 1:N) {
    lambda1_i <- rnorm(1, 0, 0.5); lambda2_i <- rnorm(1, 0, 0.5)
    ntdata$Y[ntdata$id == i] <- 2*lambda1_i*fac1 + lambda2_i*fac2 + rnorm(TT, 0, 0.3)
  }
  ntdata$D <- ifelse(ntdata$id > Nco & ntdata$time >= 20, 1, 0)

  plan_before <- class(future::plan())

  ## First call: IFE nevertreated parallel CV
  set.seed(42)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = ntdata,
      index     = c("id", "time"),
      method    = "ife",
      r         = 0:2,
      CV        = TRUE,
      k         = 3,
      cv.method = "all_units",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  plan_mid <- class(future::plan())
  expect_equal(plan_before, plan_mid,
    info = "plan should be restored to sequential after IFE nevertreated parallel CV")

  ## Second call: CFE nevertreated parallel CV
  set.seed(42)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = ntdata,
      index     = c("id", "time"),
      method    = "cfe",
      r         = 0:2,
      CV        = TRUE,
      k         = 3,
      cv.method = "all_units",
      time.component.from = "nevertreated",
      se        = FALSE,
      parallel  = "cv",
      cores     = 2
    )
  ))

  plan_after <- class(future::plan())
  expect_equal(plan_before, plan_after,
    info = "plan should be restored to sequential after CFE nevertreated parallel CV")
})
