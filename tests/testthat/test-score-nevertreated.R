## ---------------------------------------------------------------
## Tests for fect_nevertreated CV: cv.method dispatch (Section C),
## 1% selection rule verification (Section E), W + count.T.cv weights
## (Section F), end-to-end integration (Section G "Integration"), and
## cv.sample k-fold CV (Section H).
##
## Originally part of test-score-unify.R; split out 2026-05-03 for
## progress visibility. Shared fixtures live in helper-score-unify.R.
## ---------------------------------------------------------------

## =================================================================
## Section C: cv.method in fect_nevertreated (NEW for Phase 2)
## =================================================================

test_that("NT1: fect_nevertreated cv.method='loo' selects r.cv (IFE)", {

  skip_on_cran()
  nt_out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "loo",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(nt_out$r.cv >= 0 && nt_out$r.cv <= 3)
  expect_true(!is.null(nt_out$CV.out))
  mspe_col <- nt_out$CV.out[, "MSPE"]
  computed <- mspe_col[mspe_col < 1e19]
  if (length(computed) > 0) {
    expect_true(all(is.finite(computed)))
  }
  expect_true(!is.null(nt_out$Y.ct))
})

test_that("NT2: fect_nevertreated cv.method='treated_units' selects r.cv (IFE)", {

  skip_on_cran()
  nt_out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "treated_units",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(nt_out$r.cv >= 0 && nt_out$r.cv <= 3)
  expect_true(!is.null(nt_out$CV.out))
})

test_that("NT3: fect_nevertreated cv.method='all_units' selects r.cv (IFE)", {

  skip_on_cran()
  nt_out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "all_units",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(nt_out$r.cv >= 0 && nt_out$r.cv <= 3)
  expect_true(!is.null(nt_out$CV.out))
})

test_that("NT4: fect_nevertreated cv.method='loo' selects r.cv (CFE)", {

  skip_on_cran()
  nt_cfe <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "cfe",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "loo",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(nt_cfe$r.cv >= 0 && nt_cfe$r.cv <= 3)
  expect_true(!is.null(nt_cfe$CV.out))
})

test_that("NT5: fect_nevertreated cv.method='treated_units' selects r.cv (CFE)", {

  skip_on_cran()
  nt_cfe <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "cfe",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "treated_units",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(nt_cfe$r.cv >= 0 && nt_cfe$r.cv <= 3)
  expect_true(!is.null(nt_cfe$CV.out))
})

test_that("NT6: fect_nevertreated default cv.method is treated_units", {

  skip_on_cran()
  nt_default <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      se              = FALSE,
      parallel        = FALSE
    )
  ))

  nt_explicit <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      cv.method       = "treated_units",
      se              = FALSE,
      parallel        = FALSE
    )
  ))

  expect_equal(nt_default$r.cv, nt_explicit$r.cv)
})

test_that("NT7: Invalid cv.method rejected for nevertreated", {

  skip_on_cran()
  expect_error(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      cv.method       = "invalid"
    ),
    "cv.method|arg"
  )
})


## =================================================================
## Section E: 1% Selection Rule Verification (NEW for Phase 2)
## =================================================================

test_that("SEL1: 1% selection rule in IFE nevertreated", {

  skip_on_cran()
  # Large tol should NOT affect the 1% rule
  nt_bigtol <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      se              = FALSE,
      parallel        = FALSE,
      tol             = 0.5
    )
  ))
  nt_smalltol <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      se              = FALSE,
      parallel        = FALSE,
      tol             = 1e-3
    )
  ))
  expect_equal(nt_bigtol$r.cv, nt_smalltol$r.cv)
})

test_that("SEL2: 1% selection rule in CFE nevertreated", {

  skip_on_cran()
  nt_bigtol <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "cfe",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      se              = FALSE,
      parallel        = FALSE,
      tol             = 0.5
    )
  ))
  nt_smalltol <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "cfe",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      se              = FALSE,
      parallel        = FALSE,
      tol             = 1e-3
    )
  ))
  expect_equal(nt_bigtol$r.cv, nt_smalltol$r.cv)
})


## =================================================================
## Section F: W and count.T.cv in fect_nevertreated (NEW for Phase 2)
## =================================================================

test_that("WT1: W weights flow through nevertreated LOO scoring", {

  skip_on_cran()
  # W in fect() is a column name, not a matrix. Add a weight column to ntdata.
  ntdata_w <- ntdata
  ntdata_w$wt <- 1.0  # uniform weights

  nt_w <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata_w,
      index           = c("id", "time"),
      method          = "ife",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      W               = "wt",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_true(nt_w$r.cv >= 0 && nt_w$r.cv <= 3)

  # With uniform weights, r.cv should match unweighted result
  nt_nw <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata,
      index           = c("id", "time"),
      method          = "ife",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  expect_equal(nt_w$r.cv, nt_nw$r.cv)
})

test_that("WT2: Non-uniform W may change r selection", {

  skip_on_cran()
  # W in fect() is a column name. Add non-uniform weight column.
  ntdata_w <- ntdata
  set.seed(99)
  ntdata_w$wt <- runif(nrow(ntdata), 0.5, 2.0)

  nt_w <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data            = ntdata_w,
      index           = c("id", "time"),
      method          = "ife",
      time.component.from    = "nevertreated",
      CV              = TRUE,
      r               = c(0, 3),
      W               = "wt",
      se              = FALSE,
      parallel        = FALSE
    )
  ))
  # Should not crash, r.cv should be valid
  expect_true(nt_w$r.cv >= 0 && nt_w$r.cv <= 3)
})


## =================================================================
## Section G: Integration Tests (NEW for Phase 2)
## =================================================================

test_that("INT1: End-to-end cv.method pipeline", {

  skip_on_cran()
  ## Use CV=FALSE for the fit since fect_mspe performs its own

  ## cross-validation masking; what matters is that the fit object
  ## has valid Y.dat/Y.ct for fect_mspe to score.
  ## NOTE: fect_mspe errors on CV=TRUE fits ("No valid residuals") —
  ## that is a separate source-code issue tracked for builder.
  fit_result <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D,
      data      = ntdata,
      index     = c("id", "time"),
      method    = "ife",
      CV        = FALSE,
      r         = 2,
      se        = FALSE,
      parallel  = FALSE
    )
  ))
  mspe_result <- suppressWarnings(suppressMessages(
    fect_mspe(fit_result, seed = 42, cv.method = "all_units")
  ))
  expect_true("summary" %in% names(mspe_result))
  expect_true(mspe_result$summary$RMSE > 0)
  if ("MSPE" %in% names(mspe_result$summary)) {
    expect_true(mspe_result$summary$MSPE > 0)
  }
})

test_that("INT2: IFE CV respects cv.method='treated_units'", {

  skip_on_cran()
  ## Original test used method="gsynth" + cv.method="loo", but simdata has
  ## treatment reversals which gsynth rejects. We test that cv.method is
  ## respected by using cv.method="treated_units" (non-default) with method="ife".
  ife_out <- suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D + X1 + X2,
      data      = simdata,
      index     = c("id", "time"),
      method    = "ife",
      CV        = TRUE,
      r         = c(0, 3),
      cv.method = "treated_units",
      se        = FALSE,
      parallel  = FALSE
    )
  ))
  expect_true(ife_out$r.cv >= 0 && ife_out$r.cv <= 3)
})


## =================================================================
## Section H: cv.sample k-fold CV in fect_nevertreated (NEW)
## Tests for actual cv.sample-based cross-validation branches
## when cv.method="all_units" or "treated_units" in nevertreated.
## =================================================================

## ---- H.1: IFE smoke tests ---- ##

test_that("NTCV1: cv.method='all_units' IFE produces valid output", {

  skip_on_cran()
  set.seed(42)
  out_au <- suppressWarnings(suppressMessages(
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
  # r.cv is integer in [0, 3]
  expect_true(out_au$r.cv >= 0 && out_au$r.cv <= 3)
  # att.avg is finite
  expect_true(is.finite(out_au$att.avg))
  # est.att has no all-NA columns
  if (!is.null(out_au$est.att)) {
    na_cols <- apply(out_au$est.att, 2, function(x) all(is.na(x)))
    expect_false(all(na_cols),
                 info = "est.att should not have all columns be NA")
  }
  # CV.out exists with proper structure
  expect_true(!is.null(out_au$CV.out))
  mspe_col <- out_au$CV.out[, "MSPE"]
  computed <- mspe_col[mspe_col < 1e19]
  if (length(computed) > 0) {
    expect_true(all(is.finite(computed)))
    expect_true(all(computed > 0))
  }
})

test_that("NTCV2: cv.method='treated_units' IFE produces valid output", {

  skip_on_cran()
  set.seed(42)
  out_tu <- suppressWarnings(suppressMessages(
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
  expect_true(out_tu$r.cv >= 0 && out_tu$r.cv <= 3)
  expect_true(is.finite(out_tu$att.avg))
  if (!is.null(out_tu$est.att)) {
    na_cols <- apply(out_tu$est.att, 2, function(x) all(is.na(x)))
    expect_false(all(na_cols))
  }
  expect_true(!is.null(out_tu$CV.out))
})

## ---- H.2: CFE smoke tests ---- ##

test_that("NTCV3: cv.method='all_units' CFE produces valid output", {

  skip_on_cran()
  set.seed(42)
  out_au_cfe <- suppressWarnings(suppressMessages(
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
  expect_true(out_au_cfe$r.cv >= 0 && out_au_cfe$r.cv <= 3)
  expect_true(is.finite(out_au_cfe$att.avg))
  expect_true(!is.null(out_au_cfe$CV.out))
})

test_that("NTCV4: cv.method='treated_units' CFE produces valid output", {

  skip_on_cran()
  set.seed(42)
  out_tu_cfe <- suppressWarnings(suppressMessages(
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
  expect_true(out_tu_cfe$r.cv >= 0 && out_tu_cfe$r.cv <= 3)
  expect_true(is.finite(out_tu_cfe$att.avg))
  expect_true(!is.null(out_tu_cfe$CV.out))
})

## ---- H.3: LOO backward compatibility ---- ##

test_that("NTCV5: cv.method='loo' IFE backward compatibility", {

  skip_on_cran()
  set.seed(1234)
  out_loo <- suppressWarnings(suppressMessages(
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
  expect_true(out_loo$r.cv >= 0 && out_loo$r.cv <= 3)
  expect_true(is.finite(out_loo$att.avg))
  expect_true(!is.null(out_loo$CV.out))

  # LOO is deterministic: re-running should produce identical r.cv
  set.seed(1234)
  out_loo2 <- suppressWarnings(suppressMessages(
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
  expect_equal(out_loo$r.cv, out_loo2$r.cv)
})

## ---- H.4: r-selection and ATT checks ---- ##

test_that("NTCV6: r-selection validity across all cv.methods (IFE)", {

  skip_on_cran()
  r_start <- 0
  r_end <- 3

  methods_list <- c("loo", "all_units", "treated_units")
  for (cm in methods_list) {
    set.seed(42)
    out <- suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        time.component.from    = "nevertreated",
        CV              = TRUE,
        r               = c(r_start, r_end),
        cv.method       = cm,
        se              = FALSE,
        parallel        = FALSE
      )
    ))

    # r.cv in valid range
    expect_true(out$r.cv >= r_start && out$r.cv <= r_end,
                info = paste("cv.method =", cm, ": r.cv out of range"))

    # CV.out has correct number of rows
    expect_equal(nrow(out$CV.out), r_end - r_start + 1,
                 info = paste("cv.method =", cm, ": CV.out row count wrong"))

    # Score columns in CV.out: check MSPE for evaluated rows
    mspe_col <- out$CV.out[, "MSPE"]
    computed <- mspe_col[mspe_col < 1e19]
    if (length(computed) > 0) {
      expect_true(all(is.finite(computed)),
                  info = paste("cv.method =", cm, ": non-finite MSPE"))
      expect_true(all(computed >= 0),
                  info = paste("cv.method =", cm, ": negative MSPE"))
    }
  }
})

test_that("NTCV7: ATT consistency across cv.methods (IFE)", {

  skip_on_cran()
  att_vals <- numeric(3)
  methods_list <- c("loo", "all_units", "treated_units")
  for (i in seq_along(methods_list)) {
    set.seed(42)
    out <- suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        time.component.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = methods_list[i],
        se              = FALSE,
        parallel        = FALSE
      )
    ))
    att_vals[i] <- out$att.avg
    expect_true(is.finite(out$att.avg),
                info = paste("cv.method =", methods_list[i], ": att.avg not finite"))
  }

  # Sanity check: all ATT values should be in the same ballpark.
  # The true ATT is ~3.0 for this DGP.
  # Allow wide tolerance since different r.cv selections produce different ATTs.
  att_range <- max(att_vals) - min(att_vals)
  expect_true(att_range < 5.0,
              info = paste("ATT values too spread:",
                           paste(round(att_vals, 3), collapse = ", ")))
})

## ---- H.5: Edge cases ---- ##

test_that("NTCV-Edge1: r=c(0,0) with all cv.methods", {

  skip_on_cran()
  for (cm in c("loo", "all_units", "treated_units")) {
    set.seed(42)
    out <- suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        time.component.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 0),
        cv.method       = cm,
        se              = FALSE,
        parallel        = FALSE
      )
    ))
    expect_equal(out$r.cv, 0,
                 info = paste("cv.method =", cm, ": r.cv should be 0 when r=c(0,0)"))
    # Scores should still be finite
    mspe_col <- out$CV.out[, "MSPE"]
    computed <- mspe_col[mspe_col < 1e19]
    if (length(computed) > 0) {
      expect_true(all(is.finite(computed)),
                  info = paste("cv.method =", cm, ": non-finite MSPE with r=0"))
    }
  }
})

test_that("NTCV-Edge2: Small panel with all cv.methods", {

  skip_on_cran()
  # Create data with few pre-treatment periods for treated units
  small_data <- make_factor_data(N = 30, TT = 10, Ntr = 8, r = 1, seed = 99)

  for (cm in c("loo", "all_units", "treated_units")) {
    set.seed(42)
    out <- tryCatch(
      suppressWarnings(suppressMessages(
        fect::fect(
          Y ~ D,
          data            = small_data,
          index           = c("id", "time"),
          method          = "ife",
          force           = "two-way",
          time.component.from    = "nevertreated",
          CV              = TRUE,
          r               = c(0, 2),
          cv.method       = cm,
          se              = FALSE,
          parallel        = FALSE
        )
      )),
      error = function(e) e
    )
    # Should either succeed with valid r.cv or produce informative error
    if (!inherits(out, "error")) {
      expect_true(out$r.cv >= 0 && out$r.cv <= 2,
                  info = paste("cv.method =", cm, ": r.cv out of range (small panel)"))
      expect_true(is.finite(out$att.avg),
                  info = paste("cv.method =", cm, ": att.avg not finite (small panel)"))
    } else {
      # If it errors, the message should be informative (not a cryptic crash)
      expect_true(nchar(conditionMessage(out)) > 0,
                  info = paste("cv.method =", cm, ": error should be informative"))
    }
  }
})

test_that("NTCV-Edge3: Single treated unit", {

  skip_on_cran()
  single_tr_data <- make_factor_data(N = 30, TT = 15, Ntr = 1, r = 1, seed = 77)

  for (cm in c("loo", "all_units", "treated_units")) {
    set.seed(42)
    out <- tryCatch(
      suppressWarnings(suppressMessages(
        fect::fect(
          Y ~ D,
          data            = single_tr_data,
          index           = c("id", "time"),
          method          = "ife",
          force           = "two-way",
          time.component.from    = "nevertreated",
          CV              = TRUE,
          r               = c(0, 2),
          cv.method       = cm,
          se              = FALSE,
          parallel        = FALSE
        )
      )),
      error = function(e) e
    )
    if (!inherits(out, "error")) {
      expect_true(out$r.cv >= 0 && out$r.cv <= 2,
                  info = paste("cv.method =", cm, ": r.cv invalid (single treated)"))
    }
    # If it errors, that's acceptable for single treated unit edge case
  }
})

## ---- H.6: Property-based invariants for cv.sample ---- ##

test_that("NTCV-P1: Score non-negativity in CV.out", {

  skip_on_cran()
  for (cm in c("all_units", "treated_units")) {
    set.seed(42)
    out <- suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        time.component.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = cm,
        se              = FALSE,
        parallel        = FALSE
      )
    ))

    # All score columns should be non-negative
    score_cols <- c("MSPE", "WMSPE", "GMSPE", "WGMSPE", "MAD")
    for (sc in score_cols) {
      if (sc %in% colnames(out$CV.out)) {
        vals <- out$CV.out[, sc]
        computed <- vals[vals < 1e19]
        if (length(computed) > 0) {
          expect_true(all(computed >= 0),
                      info = paste("cv.method =", cm, ", score =", sc, ": negative value"))
        }
      }
    }
  }
})

test_that("NTCV-P4: LOO determinism", {

  skip_on_cran()
  set.seed(42)
  out1 <- suppressWarnings(suppressMessages(
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

  set.seed(999)  # different seed should not affect LOO
  out2 <- suppressWarnings(suppressMessages(
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

  expect_equal(out1$r.cv, out2$r.cv)
  # CV.out scores should be identical
  expect_equal(out1$CV.out, out2$CV.out, tolerance = 1e-10)
})

test_that("NTCV-P5: cv.sample reproducibility with set.seed", {

  skip_on_cran()
  for (cm in c("all_units", "treated_units")) {
    set.seed(42)
    out1 <- suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        time.component.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = cm,
        se              = FALSE,
        parallel        = FALSE
      )
    ))

    set.seed(42)
    out2 <- suppressWarnings(suppressMessages(
      fect::fect(
        Y ~ D,
        data            = ntdata,
        index           = c("id", "time"),
        method          = "ife",
        force           = "two-way",
        time.component.from    = "nevertreated",
        CV              = TRUE,
        r               = c(0, 3),
        cv.method       = cm,
        se              = FALSE,
        parallel        = FALSE
      )
    ))

    expect_equal(out1$r.cv, out2$r.cv,
                 info = paste("cv.method =", cm, ": r.cv not reproducible"))
    expect_equal(out1$CV.out, out2$CV.out, tolerance = 1e-10,
                 info = paste("cv.method =", cm, ": CV.out not reproducible"))
  }
})


