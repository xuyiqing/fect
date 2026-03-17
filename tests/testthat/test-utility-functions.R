## ---------------------------------------------------------------
## Tests for utility functions: fect_mspe,
## esplot (parameter variations), att.cumu, and effect.
## ---------------------------------------------------------------

## Shared fixtures (file-scope) — fitted once, reused across blocks.
suppressWarnings(data("simdata", package = "fect"))

## Basic fect output (no SE, fast)
out_base <- suppressWarnings(suppressMessages(
  fect::fect(
    Y ~ D + X1 + X2,
    data    = simdata,
    index   = c("id", "time"),
    method  = "ife",
    r       = 1,
    CV      = FALSE,
    se      = FALSE,
    parallel = FALSE
  )
))

## Fect output with bootstrap simulations (needed for att.cumu / effect / esplot CI)
out_boot <- suppressWarnings(suppressMessages(
  fect::fect(
    Y ~ D + X1 + X2,
    data      = simdata,
    index     = c("id", "time"),
    method    = "ife",
    r         = 1,
    CV        = FALSE,
    se        = TRUE,
    nboots    = 20,
    keep.sims = TRUE,
    parallel  = FALSE
  )
))

## Subset without treatment reversals — needed for effect() which
## refuses to compute cumulative effects with reversals.
no_rev_ids <- {
  splits <- split(simdata$D, simdata$id)
  as.integer(names(which(vapply(splits, function(x) all(diff(x) >= 0), logical(1)))))
}
simdata_norev <- simdata[simdata$id %in% no_rev_ids, ]

out_norev <- suppressWarnings(suppressMessages(
  fect::fect(
    Y ~ D + X1 + X2,
    data      = simdata_norev,
    index     = c("id", "time"),
    method    = "ife",
    r         = 1,
    CV        = FALSE,
    se        = TRUE,
    nboots    = 20,
    keep.sims = TRUE,
    parallel  = FALSE
  )
))

## -----------------------------------------------------------------
## 1. fect_mspe
## -----------------------------------------------------------------
test_that("fect_mspe returns correct structure and values", {

  ## Basic invocation
  result <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, hide_n = 5, seed = 42)
  ))

  expect_type(result, "list")
  expect_true(all(c("summary", "records", "hide_mask", "fits") %in% names(result)))
  expect_s3_class(result$summary, "data.frame")
  expect_true(all(c("Model", "Hidden_N", "RMSE", "Bias") %in% names(result$summary)))
  expect_s3_class(result$records, "data.frame")
  expect_true(all(c("Rep", "Model", "Hidden_N", "RMSE", "Bias") %in% names(result$records)))
  expect_true(result$summary$RMSE > 0)

  ## hide_mask is a matrix with correct dimensions
  TT <- nrow(out_base$Y.ct.full)
  N  <- ncol(out_base$Y.ct.full)
  expect_true(is.matrix(result$hide_mask))
  expect_equal(dim(result$hide_mask), c(TT, N))

  ## Reproducibility: same seed => same RMSE
  r1 <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, hide_n = 5, seed = 42)
  ))
  r2 <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, hide_n = 5, seed = 42)
  ))
  expect_equal(r1$summary$RMSE, r2$summary$RMSE, tolerance = 1e-10)

  ## Multiple models
  multi <- suppressWarnings(suppressMessages(
    fect_mspe(list(ife = out_base, ife2 = out_base), hide_n = 5, seed = 42)
  ))
  expect_equal(nrow(multi$summary), 2)
  expect_equal(length(multi$fits), 2)
  expect_true(all(c("ife", "ife2") %in% names(multi$fits)))

  ## Multiple reps
  mrep <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, hide_n = 5, n_rep = 2, seed = 123)
  ))
  expect_equal(nrow(mrep$records), 2)
  expect_equal(nrow(mrep$summary), 1)

  ## Custom hide_mask
  custom_mask <- matrix(FALSE, nrow = TT, ncol = N)
  ## Pick a few control-cell positions
  D_idx <- which(names(out_base) == "D")
  D_mat <- NULL
  for (k in rev(D_idx)) {
    obj <- out_base[[k]]
    if (is.matrix(obj) && all(dim(obj) == c(TT, N))) { D_mat <- obj; break }
  }
  ctrl <- which(D_mat == 0, arr.ind = TRUE)
  if (nrow(ctrl) >= 10) custom_mask[ctrl[1:10, , drop = FALSE]] <- TRUE
  res_cm <- suppressWarnings(suppressMessages(
    fect_mspe(out_base, hide_mask = custom_mask, hide_n = 5, seed = 42)
  ))
  expect_type(res_cm, "list")
  expect_true(res_cm$summary$RMSE > 0)
})


## -----------------------------------------------------------------
## 2. esplot parameter variations
## -----------------------------------------------------------------
test_that("esplot handles parameter variations correctly", {

  ## Manual test data
  es_data <- data.frame(
    time = -3:3,
    ATT  = c(-0.1, 0.05, -0.02, 0, 0.5, 0.8, 1.2),
    SE   = rep(0.2, 7),
    count = c(100, 100, 100, 100, 80, 60, 40)
  )

  ## Connected mode with SE-derived CIs
  p1 <- suppressMessages(
    esplot(es_data, Period = "time", Estimate = "ATT", SE = "SE", connected = TRUE)
  )
  expect_s3_class(p1, "ggplot")

  ## SE-derived CI (no CI.lower/CI.upper columns in data)
  p2 <- suppressMessages(
    esplot(es_data, Period = "time", Estimate = "ATT", SE = "SE")
  )
  expect_s3_class(p2, "ggplot")

  ## show.count
  p3 <- suppressMessages(
    esplot(es_data, Period = "time", Estimate = "ATT", SE = "SE",
           Count = "count", show.count = TRUE)
  )
  expect_s3_class(p3, "ggplot")

  ## highlight.periods
  p4 <- suppressMessages(
    esplot(es_data, Period = "time", Estimate = "ATT", SE = "SE",
           highlight.periods = c(-1, 1), highlight.colors = c("red", "blue"))
  )
  expect_s3_class(p4, "ggplot")

  ## fill.gap — create data with a gap (skip time = 0)
  es_gap <- es_data[es_data$time != 0, ]
  p5 <- suppressMessages(
    esplot(es_gap, Period = "time", Estimate = "ATT", SE = "SE", fill.gap = TRUE)
  )
  expect_s3_class(p5, "ggplot")
  ## fill.gap should have added a row for the missing period
  expect_true(nrow(p5$data) > nrow(es_gap))

  ## start0
  p6 <- suppressMessages(
    esplot(es_data, Period = "time", Estimate = "ATT", SE = "SE", start0 = TRUE)
  )
  expect_s3_class(p6, "ggplot")

  ## only.pre
  p7 <- suppressMessages(
    esplot(es_data, Period = "time", Estimate = "ATT", SE = "SE", only.pre = TRUE)
  )
  expect_s3_class(p7, "ggplot")

  ## only.post
  p8 <- suppressMessages(
    esplot(es_data, Period = "time", Estimate = "ATT", SE = "SE", only.post = TRUE)
  )
  expect_s3_class(p8, "ggplot")

  ## From fect output
  p9 <- suppressMessages(
    esplot(out_boot$est.att, Estimate = "ATT")
  )
  expect_s3_class(p9, "ggplot")
})


## -----------------------------------------------------------------
## 4. att.cumu and effect — cumulative treatment effects
## -----------------------------------------------------------------
test_that("att.cumu and effect produce valid cumulative ATT", {

  ## ---- att.cumu (works with reversals) ----
  c1 <- att.cumu(out_boot, period = c(1, 3), plot = FALSE)
  expect_true(is.matrix(c1))
  expect_equal(nrow(c1), 3)
  expect_true("catt" %in% colnames(c1))
  ## All start values equal period[1]
  expect_true(all(c1[, "start"] == 1))
  ## end column goes from 1 to 3
  expect_equal(as.numeric(c1[, "end"]), 1:3)
  ## Values are finite
  expect_true(all(is.finite(c1[, "catt"])))

  ## SE columns present (since out_boot has se=TRUE)
  expect_true(all(c("S.E.", "CI.lower", "CI.upper", "p.value") %in% colnames(c1)))

  ## att.cumu with different period
  c2 <- att.cumu(out_boot, period = c(1, 5), plot = FALSE)
  expect_equal(nrow(c2), 5)

  ## ---- effect (requires no treatment reversals) ----
  ## effect() warns and returns NULL with reversals
  expect_warning(
    e_rev <- effect(out_boot, cumu = TRUE, period = c(1, 3), plot = FALSE),
    "reversals"
  )
  expect_null(e_rev)

  ## Use no-reversal subset for effect() tests
  e1 <- suppressWarnings(
    effect(out_norev, cumu = TRUE, period = c(1, 3), plot = FALSE)
  )
  expect_s3_class(e1, "fect")
  expect_true(!is.null(e1$effect.est.avg))
  expect_equal(length(e1$effect.est.avg), 3)
  expect_true(all(is.finite(e1$effect.est.avg)))
  expect_true(!is.null(e1$effect.est.att))
  expect_true(all(c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value") %in%
                    colnames(e1$effect.est.att)))
  expect_true(all(is.finite(e1$effect.est.att)))

  ## effect non-cumulative
  e2 <- suppressWarnings(
    effect(out_norev, cumu = FALSE, period = c(1, 3), plot = FALSE)
  )
  expect_s3_class(e2, "fect")
  expect_true(!is.null(e2$effect.est.avg))
  expect_true(all(is.finite(e2$effect.est.avg)))

  ## Both att.cumu and effect produce finite output for the same window
  c_norev <- att.cumu(out_norev, period = c(1, 3), plot = FALSE)
  expect_true(all(is.finite(c_norev[, "catt"])))
  expect_true(all(is.finite(e1$effect.est.avg)))

  ## effect with unit subset
  treated_ids <- out_norev$id[colSums(out_norev$D.dat) > 0]
  if (length(treated_ids) >= 5) {
    e3 <- suppressWarnings(
      effect(out_norev, id = treated_ids[1:5], period = c(1, 3), plot = FALSE)
    )
    expect_s3_class(e3, "fect")
    expect_true(!is.null(e3$effect.est.avg))
  }
})
