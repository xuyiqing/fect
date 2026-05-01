## ---------------------------------------------------------------
## Tests for v2.4.2 ci.method extensions and per-type defaults.
##
## v2.4.2 adds two ci.method values --- "bc" (bias-corrected
## percentile) and "normal" (Wald) --- and switches the per-type
## default (NULL trigger): att -> normal, att.cumu -> percentile,
## aptt -> bc, log.att -> bc.
##
## See statsclaw-workspace/fect/ref/v242-vartype-cimethod-design.md
## ---------------------------------------------------------------

suppressWarnings(data("simdata", package = "fect"))


.fit_canonical <- function(nboots = 100, keep.sims = TRUE) {
  set.seed(42)
  suppressWarnings(suppressMessages(
    fect::fect(
      Y ~ D, data = simdata, index = c("id", "time"),
      method = "fe", force = "two-way",
      se = TRUE, nboots = nboots, parallel = FALSE,
      keep.sims = keep.sims
    )
  ))
}

## Positive-Y panel for log.att tests (avoids cell-drop hard-error).
.fit_positive_Y <- function(nboots = 30) {
  set.seed(7)
  N <- 30; TT <- 20
  df <- expand.grid(id = 1:N, time = 1:TT)
  treat_start <- sample(c(NA, 8:15), N, replace = TRUE)
  df$D <- ifelse(is.na(treat_start[df$id]) | df$time < treat_start[df$id],
                 0, 1)
  df$Y <- exp(2.0 + 0.05 * df$time + 0.3 * df$D + rnorm(nrow(df), sd = 0.1))
  set.seed(42)
  suppressWarnings(suppressMessages(
    fect::fect(Y ~ D, data = df, index = c("id", "time"),
               method = "fe", force = "two-way",
               se = TRUE, nboots = nboots, parallel = FALSE,
               keep.sims = TRUE)
  ))
}


## -- CI.1 ci.method enum accepts the four values --------------------

test_that("CI.1: ci.method accepts basic / percentile / bc / normal", {

  skip_on_cran()

  fit <- .fit_canonical()
  for (m in c("basic", "percentile", "bc", "normal")) {
    res <- fect::estimand(fit, "att", "overall", window = c(1, 5),
                          ci.method = m)
    expect_s3_class(res, "data.frame")
    expect_true(!is.na(res$estimate))
    expect_true(!is.na(res$se))
  }
})


## -- CI.2 ci.method = NULL triggers per-type defaults ----------------

test_that("CI.2: ci.method = NULL gives type-specific defaults", {

  skip_on_cran()

  fit <- .fit_canonical()

  ## att -> normal: should match the Wald CI from fit$est.att passthrough
  res_att <- fect::estimand(fit, "att", "event.time")
  ## Compare to explicit normal: should be identical (fast path same numbers)
  res_att_normal <- fect::estimand(fit, "att", "event.time",
                                   ci.method = "normal")
  expect_equal(res_att$estimate, res_att_normal$estimate)
  expect_equal(res_att$se,       res_att_normal$se)
  expect_equal(res_att$ci.lo,    res_att_normal$ci.lo)
  expect_equal(res_att$ci.hi,    res_att_normal$ci.hi)
})


test_that("CI.2b: att overall default = normal differs from basic", {

  skip_on_cran()

  fit <- .fit_canonical()

  res_default <- fect::estimand(fit, "att", "overall", window = c(1, 5))
  res_basic   <- fect::estimand(fit, "att", "overall", window = c(1, 5),
                                ci.method = "basic")
  res_normal  <- fect::estimand(fit, "att", "overall", window = c(1, 5),
                                ci.method = "normal")

  ## Default for att is normal; should equal the explicit normal.
  expect_equal(res_default$ci.lo, res_normal$ci.lo)
  expect_equal(res_default$ci.hi, res_normal$ci.hi)
  ## Basic and normal will generally differ.
  expect_false(isTRUE(all.equal(res_basic$ci.lo, res_normal$ci.lo)))
})


## -- CI.3 normal CI: ci.lo = est - z*SE, ci.hi = est + z*SE ----------

test_that("CI.3: normal CI is symmetric around the point estimate", {

  skip_on_cran()

  fit <- .fit_canonical()
  res <- fect::estimand(fit, "att", "overall", window = c(1, 5),
                        ci.method = "normal")

  z <- stats::qnorm(0.975)
  expect_equal(res$ci.lo, res$estimate - z * res$se, tolerance = 1e-10)
  expect_equal(res$ci.hi, res$estimate + z * res$se, tolerance = 1e-10)
})


## -- CI.4 bc CI: when bootstrap is symmetric around point, bc ~ percentile

test_that("CI.4: bc CI ~ percentile when bootstrap median = point estimate", {

  skip_on_cran()

  ## Use the internal helper with a large symmetric bootstrap. bc reduces
  ## to percentile when z0 = 0 (i.e., bootstrap median = point). With
  ## finite Monte Carlo noise the median isn't exactly the mean, so tol
  ## is loosened to 0.05 (5% of one bootstrap SD).
  set.seed(1)
  boot <- rnorm(5000, mean = 0.5, sd = 0.1)
  ci_pct <- fect:::.compute_ci(estimate = 0.5, boot = boot,
                                 ci.method = "percentile",
                                 conf.level = 0.95)
  ci_bc  <- fect:::.compute_ci(estimate = 0.5, boot = boot,
                                 ci.method = "bc",
                                 conf.level = 0.95)
  expect_equal(ci_pct$ci.lo, ci_bc$ci.lo, tolerance = 0.05)
  expect_equal(ci_pct$ci.hi, ci_bc$ci.hi, tolerance = 0.05)
})


## -- CI.5 bc CI: shifts cutoffs when bootstrap is biased --------------

test_that("CI.5: bc shifts cutoffs when bootstrap median != point estimate", {

  skip_on_cran()

  ## Bootstrap centered at 0.5, but point is at 0.7 (above bootstrap median).
  ## bc should shift cutoffs UP from raw percentile.
  set.seed(2)
  boot <- rnorm(1000, mean = 0.5, sd = 0.1)
  ci_pct <- fect:::.compute_ci(estimate = 0.7, boot = boot,
                                 ci.method = "percentile",
                                 conf.level = 0.95)
  ci_bc  <- fect:::.compute_ci(estimate = 0.7, boot = boot,
                                 ci.method = "bc",
                                 conf.level = 0.95)
  expect_true(ci_bc$ci.lo > ci_pct$ci.lo)
  expect_true(ci_bc$ci.hi > ci_pct$ci.hi)
})


## -- CI.6 vartype column reports method actually used at fit time ----

test_that("CI.6: vartype column reports the fit-time vartype regardless of ci.method", {

  skip_on_cran()

  fit <- .fit_canonical()
  for (m in c("basic", "percentile", "bc", "normal")) {
    res <- fect::estimand(fit, "att", "overall", window = c(1, 5),
                          ci.method = m)
    expect_true(res$vartype %in%
                c("bootstrap", "jackknife", "parametric"))
  }
})


## -- CI.7 cell-drop hard-error fires for log.att on negative-Y data --

test_that("CI.7: log.att on simdata triggers point-level hard-error (v2.4.2+)", {

  skip_on_cran()

  ## simdata has cells with Y <= 0 in the treated post-treatment window;
  ## v2.4.2's point-estimate-level hard-error fires before any bootstrap
  ## work. The bootstrap-level "log-ATT bootstrap is unreliable" message
  ## remains in R/po-estimands.R for strictly-positive panels where
  ## Y0_b crosses zero only in some replicates.
  fit_neg <- .fit_canonical()
  expect_error(
    fect::estimand(fit_neg, "log.att", "event.time"),
    "log\\.att requires Y > 0"
  )
})


test_that("CI.7b: log.att works on a positive-Y panel (no cell-drop pathology)", {

  skip_on_cran()

  fit_pos <- .fit_positive_Y(nboots = 30)
  res <- fect::estimand(fit_pos, "log.att", "event.time")

  expect_s3_class(res, "data.frame")
  expect_true(any(!is.na(res$estimate)))
  expect_true(any(res$n_cells > 0))
})


## -- CI.8 nboots default is 1000 (v2.4.2+) ---------------------------

test_that("CI.8: fect()'s nboots default is 1000", {
  ## Cheap arg-validation; doesn't need a fit. fect.formula and
  ## fect.default are S3 methods, accessed via getS3method().
  expect_equal(as.character(formals(fect::fect)$nboots), "1000")
  expect_equal(as.character(formals(getS3method("fect", "formula"))$nboots),
               "1000")
  expect_equal(as.character(formals(getS3method("fect", "default"))$nboots),
               "1000")
})
