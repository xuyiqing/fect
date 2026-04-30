<!-- markdownlint-disable MD025 -->
# fect 2.4.1

## New: parametric variance support in `estimand()`

* `estimand()`'s `vartype` argument now accepts `"parametric"` in
  addition to `"bootstrap"`, `"jackknife"`, and `"none"`. When the
  fit was produced with `fect(..., vartype = "parametric",
  keep.sims = TRUE)`, all four `type` values (`"att"`, `"att.cumu"`,
  `"aptt"`, `"log.att"`) work without modification, sourcing
  replicates from the parametric `fit$eff.boot` surface populated
  by the existing fit-time machinery.
* Byte-equality between `estimand(fit, "att", "event.time")` and
  `fit$est.att` is preserved under parametric (asserted by tests).
* The output `vartype` column reports the variance method actually
  used at fit time (read from `fit$vartype`), which may differ from
  the user-supplied `vartype` argument; the argument is informational
  and does not re-aggregate replicates.
* No changes to fit-time machinery, slot semantics, or the v2.4.0
  contract documented in `statsclaw-workspace/fect/ref/po-estimands-contract.md`.

# fect 2.4.0

## New: post-hoc estimands API

* New `estimand(fit, type, by, ...)` typed dispatcher computes
  alternative estimands directly from any fect imputation fit. Shipped
  types: `"att"` (default per-event-time ATT, byte-identical to
  `fit$est.att`), `"att.cumu"` (cumulative ATT, replaces `effect()`),
  `"aptt"` (average proportional treatment effect on the treated;
  Chen & Roth 2024), `"log.att"` (mean log-scale treatment effect).
  Shipped `by` axes: `"event.time"`, `"overall"`, plus reserved
  canonical values `"cohort"` / `"calendar.time"` for future commits.
  Returns a tidy data frame with consistent columns
  `(<by_key>, estimate, se, ci.lo, ci.hi, n_cells, vartype)` regardless
  of `type`. See `?estimand` and the new "Alternative estimands"
  vignette chapter for worked examples and the design rationale.
* New `imputed_outcomes(fit, cells, replicates, direction)` low-level
  accessor returns the cell-level imputed potential-outcome surface as
  a long-form data frame with documented columns
  `(id, time, event.time, cohort, treated, Y_obs, Y0_hat, eff,
  eff_debias, W.agg, [replicate])`. Use this for custom estimands the
  dispatcher does not ship; pipe to dplyr / data.table for arbitrary
  aggregation.
* `cells = ` filter argument on both `imputed_outcomes()` and
  `estimand()` accepts NULL (default), a logical vector, or a one-sided
  formula evaluated against the long-form data
  (e.g. `~ event.time %in% 1:5 & !id %in% bad_ids`). The
  `window = c(L, R)` argument on `estimand()` is sugar over
  `cells = ~ event.time >= L & event.time <= R`.
* `direction = c("on", "off")` on both functions selects the event-time
  grid for reversal panels.
* New `eff_debias` reserved slot on the fit object (NULL for plain
  imputation estimators; populated by future doubly-robust estimators)
  so DR scores can be added to the surface without breaking the
  long-form schema.

## Soft-deprecation: `effect()` and `att.cumu()`

* `effect()` and `att.cumu()` continue to work byte-identically to
  v2.3.x and emit a one-time-per-session message pointing at the
  unified `estimand()` API. Removal not before v3.0.0. Migration:
  - `effect(fit, cumu = TRUE)` → `estimand(fit, "att.cumu", "event.time")`
  - `effect(fit, cumu = FALSE)` → `estimand(fit, "att", "event.time")`
  - `att.cumu(fit, period = c(L, R))` →
    `estimand(fit, "att.cumu", "overall", window = c(L, R))`
  Numerical equality is asserted by package tests.

# fect 2.3.3

## Bug fixes

* Fix `"incorrect number of dimensions"` crash in `diagtest()` that
  surfaced intermittently with `parallel = TRUE` + small `nboots`.
  Two layers: (a) `R/diagtest.R` now uses `drop = FALSE` when filtering
  bootstrap columns by all-non-NA, so a single surviving column stays
  a matrix; (b) the bootstrap parallel path in `R/boot.R` now builds
  the PSOCK cluster via `parallelly::makeClusterPSOCK(rscript_libs =
  .libPaths())` (the same robust pattern used by the CV path),
  wrapped in a 3-attempt retry-with-backoff. If `doParallel` cluster
  init exhausts retries the bootstrap now degrades to sequential
  rather than crashing.
* `carryover.rm` is now stored on the fit object. `plot.fect()` no
  longer reads it from `as.list(x$call)$carryover.rm` (which silently
  gave the wrong K under `do.call()`, programmatic wrappers, or any
  call-rewriting code path). Behaviorally a no-op for fits built via
  named-argument `fect(...)`; correct under any other construction.

## Documentation

* Chapter 2 §Other estimands gains a worked example for post-hoc
  estimands derived from the imputed potential-outcome surface,
  showing APTT (Chen & Roth 2024) with bootstrap CIs from the
  existing fit slots. Issue #126 (ajunquera).

# fect 2.3.2

## Modern visual defaults for `plot.fect()` (visual breaking change)

* Default visual overhaul across all 14 plot types under
  `theme.bw = TRUE`: white panel, plain left-aligned title, thin grey
  reference lines, dashed treatment-onset vline, pre/post lightness
  contrast (`grey50` / `grey20`), publication-sized axis text
  (`cex.main = 11`, `cex.lab = 9`, `cex.axis = 8`, `cex.text = 3.0`),
  and compact legends.
* Placebo / carryover plots render highlighted periods as a single
  accent glyph (orange triangle for placebo, blue diamond for
  carryover, orange triangle for `carryover.rm`) instead of a stacked
  pair of circle + accent. Background rectangle behind each highlight
  period is opt-in via `highlight.fill = TRUE` --- default is glyph
  only, which keeps figures clean for print and grayscale.
* `highlight` argument extended to accept a character subset of
  `c("placebo", "carryover", "carryover.rm")` for selective
  per-test-type highlighting (e.g., `highlight = "placebo"` to
  render carryover periods as plain circles when both tests ran at
  fit time). Backward-compatible: `NULL` / `TRUE` / `FALSE` still
  behave as before.
* Stats annotation block (placebo / carryover / F / equivalence
  p-values) sits at the top-left panel corner with symmetric 2.5%
  inset and `2 × num_stat_lines` top padding so it never grazes the
  leftmost CI. Sized at `cex.text * 1.0` so it does not overpower
  the title. User-supplied `stats.pos` still wins.
* `loadings` (ggpairs) plot: correlation panel reformatted with
  overall + per-group entries; per-group label colors match the
  density-plot fills.
* Migrated off ggplot2 4.0's deprecated `fatten` / `lwd` arguments
  to the `size` / `linewidth` aesthetics. Clears the per-plot
  deprecation warnings.

## `legacy.style = TRUE` escape hatch

New `legacy.style` argument (default `FALSE`). Pass
`legacy.style = TRUE` for byte-identical reproduction of pre-2.3.1
figures (bold centered title, larger axis sizes, solid vline, blue
placebo triangles, no peach rectangle), regardless of `theme.bw`.

## `theme.bw = FALSE` soft-deprecated

Setting `theme.bw = FALSE` now emits a one-time per-session message
flagging removal in v2.5.0. Users who want the gray-panel look
should pass `legacy.style = TRUE` (which honors `theme.bw = FALSE`
exactly), or apply `+ ggplot2::theme_gray()` to the returned plot.

# fect 2.3.1

## New: `W.est` and `W.agg` arguments distinguish survey weights from IPW / balancing weights

`fect()` and `fect.formula()` gain two new arguments that control where
the weight column enters the estimator:

* `W.est` --- weight column for the outcome-model fit (the weighted least
  squares applied inside the IFE / MC / CFE solver).
* `W.agg` --- weight column for the across-treated-obs aggregation
  (`att.on`, `est.avg`, `est.att`).

Both default to `NULL` and fall back to the existing `W` argument when
left unset, so callers who pass only `W = "col"` get the same behavior as
v2.3.0 + the consistency fix below (W enters both fit and aggregation).
Pass the per-role arguments to specify finer behavior:

* Survey / sample weights: `W = "ws"` (or equivalently
  `W.est = W.agg = "ws"`). W enters both fit and aggregation.
* Robust-regression / heteroskedasticity / GLS weights: `W.est = "wr"`
  alone. The fit is weighted, the aggregation is unweighted.
* Inverse-probability / balancing / post-stratification weights:
  `W.agg = "ipw"` alone. The outcome model is fit unweighted (preserving
  doubly-robust properties), and the aggregation is weighted by IPW.

In v2.3.1, `W.est` and `W.agg` (when both supplied) must point to the
same column. Truly distinct columns for fit vs. aggregation (e.g. a
combined survey x IPW design where the outcome model uses survey weights
and the aggregation uses survey x IPW) are scheduled for v2.4.0; the
v2.3.1 design errors with an instructive message if requested.

**Caveat for IPW users.** `W.agg = "ipw"` fits the outcome model
unweighted and applies IPW only at the across-treated-obs aggregation.
This is closer to a doubly-robust estimator than v2.3.0's silent
everywhere-weighting --- but it is not a fully cross-fit doubly-robust
estimator. Residuals on never-treated controls used in any de-bias term
inherit in-sample shrinkage from the outcome fit, which DR theory
requires cross-fitting to eliminate. A fully cross-fit DR path is
scheduled for v3.0.

## Breaking change: weighted fits now have a single, consistent ATT surface

When `W = "<column>"` is supplied to `fect()`, every reported quantity on
the returned fit object now reflects those weights. Prior versions
maintained two parallel ATT pipelines on weighted fits: an unweighted one
(populated into `est.att`, `est.avg`, `att.boot`, `att.vcov`, etc.) and a
W-weighted one (populated into the parallel `est.att.W`, `est.avg.W`,
`att.W.boot`, `att.W.vcov` slots). `plot(fit)` silently substituted the
W-weighted pipeline for rendering while `print(fit)` and `fit$est.att`
returned the unweighted pipeline --- so the same fit object reported
different per-period ATTs and aggregate CIs depending on which surface
the user looked at.

As of 2.3.1, when `W` is non-NULL:

* `fit$est.att`, `fit$est.avg`, `fit$est.att90`, `fit$att.bound`,
  `fit$att.boot`, `fit$att.vcov`, `fit$est.placebo`, `fit$est.carryover`,
  and the `*.off` reverse-treatment counterparts all carry the
  W-weighted aggregations.
* `fit$att`, `fit$time`, `fit$count`, `fit$att.avg`, `fit$att.off`,
  `fit$time.off`, `fit$count.off`, `fit$att.placebo`, `fit$att.carryover`
  similarly carry the W-weighted aggregations from the per-method
  estimator.
* `print(fit)` now labels the obs-level row as `Tr obs sample-weighted (W)`
  (instead of `Tr obs equally weighted`) when W was supplied.
* The redundant `*.W` slots (`est.att.W`, `est.avg.W`, `att.W.boot`,
  `att.W.vcov`, `att.W.bound`, `att.on.W`, `time.on.W`, `count.on.W`,
  `att.avg.W`, `att.on.sum.W`, `W.on.sum`, `att.off.W`, `time.off.W`,
  `count.off.W`, `att.off.sum.W`, `W.off.sum`, `att.placebo.W`,
  `att.carryover.W`, `est.placebo.W`, `est.carryover.W`, `est.att.off.W`,
  `att.off.W.bound`, `att.off.W.vcov`) are no longer present on the
  returned fit object.

If you want the unweighted view of the same fit, refit with `W = NULL`.

The `weight` argument to `plot.fect()` is now a no-op (deprecated),
slated for removal in v2.5.0; passing it emits a deprecation warning.
Internally `plot.fect()` no longer auto-flips between two pipelines ---
it reads the canonical slots, which are already W-weighted when W was
supplied at fit time.

The C++ matrix-completion / IFE / CFE solvers (`inter_fe_mc`,
`inter_fe_ub`, `inter_fe_cfe`) already used W as a fit-time
weighted-least-squares weight in 2.3.0; this release does not change the
fit, only the result-object surface.

# fect 2.3.0

## Rolling-window cross-validation (standard ML design)

* New exported function `r.cv.rolling()`: a standalone user-facing helper
  for picking the number of factors `r` via standard rolling-window CV.
  For each of `k` folds, a fraction `cv.prop` of eligible units (controls
  plus treated pre-treatment) is sampled; only sampled units carry a
  mask in that fold. For each sampled unit, a random anchor time `t*` is
  drawn and the fold's training set excludes:
    1. `cv.nobs` observations starting at `t*` (the held-out, scored block);
    2. `cv.buffer` observations immediately before `t*` (gap buffer, attenuates
       AR-leakage at the past-side train/test boundary --- analogous to
       `cv.donut` for the existing CV strategies, but only on the past side
       since the future side is dropped by construction);
    3. all observations from `t* + cv.nobs` through the unit's
       end-of-eligible (the rolling-window step --- training cannot see
       the future of the held-out block; for treated units, end-of-eligible
       is the cell strictly before treatment onset, so post-treatment
       cells are never masked).
  MSPE is scored at the held-out block only and averaged across folds.
  This is the standard time-series CV design (cf. `forecast::tsCV`,
  `tidymodels::sliding_window`, `caret::createTimeSlices`) adapted to
  panel data.

* Per-fold unit sampling is required: masking every eligible unit at the
  same time would leave no donor data at the masked time points and break
  factor identification at the masked tails. With `cv.prop = 0.2`, every
  eligible unit lands in the holdout roughly `k * cv.prop = 2` times in
  expectation across `k` folds; unsampled units stay fully observed and
  contribute training data at every period.

* New parameters: `cv.buffer` (default 1, past-side buffer length ---
  analogous to `cv.donut` for the existing CV strategies, but only on
  the past side because the future side is dropped by construction);
  `k` (default 10 folds, matching the default for the existing CV
  strategies); `cv.prop` (default 0.2, fraction of eligible units
  sampled per fold; raised from 0.1 after small-panel stability tuning);
  `seed` (optional integer base seed for reproducible per-fold sampling
  and anchor selection).

* Closes the forward-leakage channel that the existing
  `cv.method = "all_units"` / `"treated_units"` (random contiguous-block
  masking) leaves open at `cv.donut = 0 / 1` under serially correlated
  residuals: under rolling-window CV, the train/test boundary on the
  future side is closed by construction.

* Workflow: call
  `r.cv.rolling(formula, data, index, method = "ife", cv.buffer = 1, k = 10, cv.prop = 0.2)`
  to get the chosen `r.cv`, then pass that to
  `fect(..., CV = FALSE, r = r.cv, se = TRUE)` for the inferential fit.

* Identical CV behavior across `method = "ife"` (IFE-EM, internal
  `time.component.from = "notyettreated"`) and `method = "gsynth"` (GSC,
  internal `time.component.from = "nevertreated"`). Both paths populate
  `Y.ct.full` at masked positions: the IFE-EM path via EM imputation,
  the GSC path via the model-implied factor product `F * t(lambda_co)`
  (see "GSC: Y.ct.full populated at control positions" below). Other
  methods (e.g. `"mc"`) are not yet supported.

* Return value is a list with `r.cv`, `cv.rule`, `mspe` (data.frame of
  per-r MSPE averaged across folds, plus fold-SE and held-out cell
  counts), `mspe.per.fold` (r-by-k matrix of per-fold MSPE), and the
  chosen `k`, `cv.nobs`, `cv.buffer`, `cv.prop`. Promoting the design
  to a `cv.method = "rolling"` option inside the main `fect()` CV
  dispatcher is deferred to a future release.

* Empirical motivation: on the Eibl & Hertog (2023) oil-rich panels with
  residual AR(1) of 0.56--0.93, fect's default CV (random anchors,
  `cv.donut` 0 or 1) pegs `r.cv = 5` on every (cell, estimator, rule)
  combination, even at widened `cv.nobs = 6` --- the forward-leakage hides
  the rank overfit. `r.cv.rolling()` recovers ranks consistent with the
  placebo-passing preferred rank for each outcome.

* **Behavior change vs the v2.3.0 development tip's tail-only design**:
  the prior implementation deterministically masked the LAST `cv.nobs`
  observations of each control unit (no folds, no random anchors).
  This was simpler but gave a single MSPE estimate per `r` with no
  fold-to-fold SE. Existing callers' `r.cv.rolling()` invocations will
  produce different numerical results under the new design; the
  selected `r.cv` is typically similar but no longer deterministic for
  a fixed dataset (set `seed` for reproducibility).

## CV API: cv.method = "rolling" in the main dispatcher (additive)

* `fect(CV = TRUE, cv.method = "rolling", ...)` now wires rolling-window
  cross-validation into the main `fect()` CV dispatcher for all factor-model
  methods (ife, cfe, gsynth, mc, both). The same masking logic that powers
  the standalone `r.cv.rolling()` is now available through the standard
  CV API: per-fold sampling of `cv.prop` of eligible units (controls plus
  treated pre-treatment), random anchors per sampled unit, `cv.nobs`-cell
  scored holdout, `cv.buffer`-cell past-side buffer, drop-from-anchor-to-
  end-of-eligible (rolling-window step). Treated post-treatment cells are
  never masked.

* New parameter `cv.buffer` (default 1) controls the past-side buffer for
  rolling CV; replaces the role `cv.donut` plays for block CV. `cv.donut`
  is unchanged for block strategies. Defaults are otherwise unchanged ---
  the existing `cv.method = "all_units"` / `"treated_units"` defaults
  still apply, so this is a fully additive change with no breaking
  defaults.

* `fect_mspe(out, cv.method = "rolling", cv.buffer = 1, ...)` for
  rolling-window-CV-based model comparison.

* `r.cv.rolling(method = "cfe", ...)` extends the standalone helper to
  Complex Fixed Effects. CFE-specific args (Z, gamma, Q, Q.type, kappa,
  extra index columns) are forwarded via `...` and held fixed; rolling
  CV picks `r` only.

* Default `cv.method` flips ("rolling" as the package-wide default),
  `cv.donut → cv.buffer` rename, and `(cv.method, cv.units)` API
  decomposition all remain deferred to a future "CV API unification"
  PR per the plan at statsclaw-workspace/fect/ref/cv-unification-plan.md.

## Bounded factor loadings for GSC

* New argument `loading.bound = "simplex"` (default `"none"`): constrains
  treated-unit factor loadings to the convex hull of control loadings via an
  entropy-regularized simplex projection. Solves, per treated unit `i`:

  ```text
  minimize     (1/gamma) * KL(w || uniform) + || u_pre - F_pre %*% t(Lambda_co) %*% w ||^2
  over w in Delta_{Nco}
  lambda.tr_i = t(Lambda_co) %*% w_i
  ```

  By construction, `Y_hat(0) = F %*% lambda.tr` is a convex combination of
  factor-implied control outcomes, so the counterfactual lies pointwise in
  `conv({Y_hat_co_j})` for every time period. Applies only to
  `method = "ife"` with `time.component.from = "nevertreated"` in this
  version.

* New argument `gamma.loading` (default `NULL`): scalar regularization
  strength for the new `"simplex"` projection. `NULL` triggers 5-fold
  cross-validation over a log-grid. When numeric, `gamma.loading` is used
  directly.

* New argument `gamma.loading.grid` (default `NULL`): user-supplied grid for
  `gamma.loading` CV; `NULL` uses `10^seq(-2, 2, length.out = 9)`.

* The unit-FE scalar `alpha.tr` is NOT bounded; under `loading.bound = "simplex"`
  with `force %in% c("unit", "two-way")`, it is computed as the residual-mean
  `mean(U.tr.pre - F.hat.pre %*% t(lambda.tr))` per treated unit.

* Solver: softmax reparameterization with `stats::optim(method = "L-BFGS-B")`
  and an analytic gradient; mirror-descent fallback on ill-conditioned
  inputs. No new R dependencies.

### Diagnostic outputs (under `loading.bound = "simplex"`)

* `loading.bound`: character, records the setting used.
* `gamma.loading`: the value used (CV-selected or user-supplied).
* `loading.proj.resid`: `Ntr`-vector of `|| U.tr.pre - F.hat.pre %*% lambda.tr ||`
  per treated unit. Values substantially above the control-fit RMSE flag
  treated units lying near or outside `conv(Lambda_co)` (simplex constraint
  binds).

### Semantic change to `wgt.implied` (under `loading.bound = "simplex"` only)

* `wgt.implied` becomes the `Ntr x Nco` simplex-weight matrix: each row sums
  to 1 and is non-negative. This is a direct byproduct of the solver and
  replaces the Moore-Penrose pseudo-inverse representation under the bound.
  When `loading.bound = "none"` (default), `wgt.implied` is unchanged.

### Known caveats

* Percentile bootstrap intervals may under-cover when the simplex constraint
  binds (true treated loading on the boundary of `conv(Lambda_co)`); this is
  the Andrews (1999, 2001) non-standard-limit regime for constrained
  estimators. Detect via `loading.proj.resid`. Boundary-corrected inference
  is deferred to a later release.

* v1 does not support the not-yet-treated IFE dispatch, the MC method, or
  the CFE method. `loading.bound = "simplex"` errors cleanly when combined
  with any of these.

## Parallelism cleanup (Phase A bootstrap)

* Phase A's bootstrap error simulation (`R/boot.R::draw.error`) migrated
  from `foreach %dopar%` to `future.apply::future_lapply`. The old
  `%dopar%` inherited whatever backend was registered globally; after any
  prior parallel fect call, `run_dopar_retry`'s `on.exit` left `doFuture`
  registered, so a subsequent call's Phase A inherited a backend that
  shipped heavy closures per iteration --- producing an ~8x slowdown on
  variant (iii) bootstraps in multi-fit sessions (e.g., a forest plot run).

* The `doFuture::registerDoFuture()` re-registration inside
  `run_dopar_retry`'s `on.exit` was removed; it was the source of the
  global state pollution. The function still falls back to `doParallel`
  if the future backend errors; it just no longer leaves a global
  doFuture registration behind.

* New regression test (`tests/testthat/test-phase-a-future-state.R`):
  asserts two consecutive `fect(parallel = TRUE)` calls in the same R
  process have wall-time ratio < 3x.

## GSC: Y.ct.full populated at control positions

* On the GSC path (`method = "ife"` with
  `time.component.from = "nevertreated"`), `Y.ct.full[, co]` is now
  overwritten with the model-implied factor product `F * t(lambda_co)`
  after the shared `Y.ct.full <- Y.ct` assignment (sourced from
  `est.co.best$factor` and `est.co.best$lambda`, gated on dim agreement
  and non-empty rank). Closes a gap that left `NA` at masked control
  positions because the residual recipe `Y.co - residuals` propagates
  `NA`. Enables user-space rolling CV (`r.cv.rolling()`) for
  `method = "gsynth"`.

* No change to ATT, gap, or `est.avg`: those are computed from
  treated-unit positions and do not consume `Y.ct.full[, co]`. Verified
  on the `simgsynth` anchor (set.seed(11), r=2, force="two-way"):
  ATT.avg unchanged at 4.639593 vs unmodified dev; new control-column
  contents match `F * t(lambda.co)` with max abs diff = 0.

# fect 2.2.1

Parametric-bootstrap fixes (`se = TRUE`, `vartype = "parametric"`):

* Parametric bootstrap on unbalanced panels now correctly reflects
  within-unit serial correlation. Prior versions used a diagonal residual
  covariance in the Gaussian draw, under-estimating ATT standard errors on
  serially correlated data by ≈ √((1+ρ)/(1-ρ)) at AR(1) coefficient ρ.
  Balanced panels and `vartype` ∈ {"bootstrap", "jackknife"} are unaffected.
* Fixed `Unsupported bootstrap method: fe` crash when `method = "gsynth"` or
  `"cfe"` and CV selected `r.cv = 0`. Reported against the `gsynth` wrapper,
  which delegates SE to fect.
* Fixed `one.nonpara` dispatcher so `ife+notyettreated`, `cfe+nevertreated`,
  and `cfe+notyettreated` bootstraps route to the correct Loop-2 estimator;
  previously all three produced incorrect SEs. Introduces internal helpers
  `impute_Y0()` and `valid_controls()`.
* Added hard gate erroring on `ife + notyettreated + parametric`; a coverage
  simulation showed ~80% vs 95% nominal. Use `time.component.from = "nevertreated"`
  or a non-parametric `vartype`.

Parallel cross-validation:

* IFE, MC, and CFE cross-validation now run in parallel via `future_lapply`,
  dispatching `(r, fold)` (or `(lambda, fold)`) flat across workers. Auto-engages
  above `Nco * TT > 20000` for IFE/MC, `> 60000` for CFE.
* The `parallel` argument now accepts five forms: `TRUE`, `FALSE`, `"cv"`,
  `"boot"`, `c("cv", "boot")`. Scalar forms are backward-compatible; string
  forms bypass the auto-threshold.
* MC-only tradeoff: `break_check` short-circuits in serial mode only; parallel
  MC evaluates all candidate lambdas. Prefer `parallel = FALSE` for MC if the
  search typically terminates early.
* The future plan is saved and restored on exit (including error paths), so
  caller-set `future::plan()` is unaffected.
* Internal: `(r, fold)` scoring extracted to `R/cv-helpers.R`;
  `fect_nevertreated.R` parallel CV migrated from `foreach %dopar%` (fold-only)
  to `future_lapply` (flat r × k). The migration also resolves a latent
  worker-visibility issue in the prior `foreach` path.
* Fix: vector `parallel = c("cv", "boot")` no longer errors in `fect_boot`,
  `fit_test`, `permutation`, `fect_sens`, or `did_wrapper` — five sites had
  legacy scalar `parallel == TRUE` / `if (parallel)` checks.

# fect 2.2.0

* Added CFE (Complex Fixed Effects) estimator (`method = "cfe"`)
* Added `time.component.from` parameter for latent factor estimation timing
* Added k-fold cross-validation (`cv.sample`) for nevertreated designs
* Improved plot styling: pre/post shading colors, harmonized `type = "esplot"`
* Fixed EM convergence and solver equivalence issues in CFE routines (C++)
* Fixed bootstrap and parallel setup crashes
* Restructured Quarto book with new CFE chapter

# fect 2.1.1

* Added `codetools` to Imports in DESCRIPTION (required by `trim_closure_env()` in `boot.R`)
* Added `importFrom("utils", "tail")` to NAMESPACE (used in `fect_mspe.R`)
* Bumped version to 2.1.1 and updated Date field for CRAN submission

# fect 2.0.4

* Add new plot `type = "hte"`

# fect 2.0.0

* New syntax
* Merged in **gsynth**

# fect 1.0.0

* First CRAN version
* Fixed bugs

# fect 0.6.5

* Replace fastplm with fixest for fixed effects estimation
* Added plots for heterogeneous treatment effects
* Fixed bugs

# fect 0.4.1

* Added a `NEWS.md` file to track changes to the package.
