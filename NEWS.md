# fect 2.3.0 (development)

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
