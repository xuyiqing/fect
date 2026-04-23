# fect 2.2.1

* Fixed `Unsupported bootstrap method: fe` crash when `fect(method = "gsynth", CV = TRUE, se = TRUE, ...)` (or `method = "cfe"`) selected `r.cv = 0` via cross-validation. `fect_nevertreated()` used to relabel its outgoing `$method` to `"fe"` in that case, and `fect_boot()`'s per-iteration dispatcher had no `"fe"` branch. The fix keeps `$method` as `"gsynth"` (or `"cfe"`); the existing dispatcher branches already handle `r = 0` correctly. Reported against the `gsynth` wrapper, which delegates SE to fect.
* Fixed parametric-bootstrap dispatcher so that `one.nonpara` routes to the correct estimator for `ife+notyettreated`, `cfe+nevertreated`, and `cfe+notyettreated` bootstraps. Previously Loop 2 called `fect_nevertreated()` unconditionally regardless of `method` / `time.component.from`, producing incorrect SEs for those combinations. Introduces internal helpers `impute_Y0()` and `valid_controls()`.
* Added a hard gate that errors when `se = TRUE`, `vartype = "parametric"`, and `time.component.from = "notyettreated"` are combined. A coverage simulation showed ~80% empirical coverage versus 95% nominal for `ife+notyettreated+parametric` (EM-induced residual shrinkage, Loop 2 circularity, finite-sample EM loading bias). Users should switch to `time.component.from = "nevertreated"` or `vartype = "bootstrap"` / `"jackknife"` for that combination.

### Parallel CV (Phase 1)

* The `parallel` argument now accepts five forms: `TRUE`, `FALSE`, `"cv"`,
  `"boot"`, and `c("cv", "boot")`. This allows parallelizing CV and bootstrap
  independently. Existing `TRUE`/`FALSE` usage is fully backward-compatible.
* Cross-validation for the IFE model on the `notyettreated` path (`R/cv.R`)
  now supports parallel execution via `future_lapply`. All `(r, fold)` tasks
  are dispatched as a flat task list; the 1% rank-selection rule is applied
  sequentially in the master process, preserving numerical identity with the
  serial path.
* A size threshold (`Nco * TT > 20000`) governs auto-enable when
  `parallel = TRUE`. Use `parallel = "cv"` to override and engage CV
  parallelism on any panel size.
* The future plan is saved and restored on exit (including error paths), so
  scripts that set their own `future::plan()` are unaffected.
* Internal helper functions (`.fect_cv_score_one_*`) extracted into
  `R/cv-helpers.R` as the canonical single implementation of each CV scoring
  step. `fect_nevertreated.R` fold-parallel branches migrated to use these
  helpers, preserving exact numerical output.
* Deferred to later phases: MC parallel CV (Phase 2), CFE parallel CV
  (Phase 3), `boot.R` modernization (Phase 4).

### Parallel CV (Phase 2)

* Cross-validation for the MC (matrix completion) model on the `notyettreated`
  path (`R/cv.R`) now supports parallel execution via `future_lapply`. All
  `(lambda, fold)` tasks are dispatched as a flat task list; the 1% rule for
  lambda selection is applied sequentially in the master.
* In parallel MC mode, all candidate lambda values are evaluated (no early
  `break_check` stopping). In serial mode, the existing `break_check`
  short-circuit is preserved unchanged.

### Parallel CV (Phase 3)

* Cross-validation for the CFE model (complex fixed effects) in the
  `nevertreated` path now uses the same flat `(r, fold)` parallel dispatch
  as IFE and MC, via `future_lapply`. The auto-enable threshold is
  `Nco * TT > 60000`; use `parallel = "cv"` to override.
* The IFE and CFE parallel CV branches in `R/fect_nevertreated.R` have been
  migrated from `foreach %dopar%` (k-fold only) to flat r×k `future_lapply`
  dispatch, improving load balancing when multiple rank candidates are
  evaluated.
* Centralized threshold constants (`.CV_PARALLEL_THRESH`) are now used
  consistently across both `R/cv.R` and `R/fect_nevertreated.R`.
* Future plan lifecycle uses `on.exit(..., after = FALSE)` in both blocks,
  ensuring correct LIFO restoration when IFE and CFE CV run sequentially.
* The `doFuture::registerDoFuture()` call has been removed from both
  `fect_nevertreated.R` parallel setup blocks; `future_lapply` dispatch does
  not require a foreach backend.
* Fix: the prior `foreach %dopar%` parallel CV path in `fect_nevertreated.R`
  (`cv.method = "all_units"`) was non-functional — worker processes could not
  resolve dot-prefix internal helpers by name. The migration to
  `future_lapply(future.packages = "fect")` resolves worker visibility.

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
