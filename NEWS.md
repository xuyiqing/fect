# fect 2.2.1

Parametric-bootstrap fixes (`se = TRUE`, `vartype = "parametric"`):

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
