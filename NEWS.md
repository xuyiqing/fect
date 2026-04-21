# fect 2.2.1

* Fixed `Unsupported bootstrap method: fe` crash when `fect(method = "gsynth", CV = TRUE, se = TRUE, ...)` (or `method = "cfe"`) selected `r.cv = 0` via cross-validation. `fect_nevertreated()` used to relabel its outgoing `$method` to `"fe"` in that case, and `fect_boot()`'s per-iteration dispatcher had no `"fe"` branch. The fix keeps `$method` as `"gsynth"` (or `"cfe"`); the existing dispatcher branches already handle `r = 0` correctly. Reported against the `gsynth` wrapper, which delegates SE to fect.

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
