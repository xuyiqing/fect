# fect 2.2.0

* Added CFE (Complex Fixed Effects) estimator as `method = "cfe"`
* Fixed `est.fect$sigma2` crash when bootstrap is enabled with normalization
* Fixed `detectCores()` returning NA in parallel setup causing worker count failures
* Updated `index` parameter documentation to describe CFE grouping variables

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
