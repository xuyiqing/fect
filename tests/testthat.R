library(testthat)
library(fect)

## For non-interactive / captured runs (local devtools::test() and GitHub
## Actions), use a MultiReporter that preserves per-test output.
## LocationReporter prints failures as file:line (readable); JunitReporter
## writes test-results.xml that CI tooling surfaces as PR annotations.
## Skipped on CRAN (NOT_CRAN unset, GITHUB_ACTIONS unset) to keep CRAN's
## default behavior unchanged.
if (!interactive() &&
    (identical(Sys.getenv("NOT_CRAN"), "true") ||
     nzchar(Sys.getenv("GITHUB_ACTIONS")))) {
  options(testthat.default_reporter = testthat::MultiReporter$new(list(
    testthat::LocationReporter$new(),
    testthat::JunitReporter$new(file = "test-results.xml")
  )))
}

## On CRAN (NOT_CRAN unset and not on GitHub Actions) run ONLY the minimal
## smoke tests in test-cran.R to keep the check budget small. Locally and in CI
## (NOT_CRAN=true or GITHUB_ACTIONS set) run the full regression suite.
if (!identical(Sys.getenv("NOT_CRAN"), "true") && !nzchar(Sys.getenv("GITHUB_ACTIONS"))) {
  test_check("fect", filter = "^cran$")
} else {
  test_check("fect")
}
