# Coverage validation suite

Standalone Monte-Carlo studies that validate the parametric / `para.error`
inference path produces nominal-coverage confidence intervals.  These are
**not** unit tests --- they live outside `tests/testthat/` because they take
~30 minutes wall-time and are not appropriate for routine `devtools::test()`
runs.

## When to run

Run before declaring inferential changes done.  Trigger conditions:

- Editing the parametric branch in `R/boot.R` (the `binary == FALSE & method
  %in% c("gsynth", "ife", "cfe") & vartype == "parametric"` branch around
  line 774, or the `para.error` sub-dispatch on AR / empirical / wild)
- Editing the location-shift code in `R/po-estimands.R` (search for
  `is_parametric`)
- Adding or modifying a `vartype`, `ci.method`, or `para.error` value
- Changing the jackknife dispatch or slot contract on `eff.boot`
- Changing how `att.avg.boot` or `eff.boot` is populated

Skip if the change is a refactor / docs / vignette / plotting / CV / no-op
on the bootstrap distribution.  `devtools::test()` is the right gate for
those.

## How to run

```sh
cd /path/to/fect
Rscript tests/coverage-study/run_para_error_coverage.R
```

Output: `/tmp/fect-coverage-study/coverage_para_error_<timestamp>.{csv,rds}`
plus a PASS / FAIL summary printed to stdout.

## Acceptance criteria

| Test | DGP | Reps | nboots | Threshold |
|------|-----|------|--------|-----------|
| T19 | DGP-A (additive TWFE, IID Gaussian, true ATT=3) | 100 | 1000 | coverage in [0.91, 0.99] for every (`para.error`, `ci.method`) cell |
| T20 | DGP-A8 (DGP-A + AR(1) errors at rho=0.8) | 100 | 1000 | coverage >= 0.91 for every cell |
| T21 | DGP-A | 50 | 500 | wild/empirical mean-CI-width ratio in [0.70, 1.30] across all 5 `ci.methods` |

T19 + T20 are 15 cells each (3 `para.error` modes x 5 `ci.methods`); T21 is
5 cells (2 modes x 5 ci.methods, ratios computed within ci.method).

## Why fixed treated block

DGPs hold the treatment indicator $D$ fixed at units 1:Ntr across all
replications.  The parametric pseudo-treated bootstrap targets the
conditional variance

$$V_t = \mathrm{Var}(\widehat{\mathrm{ATT}}_t - \mathrm{ATT}_t \mid \Lambda, F, X, D),$$

and by the law of total variance the marginal variance equals
$\mathbb{E}_{(\Lambda,F)}[V_t] + \mathrm{Var}_{(\Lambda,F)}[b_t]$.
Re-randomizing $D$ across reps adds a $\mathrm{Var}_D[b_t]$ term to the
marginal target the simulation measures coverage against.  Because the
bootstrap calibrates against $V_t$ alone, simulations with random $D$
under-cover by exactly that amount and the under-coverage is not a property
of the bootstrap procedure.  See gsynth-note section 2 for the full
derivation; this design choice mirrors gsynth-note's own MC framework
(footnote on Table 1, and section A.2).

## .Rbuildignore

This directory is excluded from the CRAN tarball via `.Rbuildignore`:

```
^tests/coverage-study$
```

Files here are tracked in git but are not part of the package distribution.
