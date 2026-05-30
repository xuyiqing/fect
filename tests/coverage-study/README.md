# Coverage validation suite

Standalone Monte-Carlo studies that validate the v2.4.2 inference paths
(`bootstrap`, `parametric`, `jackknife`) produce nominal-coverage
confidence intervals.  These are **not** unit tests --- they live outside
`tests/testthat/` because they take more wall-time than is appropriate for
routine `devtools::test()` runs.

## Two scripts

- **`run_minimal_coverage.R`** --- routine pre-merge gate.  Four
  scenarios, ~1.5 min wall-time on 16 cores at K=200, nboots=200.  This
  is what you should run for any inference-relevant change.
- **`run_para_error_coverage.R`** --- extended characterization of the
  parametric / `para.error` machinery on a small panel.  ~30-50 min wall
  time at cores=10.  Run when modifying the `para.error` dispatch or the
  parametric pseudo-treated bootstrap specifically.

A small follow-up script `run_minimal_coverage_tail_rerun.R` re-runs at
nboots=1000 only those scenarios whose tail-CI methods (basic /
percentile / bc / bca) came in below 0.93 in the K=200 / nboots=200
default --- consistent with the `.check_tail_ci_replicates` warning
gate (Efron 1987 §3, DiCiccio & Efron 1996 §4 recommend B >= 1000 for
tail-quantile CIs).  Normal CI is unaffected by B; jackknife is normal
only.

## When to run

Run `run_minimal_coverage.R` before declaring any of these changes done:

- Editing `R/boot.R` (any of the bootstrap / parametric / jackknife /
  cluster-bootstrap branches)
- Editing `R/po-estimands.R` (location-shift code, `.compute_ci`,
  jackknife dispatch)
- Adding or modifying a `vartype`, `ci.method`, or `para.error` value
- Changing how `att.avg.boot`, `eff.boot`, `att.avg.unit.boot`, or the
  `est.avg` row is populated

Skip if the change is a pure refactor / docs / vignette / plotting / CV
/ no-op on the bootstrap distribution.  `devtools::test()` is the right
gate for those.

## How to run

```sh
cd /path/to/fect
Rscript tests/coverage-study/run_minimal_coverage.R
```

Output: `/tmp/fect-coverage-study/minimal_<scenario>_K<K>_nb<nb>_<ts>.csv`
(per-rep raw) plus `minimal_summary_K<K>_nb<nb>_<ts>.csv` (per-cell
coverage / MC SE / mean-CI-width / SE / empirical SD / mean bias).

Conditional follow-up (only if the K=200 / nboots=200 run shows tail-CI
cells below 0.93):

```sh
Rscript tests/coverage-study/run_minimal_coverage_tail_rerun.R
```

## Acceptance --- minimal coverage

Four scenarios; coverage is nominal at K=200, MC SE about 0.015 around
0.95.

| Scenario | DGP | Inference | Cells | Threshold |
|----------|-----|-----------|-------|-----------|
| A | factor (r=2), IID, gsynth-note Xu-2017 spec    | `vartype = "parametric"`, `para.error = "auto"` (-> `empirical`), all 5 ci.methods | 5 | coverage >= 0.93 (1.4 SE below 0.95 at K=200) |
| B | factor (r=2), AR(1) rho=0.8                     | `vartype = "parametric"`, `para.error = "auto"` (-> `ar`),        all 5 ci.methods | 5 | coverage >= 0.93                              |
| C1 | additive TWFE (r=0), AR(1) rho=0.5             | `vartype = "bootstrap"` (cluster), all 5 ci.methods                                | 5 | coverage >= 0.93 at nboots = 1000             |
| C2 | additive TWFE (r=0), AR(1) rho=0.5             | `vartype = "jackknife"`, `ci.method = "normal"` only (E&T 1993 §11)                | 1 | coverage >= 0.91                              |

Tail-CI methods (`basic`, `percentile`, `bc`, `bca`) under-cover at
nboots = 200 because the relevant order statistics of the bootstrap
distribution are unstable at small B (Efron 1987 §3).  The follow-up
script reruns failing scenarios at nboots = 1000, which restores nominal
coverage.

C2 jackknife normal-only is a deliberate restriction: jackknife produces
an SE estimate via the Tukey pseudo-value formula, not a sampling
distribution, so reflection-based methods (basic / percentile),
bias-corrected methods (bc / bca), and the BCa acceleration parameter
all lack a defined input.  See `.check_jackknife_ci_method` in
`R/po-estimands.R` for the full hard-error message.

## DGP details (minimal)

A and B replicate gsynth-note's `code/sims/coverage/simulate-xu-{iid,ar1}-rfit2.R`:

- N_tr = 5, N_co = 50, T = 30, T0 = 20 (10 post-periods)
- r = 2 latent factors, no covariates
- λ_i ~ U(-√3, √3); F_t ~ N(0, 1); α_i ~ U(-√3, √3); ξ_t ~ N(0, 1); μ = 5
- λ, F, α, ξ redrawn each rep
- ATT_t = t for t = 1..10 plus N(0, D.sd = 1) per (unit, post-time)
- Coverage target = realized average treated-post effect (within rep)
- A errors: IID N(0, 1).  B errors: AR(1) ρ = 0.8, marginal variance 1
- Estimator: `method = "ife", r = 2, force = "two-way",
  time.component.from = "nevertreated", CV = FALSE` (Xu Alg 2 path)

C scales N for stable cluster-bootstrap and jackknife inference:

- N_tr = 20, N_co = 80, T = 30, T0 = 20
- No factors (r = 0); α_i ~ U(-√3, √3); ξ_t ~ N(0, 1); μ = 5
- ATT = 3 constant (D.sd = 0); coverage target = 3
- Errors: AR(1) ρ = 0.5, marginal variance 1
- Estimator: `method = "fe", force = "two-way",
  time.component.from = "notyettreated"`

## Parallelization

Outer-loop parallelism via `future::plan(future::multisession,
workers = 16)` and `future.apply::future_lapply`.  Each `fect()` call
runs sequentially (`parallel = FALSE`).  Reps are embarrassingly
parallel; this avoids the per-rep cluster fork/join overhead that inner
bootstrap parallelism pays K times.

`future.seed = TRUE` produces L'Ecuyer streams across workers for
reproducible RNG independent of rep ordering.

## Acceptance --- extended (`run_para_error_coverage.R`)

Three legacy tests on a small DGP-A panel (N=40, T=20, T0=12, IID or
AR(1) ρ=0.8).  Now superseded by the minimal suite for routine gating;
retained for deep characterization of the `para.error` dispatch.

| Test | DGP | Reps | nboots | Threshold |
|------|-----|------|--------|-----------|
| T19 | DGP-A (additive TWFE, IID Gaussian, ATT=3) | 100 | 1000 | coverage >= 0.90 for every (`para.error`, `ci.method`) cell |
| T20 | DGP-A8 (DGP-A + AR(1) ρ=0.8)               | 100 | 1000 | coverage >= 0.91 for every cell                              |
| T21 | DGP-A                                       |  50 |  500 | wild/empirical mean-CI-width ratio in [0.70, 1.30]           |

T19 / T20 are 15 cells each (3 `para.error` modes × 5 `ci.methods`); T21
is 5 cells (2 modes × 5 ci.methods, ratios computed within ci.method).

T19 threshold is 0.90 (not 0.95) because at N = 40 IID the parametric
pseudo-treated bootstrap targets the conditional variance V_t alone and
misses the finite-sample bias variance Var_{Λ,F}[b_t] in the absence of
factor structure to absorb it; the empirical SD across MC reps exceeds
the bootstrap SE by ~9%, which translates to ~0.91 coverage analytically
(see "Why fixed treated block" below for the law-of-total-variance
derivation).  At realistic factor-model DGPs (the minimal suite's
Scenario A / B), this gap shrinks and coverage returns to nominal.

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
under-cover by exactly that amount and the under-coverage is not a
property of the bootstrap procedure.  See gsynth-note section 2 for the
full derivation; this design choice mirrors gsynth-note's own MC
framework (footnote on Table 1, and section A.2).

## .Rbuildignore

This directory is excluded from the CRAN tarball via `.Rbuildignore`:

```
^tests/coverage-study$
```

Files here are tracked in git but are not part of the package
distribution.
