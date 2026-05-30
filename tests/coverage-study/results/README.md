# Canonical coverage results

CSV summaries from `run_minimal_coverage.R` (and the conditional
`run_minimal_coverage_tail_rerun.R` follow-up).  Saved here so `07-inference.Rmd`
in the Quarto book can read them directly without depending on `/tmp/`.

## Files

- `minimal_summary_K200_nb200.csv` --- 4 scenarios × {5 ci.methods for A/B/C1, 1 for C2}
  at `K = 200`, `nboots = 200`.  Run dated 2026-05-03 20:04, fect commit
  `1ac0e0c` (state: pre-coverage-suite-commit; the run uses the same R/ source
  as the final v2.4.2 surface).
- `minimal_tail_rerun_summary_K200_nb1000.csv` --- C1 (TWFE bootstrap) only,
  re-run at `nboots = 1000` because basic / percentile / bc / bca tail-CI cells
  came in below 0.93 at `nboots = 200`.  Same MC seeds as the parent run.

## DGP / inference settings (recap)

| Scenario | DGP | Inference |
|---|---|---|
| A | factor `r=2`, IID errors. N_tr=5, N_co=50, T=30, T0=20, ATT_t = t for t = 1..10, D.sd=1, λ_i ~ U(-√3,√3), F_t ~ N(0,1), α_i ~ U(-√3,√3), ξ_t ~ N(0,1), μ=5; λ, F, α, ξ redrawn each rep | `vartype = "parametric"`, `para.error = "auto"` (-> empirical), all 5 ci.methods |
| B | same as A but errors AR(1) ρ = 0.8, marginal variance 1 | same as A; `auto` -> ar |
| C1 | additive TWFE r=0, AR(1) ρ = 0.5. N_tr=20, N_co=80, T=30, T0=20, ATT = 3 constant (D.sd=0); α, ξ redrawn each rep | `vartype = "bootstrap"` (cluster), all 5 ci.methods |
| C2 | same DGP as C1 | `vartype = "jackknife"`, `ci.method = "normal"` only (E&T 1993 §11) |

Estimator for A/B: `method = "ife", r = 2, force = "two-way",
time.component.from = "nevertreated", CV = FALSE`.
Estimator for C: `method = "fe", force = "two-way",
time.component.from = "notyettreated"`.

Coverage target for A/B: realized average treated-post effect within rep
(matches gsynth-note's `simulate-xu-{iid,ar1}-rfit2.R`).  Coverage target
for C: population ATT = 3.

K=200 reps → MC SE about 0.015 around coverage 0.95.
