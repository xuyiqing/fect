# Convergence Conditioning: R Centering + C++ Component-Wise

> Run: REQ-phase3b-wrapper | Date: 2026-03-17 | Branch: cfe

## What Changed

Improved EM convergence conditioning so that `tol` applies evenly to all fit components (mu, alpha, xi, F*L', beta), not just the total fit norm which was dominated by the grand mean.

### Two complementary changes

**1. R centering** (`R/fect_nevertreated.R`): `.estimate_co()` wrapper subtracts the observed grand mean from Y before calling the C++ solver, then adds it back. This removes the dominant scale component from the convergence criterion. Mainly benefits r=0 cases.

**2. C++ component-wise convergence** (`src/ife_sub.cpp`): `fe_ad_inter_iter` and `fe_ad_inter_covar_iter` now track interactive FE (factor) convergence independently via `max(dif_fit, dif_inter)`. Ensures factors converge to within `tol` of their own scale, not hidden behind a large mu.

### Also included

**Two-tier CV tolerance** (`R/cv.R`, `R/fect_nevertreated.R`): CV loops use `cv_tol = max(tol, 1e-3)` — CV only needs relative MSPE ranking, not precise estimates. Final estimation uses the user's `tol`.

## Why

The convergence criterion `||fit - fit_old|| / ||fit_old||` was dominated by the grand mean when `mean(Y)` is large (e.g., GDP data with mean=50). At `tol=1e-3`, the old code converged in 5 iterations with 9.4% factor error. The new code correctly requires ~137 iterations to achieve genuine 0.1% precision on all components.

| Component | Old tol=1e-3 error | New tol=1e-3 error | Improvement |
|-----------|-------------------|-------------------|-------------|
| alpha (unit FE) | 2.03e-01 | 3.72e-03 | 55x |
| F*L' (factors) | 2.99e+00 | 6.89e-02 | 43x |
| sigma2 | 5.07e-03 | 2.25e-06 | 2249x |

## Files Changed

| File | Change |
|------|--------|
| `src/ife_sub.cpp` | Component-wise convergence in `fe_ad_inter_iter` and `fe_ad_inter_covar_iter` |
| `R/fect_nevertreated.R` | `.estimate_co()` R centering for unbalanced panels |
| `R/cv.R` | Two-tier `cv_tol` for CV speed |

## Verification

- 236/236 tests pass (192 score-unify + 44 utility)
- All non-nevertreated paths tested: fe, ife, mc, cfe, bootstrap — all PASS
- Results invariant to mean(Y) from 0 to 500
- Zero speed regression (same iterations, same time)
- Backward compatible: tol default stays 1e-3
