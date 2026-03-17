# Cross-Validation and Model Comparison Functions

This table documents the three CV/scoring functions in fect after the scoring unification refactor. All three share `.score_residuals()` for consistent scoring.

## Comparison Table

| | `fect_cv` (notyettreated) | `fect_nevertreated` CV (nevertreated) | `fect_mspe` (post-hoc) |
|---|---|---|---|
| **Used by** | IFE, MC, CFE with `factors.from="notyettreated"` | gsynth, IFE+nevertreated, CFE+nevertreated | User — any finished model(s) |
| **Purpose** | Automatic hyperparameter selection (r or λ) | Automatic hyperparameter selection (r) | Manual model comparison / validation |
| **Grid search** | r (IFE/CFE) or λ (MC) | r | None |
| | | | |
| **`cv.method` options** | `"all_units"`, `"treated_units"` | `"treated_units"`, `"all_units"`, `"loo"` | `"all_units"`, `"treated_units"` |
| **`cv.method` default** | `"all_units"` | `"treated_units"` | `"all_units"` |
| **`"all_units"`** | cv.sample from all units' control obs | cv.sample from never-treated (control) panel | cv.sample from all units' control obs |
| **`"treated_units"`** | cv.sample from eventually-treated units' pre-treatment obs | cv.sample from eventually-treated units' pre-treatment obs | cv.sample from eventually-treated units' pre-treatment obs |
| **`"loo"`** | Not available | LOO on treated pre-treatment, one period at a time | Not available |
| | | | |
| **cv.sample sub-options** | | *(apply to `"all_units"` and `"treated_units"`)* | |
| `k` | 5 (number of random hold-out rounds) | 5 | 5 |
| `cv.prop` | 0.1 (proportion of obs held out per round) | 0.1 | 0.1 |
| `cv.nobs` | 3 (block size for structured removal) | 3 | 3 |
| `cv.donut` | 1 (periods excluded around treatment in evaluation) | 1 | 1 |
| `min.T0` | 5 (minimum pre-treatment obs per unit) | 5 | 5 |
| | | | |
| **Scoring** | `.score_residuals()` — all 9 criteria | `.score_residuals()` — all 9 criteria | `.score_residuals()` — all 9 criteria |
| **`criterion`** | Yes (default `"mspe"`) | Yes (default `"mspe"`) | Yes (default `"mspe"`) |
| **Selection rule** | >1% improvement | >1% improvement | No selection — reports scores |
| **Observation weights (W)** | Yes | Yes | Yes |
| **Period weights (count.T.cv)** | Yes | Yes | Yes |
| **norm.para** | Yes | Yes | Yes |
| **Carryover handling** | Via `carryover.rm` upstream | Via `carryover.rm` upstream | Via `carryover.rm` upstream |
| | | | |
| **Multi-model comparison** | No | No | Yes (list of fect outputs) |

## Why `"loo"` is exclusive to nevertreated

Leave-one-period-out holds out treated units' pre-treatment obs at time t, then re-estimates and scores the prediction. In the nevertreated setting, estimation uses control units only — controls retain all periods, so estimation is never compromised.

In the notyettreated setting, all units contribute to estimation. Holding out treated obs at period t is safe (controls retain it for time FE identification), but `"treated_units"` with cv.sample achieves the same goal more robustly — especially when pre-treatment periods are few, where LOO gives too few data points for stable selection.

## Scoring Criteria

All three functions use `.score_residuals()` which computes:

| Criterion | Formula | Description |
|---|---|---|
| MSPE | `SSE / n` | Mean squared prediction error |
| WMSPE | `Σ(w_period × e²) / n` | Period-weighted MSPE |
| GMSPE | `exp(Σ log(e²) / n)` | Geometric MSPE |
| WGMSPE | `exp(Σ log(w_period × e²) / n)` | Weighted geometric MSPE |
| MAD | `median(|e² - median(e²)|)` | Median absolute deviation of squared errors |
| Moment | `Σ(w_period × |mean(e by period)|) / Σ(w_period)` | Weighted moment condition |
| GMoment | `Σ(w_period × gmean(|e| by period)) / Σ(w_period)` | Geometric moment condition |
| RMSE | `sqrt(MSPE)` | Root mean squared prediction error |
| Bias | `mean(e)` | Mean prediction error |

When observation weights `W` are provided, `e²` is replaced by `W × e²` and denominators use `Σ(W)` instead of `n`.

## Masking Strategies

| Strategy | Description | Available in |
|---|---|---|
| **`"all_units"`** | k-fold cv.sample from all units' control obs. Structured removal with donut, min pre-treatment obs constraint. Tests factor estimation quality on the full panel. | `fect_cv` (default), `fect_nevertreated` (option), `fect_mspe` (default) |
| **`"treated_units"`** | k-fold cv.sample from eventually-treated units' pre-treatment obs only. Same structured removal. Tests counterfactual prediction quality on target units. | `fect_cv` (option), `fect_nevertreated` (default), `fect_mspe` (option) |
| **`"loo"`** | Leave-one-time-period-out; hold out treated pre-treatment obs at one period, re-estimate, score prediction. Established method from gsynth. Unstable when pre-treatment periods are few. | `fect_nevertreated` (option) |

## Key Architectural Difference

- **`fect_cv` (notyettreated)**: Masks observations in the **full panel**, re-estimates with the masked II indicator. Both treated and control units contribute to estimation.

- **`fect_nevertreated` CV (nevertreated)**: Estimates factors on **control units only**, then projects onto treated units. `"treated_units"` tests projection quality by holding out treated pre-treatment obs (affects loading estimation λ_tr). `"all_units"` tests factor estimation quality by holding out control obs. `"loo"` is the legacy gsynth method.

- **`fect_mspe` (post-hoc)**: Works on finished `fect()` output objects. Re-runs `fect()` with held-out observations set to NA. Compares multiple models side-by-side using the same cv.sample masking and scoring as the internal CV functions.
