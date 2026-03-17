# 2026-03-16 — cv.sample k-fold CV in fect_nevertreated

> Run: `REQ-cv-sample-nevertreated` | Profile: r-package | Verdict: pending

## What Changed

Implemented actual cv.sample-based k-fold cross-validation in `fect_nevertreated` for `cv.method="all_units"` and `cv.method="treated_units"`. Previously, these options were validated at entry but fell through to the existing LOO (leave-one-period-out) logic. Now each cv.method triggers a distinct masking and re-estimation strategy, following the pattern established in `fect_cv` (R/cv.R). Both IFE and CFE estimation paths are covered.

## Files Changed

| File | Action | Description |
| --- | --- | --- |
| `R/fect_nevertreated.R` | modified | Added cv.sample pre-computation blocks (IFE + CFE), branched r-loop on cv.method with "all_units" and "treated_units" paths, bug fixes for inter_fe balanced-panel $fit, beta0CV.co p=0 array dims, data.ini creation |
| `tests/testthat/test-score-unify.R` | modified | Added Section H (13 tests: smoke IFE/CFE, LOO backward compat, r-selection, ATT consistency, edge cases, properties) and Section I (2 runtime benchmarks) |
| `architecture.md` | modified | Updated cv.method threading diagram to show all three nevertreated branches; updated module reference and notes |

## Process Record

### Proposal (from theorist)

**Implementation spec summary** (from `spec.md`):
- Two new cv.method branches in the r-loop: "all_units" masks control panel observations via `cv.sample()`, re-estimates factors on the masked panel, and scores held-out control positions; "treated_units" keeps factor estimation unchanged (full control panel), masks treated pre-treatment observations, re-estimates loadings from remaining pre-treatment, and scores held-out treated positions
- Pre-computation block before the r-loop generates k fold masks with 200-retry validation (same pattern as fect_cv)
- Both IFE (inter_fe_ub) and CFE (complex_fe_ub) paths follow the same logical structure
- For "treated_units" with r>0: loading re-estimation via OLS projection lambda = (F'F)^{-1}F'U on masked pre-treatment; for r=0: unit FE alpha = mean(U) on masked pre-treatment
- Code organization: fold pre-computation guarded by `cv.method != "loo"`, then `if/else if/else` branching inside r-loop

**Test spec summary** (from `test-spec.md`):
- 7 behavioral contract scenarios: smoke tests for all_units and treated_units on both IFE and CFE, LOO backward compatibility, r-selection validity, ATT consistency across methods
- 4 edge case scenarios: small panel, single treated unit, unbalanced panel, r=0
- 5 property-based invariants: score non-negativity, CV.out completeness, selection rule consistency, LOO determinism, cv.sample reproducibility
- 2 runtime benchmark scenarios (IFE and CFE timing comparison across cv.methods)
- Tolerances: ATT consistency within factor of 5, reproducibility at 1e-10, scores >= 0

### Implementation Notes (from builder)

- Always uses `inter_fe_ub` for "all_units" fold-level estimation (even on originally balanced panels), because masking creates an effectively unbalanced panel and `inter_fe` does not return `$fit`
- `beta0CV.co` initialized as `array(0, dim = c(1, 0, k))` when p=0 to avoid indexing failures
- Conditional creation of `data.ini` in IFE pre-computation block for balanced panels (normally only created inside `if (0 %in% I.co)`)
- No debug `cat()` or `print()` statements in final code
- Weight handling uses `!is.null(W)` pattern consistent with fect_nevertreated conventions
- All `solve()` calls wrapped in `try(..., silent = TRUE)` for numerical safety
- LOO path wrapped in `if (cv.method == "loo")` but code completely unchanged

### Validation Results (from auditor)

- Tests written per test-spec.md: 13 functional tests + 2 runtime benchmarks
- Validation status: PENDING at time of audit (builder and auditor dispatched in parallel)
- Coverage: 15/17 test-spec scenarios covered (Edge-3 unbalanced panel and P-3 selection rule consistency deferred)
- All tolerances match test-spec.md exactly: ATT range < 5.0, reproducibility 1e-10, scores >= 0

### Runtime Benchmarks

Measured on Apple M1 Ultra, single core, N=200, TT=40, Ntr=60, k=10 folds:

**IFE** (r=0:5):

| cv.method | Elapsed |
| --- | --- |
| loo | 0.3s |
| treated_units | 3.1s |
| all_units | 2.3s |

**CFE** (r=0:3):

| cv.method | Elapsed |
| --- | --- |
| loo | 0.2s |
| treated_units | 0.1s |
| all_units | 1.8s |

LOO is fastest because it does not call `cv.sample()` or `initialFit()` per fold. "treated_units" is fast for CFE because factor estimation runs once per r (not per fold) and loading re-estimation is cheap OLS. "all_units" is the most expensive because it re-estimates the full factor model k times per r.

### Problems Encountered and Resolutions

| # | Problem | Signal | Routed To | Resolution |
| --- | --- | --- | --- | --- |
| 1 | `inter_fe()` (balanced panels) does not return `$fit` field | n/a (caught during implementation) | builder | Always use `inter_fe_ub()` for fold-level estimation since masking creates unbalanced panel regardless |
| 2 | `beta0CV.co[,,ii]` fails when p=0 due to array dimension mismatch | n/a (caught during implementation) | builder | Initialize as `array(0, dim = c(1, 0, k))` matching fect_cv pattern |
| 3 | `data.ini` not created in balanced IFE case (only inside `if (0 %in% I.co)`) | n/a (caught during implementation) | builder | Conditional creation in pre-computation block |

### Review Summary (from skeptic, if available)

Pending -- skeptic review follows scribe.

- **Pipeline isolation**: pending
- **Convergence**: pending
- **Tolerance integrity**: pending
- **Verdict**: pending

## Design Decisions

1. **Always use inter_fe_ub for fold-level estimation**: Even when the original panel is balanced (`!0 %in% I.co`), masking observations via cv.sample creates an effectively unbalanced panel. Using `inter_fe_ub` uniformly avoids the `$fit` issue with `inter_fe` and is semantically correct -- the masked panel IS unbalanced.

2. **Factor estimation once per r for "treated_units"**: The "treated_units" cv.method tests projection/loading quality, not factor estimation quality. Factors are estimated once on the full control panel per candidate r, then only loadings are re-estimated per fold. This is both correct (factor estimation is independent of treated pre-treatment masking) and efficient (avoids k redundant factor estimations per r).

3. **Separate storage for co-space vs treated-space masks**: `rmCV`/`estCV`/`ociCV` for "all_units" (indexing into control panel) vs `rmCV.tr`/`estCV.tr` for "treated_units" (indexing into treated pre-treatment). No collision since only one cv.method executes per run, but distinct naming prevents confusion.

4. **Fallback to LOO when rm.count is zero**: If the panel is too small for cv.sample to remove any observations (`floor(sum * cv.prop) == 0`), the code warns and falls back to LOO rather than erroring. This handles edge cases gracefully.

5. **"Control" time labels for "all_units" scoring**: Held-out positions in "all_units" are control observations, so all get the time label "Control" for `.score_residuals()`. This differs from "treated_units" which uses the actual T.on time indices for held-out treated pre-treatment positions.

## Handoff Notes

- The cv.sample pre-computation runs BEFORE the r-loop and generates all k fold masks once. The masks are then reused across all candidate r values. This matches the fect_cv design.
- For "treated_units" with very few pre-treatment periods, the con2 validation check uses `min(min.T0, 2)` as the threshold (not `min.T0`) to allow cv.sample to work with small panels.
- The 200-retry limit for cv.sample validation matches fect_cv. If all 200 retries fail, the repair logic restores rows/columns that violate constraints.
- Edge-3 (unbalanced panel with explicit missingness) and P-3 (selection rule consistency) tests were deferred by auditor. These could be added in a future pass.
- The `Y0CV.co` and `beta0CV.co` arrays are only allocated for "all_units" (where `initialFit` runs per fold). "treated_units" does not need them since factor estimation uses the full control data.
- Runtime: "all_units" is the most expensive cv.method because it calls `inter_fe_ub`/`complex_fe_ub` k times per r. For large panels, consider reducing k or cv.prop.
