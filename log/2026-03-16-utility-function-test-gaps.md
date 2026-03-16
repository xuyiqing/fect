# 2026-03-16 — Utility Function Test Coverage

> Run: `REQ-test-gaps-001` | Profile: r-package | Verdict: PASS

## What Changed

Added 53 new tests covering four previously undertested utility function areas in the fect package: `fect_mspe()`, `fect_mspe_sim()`, `esplot()`, and the `att.cumu()`/`effect()` relationship. These functions had minimal or zero test coverage. The new tests validate output structure, parameter variations, edge cases, and input validation. A pre-existing bug in `fect_mspe_sim()` was discovered (`.as_mask()` flattens matrix dimensions) and documented. No source code was modified — this is a test-only change.

## Files Changed

| File | Action | Description |
| --- | --- | --- |
| `tests/testthat/test-utility-functions.R` | created | 53 expectations across 4 test_that blocks covering fect_mspe, fect_mspe_sim, esplot, att.cumu/effect |

## Process Record

This section captures the full workflow history: what was proposed, what was tested, what problems arose, and how they were resolved.

### Proposal (from theorist)

**Implementation spec summary** (from `spec.md`):
- Create a single new test file `tests/testthat/test-utility-functions.R` with file-scope shared fixtures (fect model fits computed once) and 4 `test_that` blocks
- Test `fect_mspe()`: output structure, RMSE positivity, seed reproducibility (tolerance 1e-10), multi-model lists, multi-rep, custom hide_mask
- Test `fect_mspe_sim()`: input validation (Y0 requirement, hide_mask requirement)
- Test `esplot()`: 9 parameter variations (connected, SE-derived CI, show.count, highlight, fill.gap, start0, only.pre, only.post, from fect output)
- Test `att.cumu()` and `effect()`: structure, period variation, cumulative/non-cumulative, unit subset, relationship between both functions
- Construct Y0 = Y - eff * D for fect_mspe_sim tests since simdata lacks Y0
- Create no-reversal subset for `effect()` tests since simdata has treatment reversals

**Test spec summary** (from `test-spec.md`):
- 19 test scenarios defined plus 3 edge cases
- Tolerance for RMSE reproducibility: < 1e-10
- Structural checks: type, class, dimensions, column names, finiteness, positivity
- Parameter variation checks for esplot: each variation produces a valid ggplot

### Implementation Notes (from builder)

- File-scope fixtures: `out_base` (se=FALSE), `out_boot` (se=TRUE, nboots=20, keep.sims=TRUE), `out_norev` (no-reversal subset), `out_y0` (with Y0 column) — computed once, shared across all test blocks
- `simdata_norev` created by filtering to 101 units with monotonically non-decreasing D, avoiding `effect()` refusal to compute cumulative effects with reversals
- `fect_mspe_sim` positive-path test deferred due to `.as_mask()` bug: `as.logical(mask)` strips matrix dimensions, causing downstream `cbind(rr, cc)` indexing to fail
- Unit subset for `effect()` uses treated units from `out_norev` (filtered by `colSums(out_norev$D.dat) > 0`)
- 53 total expectations (vs spec's ~15-20 estimate) due to thorough structural validation

### Validation Results (from auditor)

- **Baseline** (before new tests): `devtools::test()` — FAIL 0 | WARN 256 | SKIP 0 | PASS 208
- **New tests only**: `devtools::test(filter="utility-functions")` — FAIL 0 | WARN 0 | SKIP 0 | PASS 53
- **Full regression**: `devtools::test()` — FAIL 0 | WARN 256 | SKIP 0 | PASS 261
- Net change: +53 tests, 0 regressions, 256 pre-existing warnings unchanged
- All 19 test scenarios from test-spec.md verified (scenario 5 not tested due to source bug, documented)
- Tolerance integrity: RMSE reproducibility uses tolerance = 1e-10 as specified

### Problems Encountered and Resolutions

No BLOCK, HOLD, or STOP signals occurred during this workflow.

The `.as_mask()` bug in `R/fect_mspe.R` (line 226) was discovered during implementation. This is a pre-existing source code bug — `as.logical(mask)` strips matrix dimensions, making `fect_mspe_sim()` non-functional for its primary use case. The builder correctly worked around it by testing only the input validation paths. This bug should be addressed in a separate fix run.

### Review Summary (from skeptic, if available)

Pending — skeptic review follows scribe.

- **Pipeline isolation**: pending
- **Convergence**: pending
- **Tolerance integrity**: pending
- **Verdict**: pending

## Design Decisions

1. **File-scope fixtures over per-block fitting**: Each fect fit takes ~30 seconds. Computing fixtures once at file scope and sharing across all 4 test blocks reduces total test runtime from ~8 minutes to ~2 minutes. This is safe because tests are read-only on the fect output objects.

2. **No-reversal subset for effect() tests**: `effect()` with `cumu=TRUE` refuses to compute when treatment reversals exist (warns and returns NULL). Rather than mock or modify the function, a clean subset `simdata_norev` (101 units) was created by filtering to units with monotonically non-decreasing D. This tests `effect()` on its intended input.

3. **Deferred fect_mspe_sim positive-path test**: The `.as_mask()` bug is a source code issue, not a test design issue. Testing the input validation paths (Y0 check, hide_mask check) verifies the function's defensive behavior. The positive-path test should be added after the bug is fixed.

4. **Structural assertions over value assertions**: Most tests check output structure (type, class, dimensions, column names, finiteness) rather than exact numerical values. This is appropriate for utility functions where the exact values depend on the estimation procedure, but the structure must be stable.

## Handoff Notes

- **fect_mspe_sim is broken**: The `.as_mask()` function at `R/fect_mspe.R:226` uses `as.logical(mask)` which strips matrix dimensions. Fix: use `matrix(as.logical(mask), nrow=nrow(mask), ncol=ncol(mask))` or similar. After fixing, add positive-path tests for `fect_mspe_sim()`.
- **Test runtime**: The utility function test file takes ~2 minutes due to 4 fect fits. This is unavoidable since the functions under test operate on fect output objects.
- **Treatment reversals in simdata**: 99 of 200 units have treatment reversals. Any test of `effect(cumu=TRUE)` must use a no-reversal subset. The pattern in the test file (filter by monotonically non-decreasing D) can be reused.
- **Y0 construction**: simdata has `eff` (treatment effect) and `D` (treatment indicator) columns. Y0 = Y - eff * D gives the true counterfactual outcome for simulation-based validation.
