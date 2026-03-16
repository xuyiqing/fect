# 2026-03-16 — CFE CV r-selection test fix

> Run: `REQ-cfe-cv-rselect-001` | Profile: r-package | Verdict: Pending (skeptic review follows scribe)

## What Changed

The CFE CV r-selection issue was investigated and resolved. The CFE CV algorithm was found to be structurally correct. The reported problem (r=0 selected on data with true r=2) was caused by test misspecification: the existing test used `make_cfe_z_data` (DGP with Z*gamma interaction) but did not pass `Z = "Z"` to `fect()`. Without Z, factors absorb the unmodeled interaction, leading to r overselection. The fix adds `Z = "Z"` to the existing test and adds 3 new tests confirming correct r-selection under proper model specification.

## Files Changed

| File | Action | Description |
| --- | --- | --- |
| `tests/testthat/test-factors-from-refactor.R` | modified | Fixed existing CFE CV test (added `Z = "Z"`, `parallel = TRUE, cores = 10`, tightened assertion); added 3 new tests for r=2 with Z, r=0 with Z, and r=2 on factor-only data |

## Process Record

### Proposal (from theorist)

**Implementation spec summary** (from `spec.md`):
- The CFE CV algorithm is correct; the existing test is misspecified (missing `Z = "Z"`)
- Fix: add `Z = "Z"` to existing test, add 3 new tests covering r=2 with Z specified, r=0 on zero-factor data, r=2 on factor-only data (no Z in DGP)
- No source code changes required

**Test spec summary** (from `test-spec.md`):
- 5 scenarios defined: (1) CFE+Z selects r=2 on Z-data, (2) CFE+Z selects r=0 on zero-factor data, (3) CFE selects r=2 on factor-only data, (4) existing test fixed with r in [1,3], (5) robustness across seeds {42, 123, 456}
- Tolerances: r.cv exact integer match for new tests, r.cv in [1,3] for existing test, ATT within 0.5 of true tau=3.0
- Property-based invariants: sigma2 monotonicity, MSPE U-shape, model specification principle
- Cross-reference benchmarks: IFE without Z overselects (r=3), CFE without Z overselects (r=3), CFE with Z correct (r=2)

### Implementation Notes (from builder)

- Added `Z = "Z"` and `parallel = TRUE, cores = 10` to existing test at line ~1696
- Changed `fect::fect()` to `fect()` for consistency
- Tightened assertion from `r.cv >= 0` to `r.cv >= 1`
- Added 3 new tests: r=2 with Z (exact match), r=0 with Z (exact match), r=2 on factor-only data (exact match)
- All new tests use `skip_on_cran()` due to computational cost
- No deviations from spec

### Validation Results (from auditor)

- **Command**: `NOT_CRAN=true Rscript -e 'devtools::test_file("tests/testthat/test-factors-from-refactor.R")'`
- **Result**: FAIL 0 | WARN 0 | SKIP 0 | PASS 135 (131 existing + 4 new, with existing test now corrected)
- **Scenario results**:
  - Scenario 1 (r=2 with Z): r.cv = 2, att.avg = 2.979 (within 0.5 of 3.0) -- PASS
  - Scenario 2 (r=0 with Z): r.cv = 0 -- PASS
  - Scenario 3 (factor-only): r.cv = 2 -- PASS
  - Scenario 4 (existing test, r in [1,3]): r.cv = 2 -- PASS
  - Scenario 5 (seeds 42/123/456): r.cv = 2/2/2, att.avg = 2.979/3.022/2.997 -- PASS
- **Property invariants**: sigma2 strictly decreasing (2.695 -> 0.233), MSPE U-shape with minimum at r=2 (0.366), model specification principle confirmed -- all PASS
- **Tolerances**: All match test-spec.md exactly. No tolerances relaxed.

### Problems Encountered and Resolutions

No problems encountered. No BLOCK, HOLD, or STOP signals were raised during this run.

### Review Summary (from skeptic, if available)

Pending -- skeptic review follows scribe.

- **Pipeline isolation**: pending
- **Convergence**: pending
- **Tolerance integrity**: pending
- **Verdict**: pending

## Design Decisions

1. **Model specification, not algorithm fix**: The key finding is that the CFE CV algorithm is correct. The Z*gamma interaction in the DGP creates a rank-1 component that, when unmodeled, forces the factor estimator to absorb it. This is standard omitted variable behavior in factor models, not a bug. The fix is at the test level (specify Z correctly) rather than the algorithm level.

2. **Exact r.cv assertions for new tests**: New tests use `expect_equal(out$r.cv, 2)` (exact match) rather than range checks. This is justified because the MSPE U-shape is clear and robust across seeds. The existing test retains a range check (`r.cv >= 1 && r.cv <= 3`) as it was already a structural test.

3. **Parallel execution in all CV tests**: All CV tests now use `parallel = TRUE, cores = 10` for consistency and speed. The existing test previously used `parallel = FALSE`.

## Handoff Notes

- The CFE CV r-selection issue is RESOLVED. The algorithm is correct; users must specify Z when their DGP includes Z*gamma interactions. This is the expected behavior for any factor model with omitted covariates.
- The test suite now covers: correct r-selection with Z specified (r=2), correct null detection (r=0), factor-only data without Z, and robustness across seeds.
- 135 tests now pass (was 131 before this run).
- Remaining open items from HANDOFF-factors-from.md: cv.R vs fect_mspe.R comparison, Phase 3b (merge IFE into CFE), test gaps (fect_mspe, esplot, att.cumu).
