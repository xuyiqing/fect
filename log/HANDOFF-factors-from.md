# Handoff: factors.from Refactoring (REQ-factors-from-001)

## Date
2026-03-15

## Status
REVIEW_PASSED (Skeptic verdict: PASS WITH NOTE)

---

## What was done

Added a `factors.from` parameter to the `fect` R package that controls whether factor estimation uses all not-yet-treated observations (`"notyettreated"`, the default) or only never-treated control units (`"nevertreated"`). This is a 5-phase refactoring:

| Phase | What | Status |
|-------|------|--------|
| 1 | `factors.from` parameter in all 3 signatures + II matrix filtering + output storage | Code committed, untested |
| 2 | `method="gsynth"` auto-sets `factors.from="nevertreated"` | Code committed, untested |
| 3 | Parametric bootstrap unlocked for ife/fe/cfe (was gsynth-only) | Code committed, untested |
| 4 | Deprecation warnings for `method="polynomial"` and `method="cfe_old"` | Code committed, untested |
| 5 | Input validation guards (invalid values, mc+nevertreated, insufficient units) | Code committed, untested |

---

## Target repository
- **Repo**: xuyiqing/fect
- **Branch**: `cfe`
- **Local path**: `~/GitHub/fect` (user's Mac)
- **Key commits**:
  - `be7bbcb` — 22 acceptance tests in `tests/testthat/test-factors-from-refactor.R`
  - `53b461d` — Implementation in `R/default.R` + `R/boot.R`

---

## Files modified

### R/default.R (all 5 phases)
| Location | Change |
|----------|--------|
| Line ~38 (fect signature) | Added `factors.from = "notyettreated"` |
| Line ~110 (fect.formula signature) | Added `factors.from = "notyettreated"` |
| Line ~214 (fect.formula passthrough) | Added `factors.from = factors.from` |
| Line ~288 (fect.default signature) | Added `factors.from = "notyettreated"` |
| Line ~469 (validation) | factors.from value check |
| Line ~474 (validation) | mc + nevertreated guard (**BUG: missing "both" — see open items**) |
| Line ~480 (gsynth setup) | Auto-set factors.from="nevertreated" for gsynth |
| Line ~384 (parametric block) | Relaxed: only blocks mc/both (was all non-gsynth) |
| Line ~1574 (II construction) | Nevertreated filtering with cummax + count guards |
| Line ~2262 (dispatch) | Deprecation warnings for polynomial/cfe_old |
| Line ~2893 (output list) | Added `factors.from = factors.from` |

### R/boot.R (Phase 3 only)
| Location | Change |
|----------|--------|
| Line ~702 | Parametric condition expanded: `c("gsynth")` → `c("gsynth", "ife", "cfe")` |
| Lines ~761-831 | Method-aware `draw.error()` dispatch: gsynth/ife/cfe |
| Line ~856 | `.export` list: added `"fect_fe"`, `"fect_cfe"` |

### tests/testthat/test-factors-from-refactor.R
22 acceptance tests across all 5 phases.

---

## Architecture: how it works

The **II matrix** (estimation sample indicator, TT x N) is the single control point. It's constructed at line ~1573 of default.R:

```r
II <- I
II[which(D == 1)] <- 0  ## treated observations excluded

if (factors.from == "nevertreated") {
    D.cum <- apply(D, 2, function(vec) { cummax(vec) })
    D.unit.sum <- colSums(D.cum)
    co.never <- which(D.unit.sum == 0)
    # ... count guards ...
    ever.treated <- setdiff(1:ncol(II), co.never)
    II[, ever.treated] <- 0  ## zero out all non-never-treated columns
}
```

ALL downstream estimators (fect_fe, fect_cfe, fect_gsynth, fect_mc, fect_cv, fect_boot) receive II as an argument. No downstream changes needed.

---

## Open items (must fix before merge)

### 1. Missing "both" in nevertreated guard (skeptic Note 2)
**File**: R/default.R, line ~474
**Current**: `if (method == "mc" && factors.from == "nevertreated")`
**Should be**: `if (method %in% c("mc", "both") && factors.from == "nevertreated")`
**Also update error message** to mention "both".

### 2. Runtime validation (skeptic Note 1)
No R was available during development. Must run:
```r
devtools::load_all("~/GitHub/fect")
devtools::test_file("~/GitHub/fect/tests/testthat/test-factors-from-refactor.R")
devtools::test("~/GitHub/fect")  # full regression suite
```

### 3. CFE parametric bootstrap with extra FE (skeptic concern)
In boot.R `draw.error()`, the `X.extra.FE`, `X.Z`, `X.Q` variables are passed to `fect_cfe` without being subset to the pseudo-sample indices. This may cause dimension mismatches when CFE is used with extra fixed effects + parametric bootstrap. Test with:
```r
fect(Y ~ D, data=df, method="cfe", r=0, se=TRUE, vartype="parametric",
     nboots=30, Z="some_var", parallel=FALSE)
```

---

## Nice-to-have (non-blocking)

- Add test for B3: r=0 invariance (factors.from should have no effect when r=0)
- Add test for H2: default output records "notyettreated"
- Verify parallel=TRUE works with the new .export list
- Move deprecation warnings from dispatch block to early validation (Phase 4)

---

## Full specification

### Phase 1: factors.from parameter
- Added to `fect()`, `fect.formula()`, `fect.default()` signatures after `force`, before `r`
- Default: `"notyettreated"` (preserves existing behavior)
- II matrix filtering at line ~1574: `cummax(D)` classifies units as ever-treated vs never-treated (handles reversals). `co.never <- which(D.unit.sum == 0)` identifies never-treated units. Guard: error if 0 never-treated units or fewer than `r + 1`. `II[, ever.treated] <- 0` zeroes out all non-never-treated columns.
- Stored in output: `output$factors.from <- factors.from`

### Phase 2: gsynth merger
- At line ~480: if `method=="gsynth"` and `factors.from=="notyettreated"`, auto-set to `"nevertreated"`
- Backward compatible because gsynth already restricted to never-treated internally via `co <- which(D.sum == 0)` in fect_gsynth.R
- No changes to fect_gsynth.R needed

### Phase 3: parametric bootstrap
- default.R line ~384: relaxed hard block from all-non-gsynth to only mc/both
- boot.R line ~702: condition expanded from `c("gsynth")` to `c("gsynth", "ife", "cfe")`
- boot.R `draw.error()` (lines ~761-831): method-aware dispatch
  - gsynth → `fect_gsynth(r = out$r.cv, ...)`
  - ife → `fect_fe(r.cv = out$r.cv, ...)` (note: `r.cv` not `r`)
  - cfe → `fect_cfe(r.cv = out$r.cv, X.extra.FE=..., X.Z=..., X.Q=..., ...)`
- boot.R line ~856: `.export` list updated for foreach parallel

### Phase 4: polynomial deprecation
- Warnings fire inside dispatch block at line ~2262, just before `fect_polynomial()` is called
- Both `method="polynomial"` and `method="cfe_old"` emit warnings suggesting `method="cfe"`
- Code paths still work — warnings only

### Phase 5: validation guards
- Line ~469: `factors.from` must be `"notyettreated"` or `"nevertreated"`
- Line ~474: mc + nevertreated → error (**BUG: missing "both"**)
- Line ~1580: 0 never-treated units → error
- Line ~1584: fewer than `r + 1` never-treated units → error

---

## Design: gsynth as a predictive routine, not an estimator

### Core Insight

**gsynth is NOT an estimator. It is a predictive routine that restricts factor estimation to never-treated units only.** The actual estimator is ife, cfe, etc. gsynth defines *which observations* are used for factor estimation (never-treated units) and then projects counterfactuals onto treated units.

### Two Orthogonal Axes

| Axis | Options | Controls |
|------|---------|----------|
| **Estimator** | fe (r=0), ife (r>0), cfe, mc | The statistical model / solver |
| **Predictive routine** | notyettreated (default), nevertreated | Which observations estimate factors |

`em` is a third axis determined by feasibility:

- `em=TRUE` (EM via `inter_fe_ub`) — needed when estimation sample has missing data
- `em=FALSE` (direct SVD via `inter_fe`) — only valid when estimation sample is complete
- `notyettreated` always has missing data by construction → `em=FALSE` is invalid
- `nevertreated` + balanced controls → `em=TRUE` reduces to `em=FALSE` (no missing data to fill)
- `nevertreated` + unbalanced controls → `em=TRUE` is required

### Why the current code diverges

`fect_fe` (ife path) passes the full TT x N matrix to `inter_fe_ub` with treated columns zeroed in II. The C++ solver estimates unit FEs for all N units, including treated units with zero observations — poorly identified.

`fect_gsynth` (nevertreated path) subsets to control units first (`Y.co`, `II.co`), estimates on Nco columns only, then projects counterfactuals onto treated units via pre-treatment residuals. This is the correct approach.

Empirical divergence: ~0.04–0.17 ATT gap on N=50–500 between the two paths with the same data. The gap does not shrink with N — it is structural, not numerical noise.

### Implementation Plan

**Phase 1 (current)**: Reroute + `em` parameter

1. Add `em` parameter to `fect()` (default TRUE for backward compat)
2. Validation: error if `em=FALSE` + `factors.from="notyettreated"`
3. Validation: warn + override if `em=FALSE` + nevertreated + unbalanced controls
4. When `factors.from="nevertreated"`: route ife to `fect_gsynth` (the nevertreated function)
5. `method="gsynth"` auto-sets `factors.from="nevertreated"` (already done)
6. Acceptance test: original gsynth vs rerouted ife+nevertreated must be numerically identical

**Phase 2 (later)**: Rename `fect_gsynth` → `fect_nevertreated` (pure refactoring, separate commit)

**Phase 3 (future)**: Extend `fect_nevertreated` to accept cfe as the estimator, completing the separation: predictive routine (nevertreated) × estimator (ife/cfe)

### Comparison Tests

| # | Comparison | Expected |
|---|-----------|----------|
| A | fe em=TRUE vs fe+nevertreated em=FALSE | Different (different obs) |
| B | ife(r>0) em=TRUE vs ife+nevertreated em=FALSE | Different (different factor space) |
| C | cfe em=TRUE vs cfe+nevertreated em=FALSE | Future work |
| **D** | **original gsynth vs rerouted ife+nevertreated+em=FALSE** | **Must be identical** |

---

## Roadmap

| Order | What | Tests affected | Status |
|-------|------|----------------|--------|
| **Phase 1** | Reroute ife+nevertreated → `fect_gsynth`, add `em` param | Fixes 9 nevertreated failures (1b,1c,1e,1f,1g,2a,2b,2c,5e) + new validation/comparison tests | Next |
| **Bug fix: boot.R** | Fix `draw.error()` dispatch for ife/cfe parametric bootstrap | Fixes 3a, 3b, 3c | After Phase 1 |
| **Bug fix: cfe_old** | Add deprecation warning in cfe_old dispatch path | Fixes 4b | After Phase 1 |
| **Phase 2** | Rename `fect_gsynth` → `fect_nevertreated` | 0 (pure refactoring) | After bug fixes |
| **Phase 3** | cfe + nevertreated support in `fect_nevertreated` | New tests | Future |

### Phase 1 Steps (detail)

1. Write tests first (test-driven):
   - **Test D**: `method="gsynth"` vs `method="ife", factors.from="nevertreated"` — must produce identical `att.avg`
   - **Validation**: `em=FALSE` + `factors.from="notyettreated"` → error
   - **Validation**: `method="gsynth"` auto-sets `em=FALSE`; `method="ife"` defaults `em=TRUE`
   - **Tests A/B**: informational comparison, not equality — both close to true tau
2. Add `em = TRUE` to all 3 `fect()` signatures + passthrough
3. Add validation guards for `em` + `factors.from` compatibility
4. `method="gsynth"` auto-sets `em = FALSE`
5. Reroute: when `method="ife"` + `factors.from="nevertreated"`, call `fect_gsynth()` instead of `fect_fe()`
6. Run test D — must be exact match
7. Run full test suite — expect 9 nevertreated tests to pass, 0 regressions

### Bug fix steps (after Phase 1)

**boot.R (3a, 3b, 3c)**: `draw.error()` passes wrong args to `fect_fe`/`fect_cfe`. `X.extra.FE`, `X.Z`, `X.Q` not subset to pseudo-sample indices. Fix: subset these matrices before calling the estimator.

**cfe_old (4b)**: Deprecation warning not emitted because `method="cfe_old"` is silently redirected before hitting the warning. Fix: move warning earlier or add it in the redirect path.

---

## Context for new conversation

> I'm working on the fect R package (`~/GitHub/fect`, branch `cfe`). We added a `factors.from` parameter to control whether factor estimation uses not-yet-treated or never-treated units. The implementation is committed (commits `be7bbcb` for tests, `53b461d` for code) but untested. Changes are in `R/default.R` (all 5 phases) and `R/boot.R` (parametric bootstrap dispatch). Read `~/GitHub/fect/HANDOFF-factors-from.md` for full context. First priority: fix the missing "both" guard in default.R line ~474 (`method == "mc"` should be `method %in% c("mc", "both")`), then run `devtools::test_file("tests/testthat/test-factors-from-refactor.R")` and fix any failures.
