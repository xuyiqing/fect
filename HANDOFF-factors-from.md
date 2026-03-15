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

## Context for new conversation

> I'm working on the fect R package (`~/GitHub/fect`, branch `cfe`). We added a `factors.from` parameter to control whether factor estimation uses not-yet-treated or never-treated units. The implementation is committed (commits `be7bbcb` for tests, `53b461d` for code) but untested. Changes are in `R/default.R` (all 5 phases) and `R/boot.R` (parametric bootstrap dispatch). Read `~/GitHub/fect/HANDOFF-factors-from.md` for full context. First priority: fix the missing "both" guard in default.R line ~474 (`method == "mc"` should be `method %in% c("mc", "both")`), then run `devtools::test_file("tests/testthat/test-factors-from-refactor.R")` and fix any failures.
