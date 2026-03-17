# Handoff: factors.from Refactoring (REQ-factors-from-001)

## Date
2026-03-16 (updated)

## Status
All major phases COMPLETE and PUSHED to `cfe`. Score unification (Phases 1-2), cv.sample k-fold CV in fect_nevertreated (all three cv.method paths), and parallel CV with auto-threshold are all committed. 259/259 tests pass (215 score-unify including 15 parallel/CRAN-guarded + 44 utility). Phase 3b (merge IFE into CFE) and fect_mspe+CV=TRUE bug remain.

---

## What was done

Added `factors.from` and `em` parameters to the `fect` R package, rerouted `ife+nevertreated` through `fect_nevertreated` (renamed from `fect_gsynth`), removed obsolete methods (`polynomial`, `bspline`, `cfe_old`), implemented CFE bifurcation in `fect_nevertreated`, fixed CV routing, unified scoring and CV methods, implemented cv.sample k-fold CV with parallel support, and updated the Quarto user manual.

| Phase | What | Status |
|-------|------|--------|
| 1 | `factors.from` parameter in all 3 signatures + II matrix filtering + output storage | **Done** (`53b461d`, `cafb66d`) |
| 2 | `method="gsynth"` auto-sets `factors.from="nevertreated"` | **Done** (`37d516b`) |
| 3 | Parametric bootstrap unlocked for ife/fe/cfe (was gsynth-only) | **Done** (`070d574`) |
| 4 | Removed `method="polynomial"`, `"bspline"`, `"cfe_old"` entirely | **Done** (`953da8a`, `efa2af0`) |
| 5 | Input validation guards | **Done** (`cafb66d`) |
| 6 | `em` parameter + reroute ife+nevertreated → `fect_nevertreated` | **Done** (`37d516b`) |
| 7 | Rename `fect_gsynth` → `fect_nevertreated` | **Done** (`fa9baf2`) |
| 3a | CFE bifurcation in `fect_nevertreated` | **Done** (`c0cab63`) |
| — | Category I bootstrap tests (I1-I11) + est.avg bug fix | **Done** (`5136a03`) |
| — | r=0 invariance test (I12) | **Done** (`e5f098e`) |
| — | CV routing fix for cfe+nevertreated | **Done** (`4cfe25c`) |
| — | Quarto book: two research settings, chapter reorder, factors.from/em docs | **Done** (`219510d`, `0fa790c`, `fd2abd3`) |
| — | CFE CV r-selection: test misspecification fix | **Done** |
| — | Utility function tests: fect_mspe, esplot, att.cumu, effect | **Done** (`5460b13`) |
| — | Score unification Phase 1: `.score_residuals()` extraction | **Done** (`caa456f`) |
| — | Score unification Phase 2: `cv.method` unification + fect_mspe simplification | **Done** (`eac0187`) |
| — | cv.sample k-fold CV in fect_nevertreated (all_units + treated_units) | **Done** (`cea6966`) |
| — | Parallel CV with auto-threshold (Nco*TT > 20k for all_units) | **Done** (`cea6966`) |

---

## Target repository
- **Repo**: xuyiqing/fect
- **Branch**: `cfe`
- **Local path**: `~/GitHub/fect`
- **HEAD**: `cea6966`

### Commit history (oldest → newest)
| Commit | Description |
|--------|-------------|
| `be7bbcb` | Add acceptance tests for factors.from refactoring (5 phases) |
| `53b461d` | Add factors.from parameter and refactor estimator architecture |
| `cafb66d` | Fix factors.from='nevertreated' validation and runtime bugs |
| `070d574` | Fix parametric bootstrap for ife/cfe and cfe_old deprecation warning |
| `37d516b` | Reroute ife+nevertreated through fect_gsynth and add em parameter |
| `efa2af0` | Remove all method='cfe_old' code paths |
| `953da8a` | Remove all method='polynomial' and 'bspline' code paths |
| `25cd5f8` | Move handoff notes to log/ directory |
| `fa9baf2` | Rename fect_gsynth to fect_nevertreated |
| `c0cab63` | Phase 3a: CFE bifurcation in fect_nevertreated with BCD, 3-layer projection, boot.R fixes |
| `97818bd` | docs: update handoff with Phase 3a details and Category I test plan |
| `5136a03` | test: add Category I bootstrap tests (I9-I11) and fix est.avg field name |
| `8d69239` | docs: update handoff with Category I completion |
| `e5f098e` | test: add r=0 invariance test (I12) |
| `4cfe25c` | feat: fix CV routing for cfe+nevertreated via fect_cv delegation |
| `cdb51ae` | docs: update handoff with CV fix |
| `219510d` | docs: Quarto book — two research settings, factors.from, esplot, fect_mspe |
| `0fa790c` | docs: remove em=FALSE from gsynth equivalence (em is no-op for nevertreated) |
| `fd2abd3` | docs: reorder chapters, integrate two research settings, remove dead R script links |
| `5460b13` | test: add utility function tests for fect_mspe, esplot, att.cumu, effect |
| `caa456f` | refactor: unify CV scoring with shared .score_residuals() |
| `eac0187` | refactor: unify cv.method parameter across CV functions (Phase 2) |
| `cea6966` | feat: cv.sample k-fold CV + parallel support in fect_nevertreated |

---

## Test results

**259/259 tests pass** (215 in `test-score-unify.R` + 44 in `test-utility-functions.R`; 77 in `test-factors-from-refactor.R` are CRAN-guarded/skipped)

---

## CV architecture after Phase 2 + cv.sample

### cv.method parameter

| Function | Valid values | Default |
|---|---|---|
| `fect_cv` | `"all_units"`, `"treated_units"` | `"all_units"` |
| `fect_nevertreated` | `"treated_units"`, `"all_units"`, `"loo"` | `"treated_units"` |
| `fect_mspe` | `"all_units"`, `"treated_units"` | `"all_units"` |

### cv.sample k-fold CV in fect_nevertreated

| cv.method | Masks what | Re-estimates what | Tests what |
|---|---|---|---|
| `loo` | 1 treated pre-treatment period | Loadings only | Projection quality (legacy gsynth) |
| `all_units` | k-fold on control panel | Factors (full re-estimation via inter_fe_ub) | Factor estimation quality |
| `treated_units` | k-fold on treated pre-treatment | Loadings only | Projection quality (structured) |

### Parallel CV auto-threshold

Auto-enables for `all_units` when `Nco * TT > 20000` (default `parallel=TRUE` from fect()). `treated_units` and `loo` always sequential (per-fold cost too low). Uses `doFuture::registerDoFuture()` + `future::multisession`. Default 4 cores. `parallel=FALSE` forces sequential.

**Benchmarks** (k=10, 5 cores, Apple M1 Ultra, 20% missing):
- N=500, TT=60: seq=15.6s, par=10.3s, **1.5x**
- N=1000, TT=80: seq=38.9s, par=21.0s, **1.9x**
- N=2000, TT=80: seq=69.6s, par=35.2s, **2.0x**

---

## Open items

### fect_mspe + CV=TRUE bug
`fect_mspe()` throws "No valid residuals collected from cv.sample folds" when called on a model fitted with `CV=TRUE`. The cv.sample feasibility check fails because the re-fit changes the II matrix. Workaround: fit with `CV=FALSE` before passing to `fect_mspe()`.

### Phase 3b: Merge IFE into CFE
- Verify test E0: `complex_fe_ub` with empty CFE arrays ≡ `inter_fe_ub`
- Verify test E4: `ife+nevertreated` ≡ `cfe+nevertreated` (no extras)
- Only then: replace `inter_fe_ub` calls with `complex_fe_ub` in `fect_nevertreated`
- Re-run full test suite to confirm no regressions

---

## Context for new conversation

> I'm working on the fect R package (`~/GitHub/fect`, branch `cfe`). All major refactoring phases are complete and pushed. 259/259 tests pass.
>
> **Recent commits** (all pushed to `cfe`):
> - `caa456f` — Score unification Phase 1: `.score_residuals()` shared scoring
> - `eac0187` — Phase 2: `cv.method` unification, fect_mspe simplification, 1% rule
> - `cea6966` — cv.sample k-fold CV (all_units + treated_units) + parallel auto-threshold
>
> **Open tasks** (in priority order):
>
> 1. **fect_mspe + CV=TRUE bug** — fect_mspe errors on CV-fitted models ("No valid residuals"). Likely cv.sample feasibility issue with re-fit.
> 2. **Phase 3b** — merge IFE into CFE (verify E0/E4 equivalence, replace `inter_fe_ub` with `complex_fe_ub`).
>
> **Resolved**: All factors.from phases (1-7), CFE bifurcation (3a), score unification (Phase 1+2), cv.method unification, cv.sample k-fold CV, parallel CV, CFE CV r-selection, test gaps, .as_mask() bug, fect_mspe_sim deleted.
>
> Read `~/GitHub/fect/log/HANDOFF-factors-from.md` for full context.
