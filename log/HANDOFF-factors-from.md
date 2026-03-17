# Handoff: factors.from Refactoring (REQ-factors-from-001)

## Date
2026-03-17 (updated)

## Status
All phases COMPLETE and PUSHED to `cfe`. CFE convergence conditioning done (component-wise + denominator safety + stop-burnin). 565/565 tests pass (215 score-unify + 44 utility + 98 book-claims + 77 factors-from + 131 others). Quarto book restructured: 13 chapters, all render. Phase 2 (new content for model-selection chapter, factors.from deep dive) pending.

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

### ~~fect_mspe + CV=TRUE bug~~ RESOLVED
Fixed in `c2db08c`: force `CV=FALSE`, `se=FALSE`, and `r=r.cv` when re-fitting masked data inside fect_mspe.

### Convergence conditioning — COMPLETE (IFE + CFE)
- R centering in `.estimate_co()` removes grand mean before solver — tol applies to variation, not level
- C++ component-wise convergence in `fe_ad_inter_iter` / `fe_ad_inter_covar_iter` (ife_sub.cpp) — tracks factor convergence independently
- C++ component-wise convergence in `cfe_iter` (cfe_sub.cpp) — `max(dif, dif_inter)` for interactive FE, `+1e-10` denominator safety on all convergence ratios, stop-burnin transition (`7b33c98`)
- Two-tier CV tolerance (`cv_tol = max(tol, 1e-3)`) for CV speed
- Result: tol=1e-3 now gives 43-2249x better component accuracy. Zero speed regression (confirmed for both IFE and CFE paths).

### Phase 3b: R wrapper + solver equivalence — COMPLETE
- `.estimate_co()` wrapper added: dispatches to `inter_fe`/`inter_fe_ub` (IFE) or `complex_fe_ub` (CFE)
- C++ fix #1: guard alpha/xi access in `cfe_sub.cpp` — prevents crash with empty CFE arrays
- C++ fix #2: joint `ife()` fallback in `cfe_iter` — when no CFE components (Z, Q, gamma, kappa, extra FE all empty), delegates to `ife()` for joint FE+factor estimation instead of separate BCD steps. Achieves **exact numerical equivalence** (`max diff = 0.00`) with `inter_fe_ub` across all `force` values (0,1,2,3), balanced and unbalanced, with and without covariates. Iteration counts also match.
- **Root cause of prior divergence**: `cfe_iter` decomposed FE+factors into two separate BCD substeps with two imputation calls, while `inter_fe_ub` solved them jointly via `ife()`.
- **Speed overhead** (complex_fe_ub vs inter_fe_ub, after fix): 2% for N=1000 unbalanced, 7% for N=500 unbalanced. 4-5x for tiny balanced panels (N=100) due to CFE scaffolding.
- **Decision**: keep wrapper. Merge to single solver is now technically feasible (perfect equivalence) but the 2-7% overhead on typical panels is not justified given zero user-facing benefit.

---

## Context for new conversation

> I'm working on the fect R package (`~/GitHub/fect`, branch `cfe`). 590/590 tests pass. All code work complete. Currently in Quarto book restructure — user is reviewing.
>
> **Key architecture**: `.estimate_co()` wrapper dispatches to `inter_fe`/`inter_fe_ub` (IFE) or `complex_fe_ub` (CFE). Component-wise convergence applied to both IFE and CFE paths. Default parallel cores capped at 8 with boxed runtime message.
>
> **Plot refactor**: Gap plots now distinguish pre-treatment (gray, dashed in connected mode) from post-treatment (black, solid). New params `pre.color`/`post.color` in both `esplot()` and `plot.fect()`. `esplot()` accepts fect objects directly. Defaults harmonized between the two functions.
>
> **Quarto book**: 13 chapters, all render. Restructured: Ch2 (FE/imputation, uses `simdata1`), Ch3 (IFE/MC, uses `simdata2`), Ch4 (CFE), Ch5 (diagnostics), Ch6 (model selection — new, full chapter), Ch7 (plots), Ch8 (gsynth + CFE nevertreated). Nevertreated content removed from Ch2-4. New datasets: `simdata1` (parallel trends), `simdata2` (with factors, = original `simdata`).
>
> **Uncommitted changes**: Everything from today's session. Includes: cfe_sub.cpp convergence fix (already pushed as `7b33c98`), test-book-claims.R (pushed as `7418aa0`), plus plot refactor, core cap, boxed message, all Quarto changes, simdata1/simdata2. Old chapter files (03-plots, 04-gsynth, 05-panel, 06-sens, 07-cfe) still exist — delete after user confirms.
>
> **Next steps**: User manual review of rendered book (Phase 3). User may have more Quarto instructions. Commit and push when ready.
>
> Read `~/GitHub/fect/log/HANDOFF-factors-from.md` for full context. See `log/update-20260317.md` for today's session.
