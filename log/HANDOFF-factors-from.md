# Handoff: factors.from Refactoring (REQ-factors-from-001)

## Date
2026-03-16 (updated)

## Status
COMPLETE through Phase 3a + Category I bootstrap tests + CV gap fix + r=0 invariance test + Quarto book update + CFE CV r-selection fix — all committed and pushed to `cfe`. 135/135 tests pass. Quarto book renders (10/10 chapters).

---

## What was done

Added `factors.from` and `em` parameters to the `fect` R package, rerouted `ife+nevertreated` through `fect_nevertreated` (renamed from `fect_gsynth`), removed obsolete methods (`polynomial`, `bspline`, `cfe_old`), implemented CFE bifurcation in `fect_nevertreated`, fixed CV routing, added comprehensive tests, and updated the Quarto user manual.

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
| — | CFE CV r-selection: diagnosed as test misspecification (missing Z="Z"), fixed test + added 3 new tests | **Done** (REQ-cfe-cv-rselect-001) |

---

## Target repository
- **Repo**: xuyiqing/fect
- **Branch**: `cfe`
- **Local path**: `~/GitHub/fect`
- **HEAD**: `fd2abd3`

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

---

## Test results

**135/135 tests pass** (`test-factors-from-refactor.R`)

| Category | Tests | What |
|----------|-------|------|
| Phase 6 | em + reroute | gsynth ≡ ife+nevertreated |
| Phase 1 | factors.from | Parameter acceptance, defaults, CV, bootstrap combos |
| Phase 2 | gsynth merger | Backward compat |
| Phase 3 | parametric bootstrap | fe/ife/cfe; mc error |
| Phase 5 | validation | Invalid values, guards |
| A1-A3 | Solver equivalence | complex_fe_ub vs inter_fe_ub |
| B1-B4 | Specification equivalence | ife ≡ cfe (no extras) |
| C1-C7 | Accuracy | ATT within tolerance for Z, Q, FE, full |
| D1-D5 | Validation guards | Error checks |
| E1-E5 | Output completeness | Fields, plot, print |
| F1-F4 | Bootstrap/inference | Jackknife, parametric, existing paths |
| G1-G3 | Cross-validation | r-selection, regression |
| H1-H5 | Edge cases | Single unit, no covariates, reversals |
| I1-I4 | Bootstrap basic | ife/cfe × em balanced |
| I5-I6 | Seed reproducibility | Same/different seeds |
| I7 | Different seeds SE | ife |
| I8 | ATT accuracy + CI | cfe+nevertreated |
| I9 | quantile.CI=TRUE | Bias-corrected CI |
| I10 | em no-op | em=TRUE ≡ em=FALSE for nevertreated |
| I11 | Unbalanced bootstrap | 5% dropped, _ub/EM path |
| I12 | r=0 invariance | factors.from no-op at r=0 |
| CV test | cfe+nevertreated CV | CV routing + r-selection |

---

## Quarto book update

**Chapter order** (new):
1. Welcome (index.qmd) — two research settings, estimator families
2. Get Started (01-start.Rmd)
3. Fect Main Program (02-fect.Rmd) — factors.from/em callouts, estimation vs inference separation
4. Complex Fixed Effects (07-cfe.Rmd) — nevertreated section, fect_mspe example
5. Fect Plot Options (03-plots.Rmd) — esplot() section added
6. Gsynth Program (04-gsynth.Rmd) — equivalence callout, original narrative preserved
7. New DID Methods (05-panel.Rmd)
8. Sensitivity Analysis (06-sens.Rmd)
9. Cheatsheet (aa-cheatsheet.Rmd) — factors.from/em rows, utility functions section

**Key content added**:
- "Two Research Settings" in Welcome: synth setting (never-treated controls, parametric bootstrap) vs DID setting (EM imputation, nonparametric/jackknife)
- `factors.from` parameter documented with examples
- `em` parameter: explained as no-op for nevertreated
- gsynth equivalence: `method="gsynth"` ≡ `method="ife"` + `factors.from="nevertreated"`
- `esplot()` standalone event study plotting
- `fect_mspe()` model comparison
- Utility functions: esplot, fect_mspe, fect_mspe_sim, get.cohort, att.cumu, effect
- "Download R code" dead links removed from all chapters
- Typos fixed: "estiamtes" → "estimates", "Weighitng" → "Weighting"

---

## Open items

### Compare cv.R vs fect_mspe.R
User requested comparing and potentially consolidating these two functions. `fect_cv` (1570 lines) is internal CV for factor number selection; `fect_mspe` (382 lines) is post-estimation MSPE validation (hide-and-refit). Different purposes but may share logic.

### ~~CFE CV r-selection issue~~ RESOLVED
**Finding**: The CFE CV algorithm is correct. The original test was misspecified -- it used `make_cfe_z_data` (DGP with Z*gamma) but did not pass `Z = "Z"` to `fect()`. Without Z, factors absorb the unmodeled Z*gamma interaction (a rank-1 component), inflating r.cv by ~1. With Z properly specified, CFE CV correctly selects r=2 with clear MSPE U-shape (minimum at true r). Fixed test + added 3 new tests. 135/135 pass. See `log/2026-03-16-cfe-cv-rselect-fix.md` for full process record.

### Phase 3b: Merge IFE into CFE
- Verify test E0: `complex_fe_ub` with empty CFE arrays ≡ `inter_fe_ub`
- Verify test E4: `ife+nevertreated` ≡ `cfe+nevertreated` (no extras)
- Only then: replace `inter_fe_ub` calls with `complex_fe_ub` in `fect_nevertreated`
- Re-run full test suite to confirm no regressions

### Test gaps discovered from book review
1. `fect_mspe()` and `fect_mspe_sim()` have no tests
2. `esplot()` only has basic invocation tests
3. `att.cumu()` relationship to `effect()` needs clarification
4. Verify parallel=TRUE works with the updated .export list

### Nice-to-have
- ~~Add test for CFE parametric bootstrap with Z/Q/sfe parameters~~ **Done** (F1-F3)
- ~~Add test for r=0 invariance~~ **Done** (I12)
- Verify parallel=TRUE works with the updated .export list

---

## Architecture: how it works

### Two research settings

| | Synth Setting | DID Setting |
|---|---|---|
| **Origin** | Synthetic control literature (gsynth) | DID/TWFE literature (fect) |
| **Estimation** | Never-treated controls learn time effects, project to treated | All not-yet-treated obs, EM imputation |
| **`factors.from`** | `"nevertreated"` | `"notyettreated"` (default) |
| **Inference** | Parametric bootstrap (small N_tr) | Nonparametric bootstrap / jackknife (large N) |
| **Data requirements** | Sufficient never-treated units, no treatment reversal | — |

### Estimator families

| Estimator | Family | `factors.from` | Parametric bootstrap | Nonparametric/Jackknife |
|---|---|---|---|---|
| ife (r≥1) | IFE | either | Yes | Yes |
| fe (ife r=0) | IFE | either | Yes | Yes |
| cfe | IFE | either | Yes | Yes |
| mc | MC | `"notyettreated"` only | No | Yes |
| gsynth | wrapper | `"nevertreated"` (auto) | Yes | Yes |

### The II matrix

The **II matrix** (estimation sample indicator, TT × N) is the single control point, constructed at line ~1479 of default.R:

```r
II <- I
II[which(D == 1)] <- 0  ## treated observations excluded

if (factors.from == "nevertreated") {
    D.cum <- apply(D, 2, function(vec) { cummax(vec) })
    D.unit.sum <- colSums(D.cum)
    co.never <- which(D.unit.sum == 0)
    ever.treated <- setdiff(1:ncol(II), co.never)
    II[, ever.treated] <- 0  ## zero out all non-never-treated columns
}
```

### Rerouting logic

When `method="ife"` + `factors.from="nevertreated"`, dispatch routes to `fect_nevertreated()` instead of `fect_fe()`. This function subsets to control units, estimates factors on never-treated only, then projects counterfactuals onto treated units.

When `method="gsynth"`, lines 488-493 auto-set `factors.from="nevertreated"`, then dispatch routes to `fect_nevertreated()`.

When `method="cfe"` + `factors.from="nevertreated"`, dispatch routes to `fect_nevertreated(..., method="cfe")` which uses `complex_fe_ub` with three-layer projection and block coordinate descent.

### CV routing

`fect_cv` delegates `cfe+nevertreated` to `fect_nevertreated(..., method="cfe", CV=TRUE)`, mirroring the gsynth delegation pattern. Both `fect_boot` and the direct `se=FALSE+CV=TRUE` path in `default.R` forward `factors.from` and all CFE-specific params to `fect_cv`.

---

## Context for new conversation

> I'm working on the fect R package (`~/GitHub/fect`, branch `cfe`). All code changes and Quarto book updates are complete. 135/135 tests pass, 10/10 chapters render.
>
> **Open tasks** (in priority order):
>
> 1. **Compare cv.R vs fect_mspe.R** — user wants to explore consolidation. `fect_cv` (internal CV for r-selection) vs `fect_mspe` (post-estimation hide-and-refit MSPE). Different purposes but may share logic.
> 2. **Test gaps** — fect_mspe/fect_mspe_sim have no tests; esplot has basic tests only; att.cumu needs clarification.
> 3. **Phase 3b** — merge IFE into CFE (verify E0/E4 equivalence, replace `inter_fe_ub` with `complex_fe_ub`).
>
> **Resolved**: CFE CV r-selection issue (was test misspecification, not algorithm bug).
>
> Read `~/GitHub/fect/log/HANDOFF-factors-from.md` for full context.
