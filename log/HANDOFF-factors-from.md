# Handoff: factors.from Refactoring (REQ-factors-from-001)

## Date
2026-03-16 (updated)

## Status
COMPLETE through Phase 3a + Category I bootstrap tests + CV gap fix + r=0 invariance test + Quarto book update + CFE CV r-selection fix + utility function tests — all committed and pushed to `cfe`. 261/261 tests pass (208 existing + 53 new). Quarto book renders (10/10 chapters).

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
| — | Utility function tests: fect_mspe, esplot, att.cumu, effect (53 new tests) | **Done** (`5460b13`) |

---

## Target repository
- **Repo**: xuyiqing/fect
- **Branch**: `cfe`
- **Local path**: `~/GitHub/fect`
- **HEAD**: `5460b13`

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

---

## Test results

**261/261 tests pass** (208 in `test-factors-from-refactor.R` + 53 in `test-utility-functions.R`)

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

### ~~Compare cv.R vs fect_mspe.R~~ RESOLVED — no consolidation
**Finding**: `fect_cv` (1570 lines, internal CV for r/λ selection) and `fect_mspe` (382 lines, post-estimation hide-and-refit) operate at different abstraction levels and share zero reusable code. `fect_cv` calls C++ directly (`inter_fe_ub`, `inter_fe_mc`), manipulates the II indicator matrix, and searches a hyperparameter grid. `fect_mspe` operates on user-facing `fect` output objects, sets Y to NA, and re-runs the full R pipeline. No consolidation needed.

### Refactor fect_mspe.R (NEW — approved plan)
`fect_mspe` and `fect_mspe_sim` share the hide→refit→score pattern but duplicate ~200 lines of code (helpers, loop, scoring). Refactor into shared core + thin wrappers.

**Architecture**:
```
# file-scope shared helpers (defined once, not copy-pasted)
.is_fect_output()
.get_last_matrix()
.as_mask()
.build_rerun_args()

# shared core
._refit_and_score(out_list, mask_mat, actual_mat, hide_n, n_rep, control.only, ...)
  → for each rep: hide Y at masked positions → do.call(fect) → compare pred vs actual → RMSE/Bias

# thin wrappers (own the masking strategy, delegate scoring to core)
fect_mspe(out.fect, hide_mask, hide_n, seed, n_rep,
          pre.trend, pre.trend.n, control.only=TRUE)
  → mask via: random control / pre-trend / user-provided
  → actual = observed Y matrix
  → calls ._refit_and_score()

fect_mspe_sim(out.fect, hide_mask, hide_mask_y0, hide_n, seed, n_rep,
              control.only=TRUE)
  → mask from user input
  → actual = Y0 matrix
  → calls ._refit_and_score()
```

**Consistency improvements**:
1. **Unified signature pattern** — both take same base args `(out.fect, hide_mask, hide_n, seed, n_rep)`, plus strategy-specific params
2. **Explicit ground truth** — `fect_mspe` uses observed Y, `fect_mspe_sim` uses Y0; both pass `actual_mat` to shared core
3. **Single `.build_rerun_args()`** — `fect_mspe_sim` currently duplicates this inline (lines 341-349)
4. **`control.only=TRUE` param** — both default to masking only `D==0` cells; `fect_mspe_sim` currently skips this filter, which could bite users
5. **Single formula recovery path** — use `.build_rerun_args()` in both (tries `eval(call$formula)`, falls back to `reformulate`)

**Masking strategies preserved**:
| Strategy | Wrapper | How |
|---|---|---|
| Random control | `fect_mspe`, `pre.trend=FALSE` | Sample from `D==0` cells |
| Pre-trend | `fect_mspe`, `pre.trend=TRUE` | Last `n` pre-treatment periods of treated units |
| User mask | both, `hide_mask` arg | User-supplied matrix, filtered by `D==0` when `control.only=TRUE` |
| Sim mask | `fect_mspe_sim` | External mask + Y0 ground truth |

**Prerequisite**: .as_mask() bug fix (already done — `028beca`).

### ~~CFE CV r-selection issue~~ RESOLVED
**Finding**: The CFE CV algorithm is correct. The original test was misspecified -- it used `make_cfe_z_data` (DGP with Z*gamma) but did not pass `Z = "Z"` to `fect()`. Without Z, factors absorb the unmodeled Z*gamma interaction (a rank-1 component), inflating r.cv by ~1. With Z properly specified, CFE CV correctly selects r=2 with clear MSPE U-shape (minimum at true r). Fixed test + added 3 new tests. 135/135 pass. See `log/2026-03-16-cfe-cv-rselect-fix.md` for full process record.

### Phase 3b: Merge IFE into CFE
- Verify test E0: `complex_fe_ub` with empty CFE arrays ≡ `inter_fe_ub`
- Verify test E4: `ife+nevertreated` ≡ `cfe+nevertreated` (no extras)
- Only then: replace `inter_fe_ub` calls with `complex_fe_ub` in `fect_nevertreated`
- Re-run full test suite to confirm no regressions

### ~~Test gaps discovered from book review~~ RESOLVED
1. ~~`fect_mspe()` and `fect_mspe_sim()` have no tests~~ **Done** (`5460b13`) — 18 tests for fect_mspe; 2 input validation tests for fect_mspe_sim (positive path blocked by .as_mask() bug, see below)
2. ~~`esplot()` only has basic invocation tests~~ **Done** (`5460b13`) — 10 tests covering connected, SE-derived CI, show.count, highlight, fill.gap, start0, only.pre/post
3. ~~`att.cumu()` relationship to `effect()` needs clarification~~ **Done** (`5460b13`) — 23 tests; att.cumu works on aggregate ATT time series (x$att), effect() works on unit-level matrices (x$eff, requires keep.sims=TRUE)
4. ~~Verify parallel=TRUE works with the updated .export list~~ **Done** — .export entries are redundant when .packages=c("fect") loads namespace; no change needed

### ~~fect_mspe_sim .as_mask() bug~~ FIXED (`028beca`)
Both copies of `.as_mask()` (in `fect_mspe` and `fect_mspe_sim`) fixed: `as.logical(mask)` → `matrix(as.logical(mask), nrow=TT, ncol=N)` to preserve matrix dimensions. Will be eliminated entirely by the refactor above (single shared helper).

### Nice-to-have
- ~~Add test for CFE parametric bootstrap with Z/Q/sfe parameters~~ **Done** (F1-F3)
- ~~Add test for r=0 invariance~~ **Done** (I12)
- ~~Verify parallel=TRUE works with the updated .export list~~ **Done** (no change needed)

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

> I'm working on the fect R package (`~/GitHub/fect`, branch `cfe`). All code changes and Quarto book updates are complete. 261/261 tests pass, 10/10 chapters render.
>
> **Open tasks** (in priority order):
>
> 1. **Refactor fect_mspe.R** — extract shared core from `fect_mspe` and `fect_mspe_sim` (hide→refit→score), unify signatures, add `control.only=TRUE` param, deduplicate helpers. Plan approved — see "Refactor fect_mspe.R" in open items.
> 2. **Phase 3b** — merge IFE into CFE (verify E0/E4 equivalence, replace `inter_fe_ub` with `complex_fe_ub`).
>
> **Resolved**: CFE CV r-selection issue, test gaps (53 new tests), parallel .export (no change needed), cv.R vs fect_mspe.R comparison (no consolidation), .as_mask() bug fix.
>
> Read `~/GitHub/fect/log/HANDOFF-factors-from.md` for full context.
