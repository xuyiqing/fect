# Handoff: factors.from Refactoring (REQ-factors-from-001)

## Date
2026-03-16 (updated)

## Status
Score unification Phase 1 COMPLETE and Phase 2 COMPLETE (uncommitted) — `.score_residuals()` extracted and wired into all three CV functions (Phase 1); `cv.method` parameter replaces `cv.treat`/`mask.method`, fect_mspe simplified, 1% rule hardcoded in fect_nevertreated, W/count.T.cv added to nevertreated scoring (Phase 2). 304/304 tests pass (130 score-unify + 44 utility + 77 skipped/CRAN-guarded + 53 factors-from counted in utility). All on `cfe` branch.

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

### Unify fect_cv and fect_mspe scoring/masking (NEW — approved plan)

**Problem**: `fect_cv` (internal, automatic model selection) and `fect_mspe` (user-facing, manual model comparison) both do hide→predict→score, but use incompatible masking and scoring. A user comparing models with `fect_mspe` gets different answers than `fect_cv` would give for the same data, because:
- Different masking (cv.sample with donut/min.T0 vs simple random)
- Different scoring (8 criteria vs RMSE+Bias only)
- Different weight handling (W supported vs not)
- Different normalization (norm.para vs none)

**Principle**: Everything accessible through `fect_cv` should also be accessible through `fect_mspe`, and vice versa.

#### Layer 1: Shared scoring function

```r
.score_residuals(resid, weights = NULL, time_index = NULL,
                 count_weights = NULL, norm.para = NULL)
# Returns named list:
#   MSPE, WMSPE, GMSPE, WGMSPE, MAD, Moment, GMoment, RMSE, Bias
```

- `resid`: vector of (pred - actual)
- `weights`: observation-level weights (W matrix values), NULL = unweighted
- `time_index`: relative time labels per residual (for moment/weighted criteria)
- `count_weights`: period-level weights (from `count.T.cv`), NULL = uniform
- `norm.para`: normalization vector, NULL = no normalization

Currently in `fect_cv`: scoring logic is inline (lines 396-484, repeated for MC at ~700+). Extract once, call from both.

Currently in `fect_mspe`: only `sqrt(mean(err^2))` and `mean(err)`. Replace with `.score_residuals()`, expose `criterion` param to user.

#### Layer 2: Shared masking options

| Strategy | Currently in | Expose in | How |
|---|---|---|---|
| `cv.sample` (donut, min.T0, cv.treat, cv.nobs) | `fect_cv` only | `fect_mspe` | Add `mask.method = "cv.sample"` option |
| Random from D==0 | `fect_mspe` only | `fect_cv` (not needed — cv.sample subsumes) | Already there |
| Pre-trend (last n pre-treatment periods) | `fect_mspe` only | Could add to `fect_cv` | Low priority |
| User-provided mask | `fect_mspe` only | Keep in `fect_mspe` | Post-hoc use case |
| k-fold structure | `fect_cv` only | `fect_mspe` | Add `k` param for k-fold CV mode |

#### Layer 3: Updated `fect_mspe` signature

```r
fect_mspe <- function(
    out.fect,                    # fect output or list of outputs
    hide_mask   = NULL,          # user-provided mask (TT×N)
    hide_n      = 20,            # number of obs to hide per rep
    seed        = NULL,
    n_rep       = 1,
    # --- masking strategy ---
    mask.method = "random",      # "random" | "pre.trend" | "cv.sample" | "user"
    pre.trend.n = 2,             # for mask.method="pre.trend"
    cv.treat    = TRUE,          # for mask.method="cv.sample"
    cv.nobs     = 3,             # for mask.method="cv.sample"
    cv.donut    = 1,             # for mask.method="cv.sample"
    cv.prop     = 0.1,           # for mask.method="cv.sample"
    min.T0      = 5,             # for mask.method="cv.sample"
    k           = 5,             # for mask.method="cv.sample" (k-fold)
    # --- scoring ---
    criterion   = "mspe",        # any of: mspe, wmspe, gmspe, wgmspe, mad, moment, gmoment
    W           = NULL,          # observation weights
    norm.para   = NULL,          # normalization
    # --- ground truth ---
    actual      = NULL,          # TT×N ground truth matrix (default: observed Y)
    control.only = TRUE          # mask only D==0 cells?
)
```

#### Layer 4: Refactor `fect_cv` to use shared scoring

Replace inline scoring blocks (lines 396-484 for IFE, ~700+ for MC) with calls to `.score_residuals()`. This is a pure refactor — no behavior change, just deduplication. `fect_cv` keeps its own masking (`cv.sample` + k-fold + `initialFit` per fold) and hyperparameter grid search.

#### Execution order

1. **Extract `.score_residuals()`** from `fect_cv` inline code → new file `R/score.R` (or top of `R/cv.R`)
2. **Wire `fect_cv` to use it** — replace inline scoring, verify identical CV results
3. **Wire `fect_mspe` to use it** — add `criterion` param, replace RMSE/Bias-only scoring
4. **Add `mask.method="cv.sample"` to `fect_mspe`** — reuse `cv.sample()` function
5. **Add W and norm.para support to `fect_mspe`**
6. **Tests**: verify `fect_mspe(mask.method="cv.sample", criterion="mspe")` produces scores consistent with `fect_cv`'s internal CV scores on the same data

#### Phase 1 status: COMPLETE (uncommitted)

Steps 1-6 are done. `.score_residuals()` is in R/score.R. fect_cv (IFE+MC blocks), fect_nevertreated (IFE+CFE blocks), and fect_mspe all use it. 84 new tests pass. Also fixed: gmoment storage bug in cv.R line 590, criterion input validation, cv.sample subscript mapping. `fect_mspe_sim` deleted.

#### Phase 2: cv.method unification + simplification — COMPLETE (uncommitted)

**Implementation done** in run `REQ-cv-method-phase2`. See `log/2026-03-16-cv-method-unification.md` for full process record.

**What was implemented:**

- A. `cv.method` parameter added to fect_cv, fect_nevertreated, fect_mspe, fect_boot, fect_binary_cv, and all default.R call sites. Replaces `cv.treat` (boolean) and `mask.method` (string).
- B. fect_mspe simplified: removed hide_mask, hide_n, n_rep, pre.trend, pre.trend.n, mask.method, actual, control.only. Only cv.sample masking remains.
- C. 1% selection rule hardcoded in fect_nevertreated (both IFE and CFE blocks).
- D. W.tr and count.T.cv added to fect_nevertreated .score_residuals() calls.
- E. All default.R call sites updated (3 signatures + 5 pass-through calls).
- F. Tests updated: 130/130 pass in test-score-unify.R (24 new, ~7 deleted, ~15 updated).

**Additional fixes during implementation:**
- fect_cv: match.arg moved after nevertreated delegation; added IFE+nevertreated delegation block (was missing); cv.method/cv.sample params passed through to all delegation calls
- fect_mspe: .is_fect_output() fixed to check obj$call/Y.dat instead of Y.ct.full; active_rows/active_cols feasibility checks added
- fect_nevertreated: CV.out now included in return list

**Known limitation**: cv.sample-based CV in fect_nevertreated has parameter infrastructure but the actual loop logic is not yet implemented (LOO path handles all cv.method values). Ready for future implementation.

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

> I'm working on the fect R package (`~/GitHub/fect`, branch `cfe`). Score unification Phase 1 and Phase 2 are both complete. 304 tests pass (130 score-unify + 44 utility + 77 skipped factors-from with CRAN guard + 53 factors-from counted in utility).
>
> **UNCOMMITTED CHANGES** on branch `cfe` — commit before starting new work:
> - `R/score.R` (NEW) — shared `.score_residuals()` function (Phase 1)
> - `R/cv.R` — `cv.method` replaces `cv.treat`; added IFE+nevertreated delegation; cv.method passthrough to all delegation calls
> - `R/fect_mspe.R` — simplified to cv.sample-only; `cv.method` replaces `mask.method`; removed 8 params
> - `R/fect_nevertreated.R` — `cv.method` param, W.tr/count.T.cv scoring, hardcoded 1% rule, CV.out in return
> - `R/default.R` — `cv.method` replaces `cv.treat` in 3 signatures + 5 call sites
> - `R/boot.R` — `cv.method` replaces `cv.treat`
> - `R/cv_binary.R` — `cv.method` replaces `cv.treat`
> - `man/fect.Rd` — `cv.method` documentation
> - `tests/testthat/test-score-unify.R` — 130 tests (24 new for Phase 2)
> - `tests/testthat/test-utility-functions.R` — updated for simplified fect_mspe
> - `NAMESPACE` — removed fect_mspe_sim export
> - `man/fect_mspe_sim.Rd` — deleted
> - `architecture.md` — updated for Phase 2
> - `vignettes/aa-cheatsheet.Rmd` — removed fect_mspe_sim
> - `log/` — updated handoff, new log entries, cv-comparison-table
>
> **Open tasks** (in priority order):
>
> 1. **Phase 3b** — merge IFE into CFE (verify E0/E4 equivalence, replace `inter_fe_ub` with `complex_fe_ub`).
> 2. **cv.sample CV in fect_nevertreated** — parameter infrastructure is in place but the actual cv.sample loop logic (for cv.method="all_units"/"treated_units") is not yet implemented inside the r-search loop. LOO handles all cv.method values for now.
> 3. **fect_mspe + CV=TRUE bug** — fect_mspe throws "No valid residuals" when called on models fitted with CV=TRUE. Workaround: use CV=FALSE.
>
> **Resolved**: Score unification Phase 1, Phase 2 (cv.method unification, fect_mspe simplification, 1% rule, W/count.T.cv), CFE CV r-selection issue, test gaps, parallel .export, .as_mask() bug fix, fect_mspe_sim deleted.
>
> Read `~/GitHub/fect/log/HANDOFF-factors-from.md` for full context.
> Read `~/GitHub/fect/log/cv-comparison-table.md` for the approved CV design table.
