# Handoff: factors.from Refactoring (REQ-factors-from-001)

## Date
2026-03-15 (updated)

## Status
COMPLETE through Phase 3a — CFE bifurcation in `fect_nevertreated` implemented with block coordinate descent, all 92 tests pass. draw.error() subsetting bug fixed.

---

## What was done

Added `factors.from` and `em` parameters to the `fect` R package, rerouted `ife+nevertreated` through `fect_nevertreated` (renamed from `fect_gsynth`), removed obsolete methods (`polynomial`, `bspline`, `cfe_old`), and added comprehensive tests.

| Phase | What | Status |
|-------|------|--------|
| 1 | `factors.from` parameter in all 3 signatures + II matrix filtering + output storage | **Done** (`53b461d`, `cafb66d`) |
| 2 | `method="gsynth"` auto-sets `factors.from="nevertreated"` + `em=FALSE` | **Done** (`37d516b`) |
| 3 | Parametric bootstrap unlocked for ife/fe/cfe (was gsynth-only) | **Done** (`070d574`) |
| 4 | Removed `method="polynomial"`, `"bspline"`, `"cfe_old"` entirely | **Done** (`953da8a`, `efa2af0`) |
| 5 | Input validation guards (invalid values, mc/both+nevertreated, em+notyettreated, insufficient units) | **Done** (`cafb66d`) |
| 6 | `em` parameter + reroute ife+nevertreated → `fect_nevertreated` | **Done** (`37d516b`) |
| 7 | Rename `fect_gsynth` → `fect_nevertreated` | **Done** (`fa9baf2`) |
| 3a | CFE bifurcation in `fect_nevertreated`: method param, `complex_fe_ub` solver, three-layer projection, block coordinate descent, boot.R fixes, 36 new tests | **Done** (uncommitted on `tianzhu`) |

---

## Target repository
- **Repo**: xuyiqing/fect
- **Branch**: `cfe`
- **Local path**: `~/GitHub/fect`
- **Working branch**: `tianzhu` (Phase 3a changes uncommitted)
- **Prior HEAD**: `fa9baf2` (Rename fect_gsynth to fect_nevertreated)

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

---

## Files modified

### R/default.R
| Location | Change |
|----------|--------|
| Line ~38 (fect signature) | Added `factors.from = "notyettreated"`, `em = TRUE` |
| Line ~111 (fect.formula signature) | Added `factors.from = "notyettreated"`, `em = TRUE` |
| Line ~217 (fect.formula passthrough) | Added `factors.from = factors.from`, `em = em` |
| Line ~291 (fect.default signature) | Added `factors.from = "notyettreated"`, `em = TRUE` |
| Lines ~470-493 (validation) | factors.from value check, mc/both+nevertreated guard, em+notyettreated guard, gsynth auto-override |
| Lines ~1479-1499 (II construction) | Nevertreated filtering with cummax + count guards |
| Lines ~2010-2055 (dispatch) | ife+nevertreated → `fect_nevertreated()` reroute |
| Line ~2105 (dispatch) | gsynth → `fect_nevertreated()` |
| Line ~2745 (output list) | Added `factors.from`, `em` to output |
| Throughout | Removed polynomial, bspline, cfe_old dispatch branches |

### R/boot.R
| Location | Change |
|----------|--------|
| Lines ~686-808 | `draw.error()` method-aware dispatch: gsynth → `fect_nevertreated`, ife → `fect_fe`, cfe → `fect_cfe` |
| Line ~856 | `.export` list: added `"fect_fe"`, `"fect_cfe"` |

### R/fect_nevertreated.R (renamed from R/fect_gsynth.R)
| Change | Description |
|--------|-------------|
| Renamed | `fect_gsynth.R` → `fect_nevertreated.R`, function name `fect_gsynth` → `fect_nevertreated` |

### R/fe.R, R/cfe.R
| Change | Description |
|--------|-------------|
| `eff.tr` added | Added `eff.tr` to return list so `draw.error()` can access treatment effect submatrix |

### R/cv.R
| Change | Description |
|--------|-------------|
| Updated references | `fect_gsynth` → `fect_nevertreated` |

### tests/testthat/test-factors-from-refactor.R
6 phases of acceptance tests (see Tests section below).

---

## Architecture: how it works

### Two orthogonal axes

| Axis | Options | Controls |
|------|---------|----------|
| **Estimator** | fe (r=0), ife (r>0), cfe, mc | The statistical model / solver |
| **Predictive routine** | notyettreated (default), nevertreated | Which observations estimate factors |
| **EM** | TRUE (default), FALSE | Missing data handling in estimation sample |

### The II matrix

The **II matrix** (estimation sample indicator, TT × N) is the single control point, constructed at line ~1479 of default.R:

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

### Rerouting logic

When `method="ife"` + `factors.from="nevertreated"`, dispatch at line ~2010 routes to `fect_nevertreated()` instead of `fect_fe()`. This function subsets to control units, estimates factors on never-treated only, then projects counterfactuals onto treated units.

When `method="gsynth"`, lines 488-493 auto-set `factors.from="nevertreated"` and `em=FALSE`, then dispatch at line ~2105 routes to `fect_nevertreated()`. Result: `method="gsynth"` and `method="ife", factors.from="nevertreated"` are now identical code paths.

---

## Tests

**File**: `tests/testthat/test-factors-from-refactor.R` (542 lines)

| Phase | Lines | What's tested |
|-------|-------|---------------|
| 6 (em + reroute) | 80-176 | gsynth ≡ ife+nevertreated, em=FALSE validation, em auto-setting, accuracy |
| 1 (factors.from) | 182-292 | Parameter acceptance, defaults, combinations with CV and bootstrap |
| 2 (gsynth merger) | 315-375 | Backward compat, gsynth as ife+nevertreated |
| 3 (parametric bootstrap) | 380-440 | Parametric bootstrap for fe/ife/cfe; mc error |
| 5 (validation) | 451-541 | Invalid values, mc/both+nevertreated, insufficient units, output recording |

**Test results** (Phase 1-7, from commit `070d574`): 38/39 pass; Phase 2a numerical tolerance is the sole remaining issue.

**Test results** (Phase 3a): 92/92 pass (36 new tests across categories A-H, plus 56 existing). All accuracy tests (C1-C7) within tolerance. Bootstrap tests (F1-F4) all pass. Edge cases (H1-H5) all handled.

---

## ~~Remaining bug: draw.error() CFE subsetting~~ FIXED in Phase 3a

**Status**: **FIXED** — Phase 3a added CFE array subsetting in `draw.error()` (boot.R) for all 5 arrays (X.extra.FE, X.Z, X.Q, X.gamma, X.kappa) before dispatch. Also added `cfe+nevertreated` and `ife+nevertreated` dispatch branches. Tests F1-F4 now cover bootstrap paths.

---

## Nice-to-have (non-blocking)

- ~~Add test for CFE parametric bootstrap with Z/Q/sfe parameters~~ **Done** (tests F1-F3)
- Add test for r=0 invariance (factors.from should have no effect when r=0)
- Verify parallel=TRUE works with the updated .export list

---

## Roadmap: CFE + nevertreated (Phase 3)

### Goal

Extend `fect_nevertreated` with a CFE bifurcation so that `method="cfe"` + `factors.from="nevertreated"` routes through the nevertreated estimator using `complex_fe_ub` on never-treated controls. This completes the separation of estimator × predictive routine.

### The CFE Model

```
Y_it = mu + alpha_i + xi_t + X_it·beta + Z_i·gamma_g(t) + Q_t·kappa_g(i)
        + FE_extra + F_t·lambda_i + epsilon_it
```

### Three-Layer Projection for Nevertreated

Parameters split into three layers based on how they're estimated for treated units:

| Layer | Parameters | Estimated from | Projection |
|-------|-----------|----------------|------------|
| 1: Shared | mu, xi, beta, gamma, Type-B extra FE, F (factors) | co estimation | Apply directly to tr |
| 2: Unit-specific | alpha, kappa, lambda, Type-A extra FE | tr pre-treatment | OLS/demeaning on pre-treatment residuals after removing Layer 1 |
| 3: Counterfactual | Y.ct.tr | — | Layer 1 + Layer 2 |

### Two Types of Extra Fixed Effects

**Type A (unit-nesting)**: Every level contains exclusively treated OR exclusively control units. No overlap. Structurally like alpha — must estimate from pre-treatment of treated.

**Type B (shared)**: Some levels contain both treated and control units. Can estimate from co and apply to tr. Must validate all tr levels exist in co; stop() with clear message if not.

Detection: at runtime, compare `unique(X.extra.FE[1, tr, k])` vs `unique(X.extra.FE[1, co, k])`. If `intersect` is empty → Type A. Otherwise → Type B with missing-level check.

### kappa Is Unit-Specific

kappa (unit-specific time trends via Q basis) is indexed by units (like alpha). In the nevertreated setting, kappa for treated units cannot be learned from controls — it must be estimated from pre-treatment residuals of treated units, after removing Layer 1 effects:

```
kappa.tr = (Q.pre' Q.pre)^{-1} Q.pre' U.tr.pre
```

where U.tr.pre is the pre-treatment residual after subtracting mu, xi, beta·X, gamma·Z, and Type-B extra FE.

### CV for CFE

Same leave-one-out structure as IFE nevertreated. For each candidate r, run `complex_fe_ub` on co-only data, project onto treated pre-treatment, compute MSPE. In the nevertreated setting the co panel may be complete, so `complex_fe_ub` effectively runs without EM.

### Phased Approach

**Phase 3a**: Implement CFE bifurcation in `fect_nevertreated` — **DONE** (uncommitted on `tianzhu`)
- Added `method` parameter and 7 CFE params to signature (line 39)
- Added CFE array co/tr subsetting via `.split_array` helper (lines 106-162)
- Added Type A/B extra FE classification with validation (stop if missing levels)
- Added `complex_fe_ub` solver path with CV loop (lines 641-1002)
- Added Layer 1 subtraction (mu, xi, beta, gamma, Type-B extra FE) in projection
- Added Layer 2 block coordinate descent for unit-specific parameters (lines 1065-1285):
  - Iterates alpha, kappa, Type-A FE, and lambda jointly to convergence
  - When `r.cv > 0`, alpha is embedded in augmented factor matrix (avoids double subtraction)
  - Convergence: max |change| < 1e-8, typically 2-5 iterations
- Added dispatch in `default.R` for `cfe + nevertreated` (SE and non-SE paths)
- Fixed `boot.R` draw.error() subsetting for all 5 CFE arrays
- Added 7 helper functions: `.estimate_kappa_fit`, `.estimate_alpha`, `.estimate_typeA_fit`, `.estimate_lambda_fit`, `.reconstruct_gamma_fit_tr`, `.reconstruct_kappa_fit`, `.extract_and_apply_typeB_fe`
- 36 new tests (92 total assertions), all pass

**Phase 3b**: Merge IFE into CFE (only after Phase 3a validates)
- Verify test E0: `complex_fe_ub` with empty CFE ≡ `inter_fe_ub`
- Verify test E4: `ife+nevertreated` ≡ `cfe+nevertreated` (no extras)
- Only then: replace `inter_fe_ub` calls with `complex_fe_ub` in `fect_nevertreated`
- Re-run full test suite to confirm no regressions

### Test Plan (35 tests across 8 categories)

Tests are designed so that each catches a specific failure mode. Priority order for implementation: B1 → C1 → C2 → C3 → D1 → D4 → F2 → A1.

#### DGP Helpers Needed

```r
## DGP with Z/gamma (time-invariant covariates + grouped time coefficients)
make_cfe_z_data(N, TT, Ntr, tau, r, seed)
  # Y_it = alpha_i + xi_t + Z_i * gamma_g(t) + F_t * lambda_i + tau*D + eps

## DGP with Q/kappa (unit-specific time trends)
make_cfe_q_data(N, TT, Ntr, tau, r, seed)
  # Y_it = alpha_i + xi_t + Q_t * kappa_i + F_t * lambda_i + tau*D + eps

## DGP with extra FE (shared/Type-B: industry spans treated and control)
make_cfe_fe_data(N, TT, Ntr, tau, r, seed)
  # Y_it = alpha_i + xi_t + industry_FE + F_t * lambda_i + tau*D + eps

## DGP with extra FE (unit-nesting/Type-A: industry perfectly partitions tr/co)
make_cfe_fe_nesting_data(N, TT, Ntr, tau, r, seed)
  # Same but each industry level is exclusively treated or control

## DGP with all CFE components
make_cfe_full_data(N, TT, Ntr, tau, r, seed)
  # Y_it = alpha_i + xi_t + Z_i*gamma_g(t) + Q_t*kappa_i + industry_FE
  #         + F_t*lambda_i + tau*D + eps
```

#### Category A: Solver Equivalence

Catches: regression in unified solver, mismatch between CFE and IFE backends.

| Test | DGP | What | Expected |
|------|-----|------|----------|
| A1 | `make_factor_data` (r=2) | `complex_fe_ub` with empty CFE arrays vs `inter_fe_ub` on same co-only data | att.avg identical (tol ~1e-3) |
| A2 | `make_factor_data` (r=0) | Same as A1 but r=0 (FE only) | att.avg identical |
| A3 | `make_staggered_data` (unbalanced) | Same as A1 with missing observations | att.avg identical |

#### Category B: Specification Equivalence

Catches: projection bugs where CFE path doesn't reduce to IFE when extras are empty.

| Test | DGP | What | Expected |
|------|-----|------|----------|
| B1 | `make_factor_data` | `ife+nevertreated` ≡ `cfe+nevertreated` (no Z/Q/extra FE) | att.avg identical |
| B2 | `make_factor_data` (r=0) | `fe+nevertreated` ≡ `cfe+nevertreated` (r=0, no extras) | att.avg identical |
| B3 | `make_cfe_z_data` | `cfe+nevertreated` run twice (same seed, same params) | att.avg identical |
| B4 | `make_factor_data` | `gsynth` ≡ `cfe+nevertreated` (no extras) | att.avg identical |

#### Category C: Accuracy

Catches: systematic bias in projection — wrong gamma, kappa, extra FE, or compound errors.

| Test | DGP | What | Expected |
|------|-----|------|----------|
| C1 | `make_cfe_z_data` (large N) | `cfe+nevertreated` with Z, r>0 | att.avg within 0.5 of true tau |
| C2 | `make_cfe_q_data` (large N) | `cfe+nevertreated` with Q, r>0 | att.avg within 0.5 of true tau |
| C3 | `make_cfe_fe_data` (Type-B, large N) | `cfe+nevertreated` with shared extra FE | att.avg within 0.5 of true tau |
| C4 | `make_cfe_fe_nesting_data` (Type-A, large N) | `cfe+nevertreated` with unit-nesting extra FE | att.avg within 0.5 of true tau |
| C5 | `make_cfe_full_data` (all components, large N) | `cfe+nevertreated` with everything | att.avg within 0.5 of true tau |
| C6 | `make_cfe_z_data` | `cfe+notyettreated` vs `cfe+nevertreated` | Both within 1.0 of true tau |
| C7 | `make_cfe_z_data` (r=0) | `cfe+nevertreated` with Z, r=0 | att.avg within 0.5 of tau |

#### Category D: Validation Guards

Catches: missing error checks that lead to silent garbage or crashes.

| Test | DGP | What | Expected |
|------|-----|------|----------|
| D1 | Custom (Type-B FE, missing level) | Treated units have an industry level absent from co | Informative stop() naming the missing level |
| D2 | Custom (all treated) | 0 never-treated units | Error |
| D3 | Custom (too few co) | Nco < r+1 | Error |
| D4 | Valid input | `method="cfe"` + `factors.from="nevertreated"` | Accepted, no error |
| D5 | Invalid combo | `method="mc"` + `factors.from="nevertreated"` | Error |

#### Category E: Output Completeness

Catches: missing output fields that break downstream consumers (plot, print, effect).

| Test | DGP | What | Expected |
|------|-----|------|----------|
| E1 | `make_cfe_z_data` | Output has gamma, kappa fields | Not NULL |
| E2 | `make_cfe_z_data` | Output has factors.from field | `== "nevertreated"` |
| E3 | `make_cfe_z_data` | Output has Y.ct, eff, att.avg, att (dynamic) | All non-NULL, correct dimensions |
| E4 | `make_cfe_z_data` | `plot(out)` doesn't error | No error |
| E5 | `make_cfe_z_data` | `print(out)` doesn't error | No error |

#### Category F: Bootstrap / Inference

Catches: draw.error() dimension mismatch, bootstrap incompatibility with nevertreated+CFE.

| Test | DGP | What | Expected |
|------|-----|------|----------|
| F1 | `make_cfe_z_data` | `cfe+nevertreated`, se=TRUE, vartype="jackknife", nboots=30 | Runs, SE not NA |
| F2 | `make_cfe_z_data` | `cfe+nevertreated`, se=TRUE, vartype="parametric", nboots=30 | Runs, SE not NA |
| F3 | `make_cfe_q_data` | `cfe+notyettreated`, se=TRUE, vartype="parametric", Z/Q params | Runs, SE not NA |
| F4 | `make_factor_data` | `ife+nevertreated`, se=TRUE, nboots=30 | Same SE as before change |

#### Category G: Cross-Validation

Catches: CV loop crashes or selects wrong r with complex_fe_ub.

| Test | DGP | What | Expected |
|------|-----|------|----------|
| G1 | `make_cfe_z_data` (r=2) | `cfe+nevertreated`, CV=TRUE, r=0 | Selects r.cv >= 0, runs |
| G2 | `make_cfe_z_data` (r=0) | `cfe+nevertreated`, CV=TRUE, r=0 on no-factor data | r.cv = 0 |
| G3 | `make_factor_data` | `ife+nevertreated` CV result unchanged | Same r.cv as before |

#### Category H: Edge Cases

Catches: corner cases with small panels, single units, no covariates, reversals.

| Test | DGP | What | Expected |
|------|-----|------|----------|
| H1 | Custom (Ntr=1) | `cfe+nevertreated`, single treated unit | Runs |
| H2 | Custom (Nco=1) | `cfe+nevertreated`, single control unit, r=0 | Runs or clear error |
| H3 | `make_cfe_z_data` | `cfe+nevertreated`, r=0 (no factors) | Runs, reasonable ATT |
| H4 | Custom (no covariates) | `cfe+nevertreated`, X=NULL, Z=NULL, Q=NULL | Runs |
| H5 | `make_cfe_full_data` (reversals) | `cfe+nevertreated` with treatment reversals | Runs |

### Files to Modify

| File | Changes |
|------|---------|
| `R/fect_nevertreated.R` | Add method param, CFE params, subsetting, extra FE classification, solver bifurcation, projection layers, output fields |
| `R/default.R` | Add dispatch for `cfe + nevertreated` → `fect_nevertreated(..., method="cfe")` |
| `R/boot.R` | Fix draw.error() subsetting (existing bug); add cfe+nevertreated path |
| `tests/testthat/test-factors-from-refactor.R` | Add Phase 3 DGP helpers + 35 tests |

### Risks (updated with resolutions)

1. **Extra FE projection for Type B**: ~~Need to extract group coefficients.~~ **Resolved** — used "reconstruct non-FE fit and subtract from total fit" approach. Helper `.extract_and_apply_typeB_fe` computes FE contribution as `est.co$fit - fit_no_fe.co`, then applies group means to treated.

2. **kappa.tr estimation order**: ~~Sequential estimation may introduce biases.~~ **Resolved** — replaced sequential estimation with block coordinate descent (lines 1065-1285). All unit-specific parameters (alpha, kappa, Type-A FE, lambda) are iterated jointly to convergence. Key insight: when `r.cv > 0`, alpha must be embedded in the augmented factor matrix to avoid double subtraction (discovered and fixed during implementation).

3. **Numerical equivalence in Phase 3b**: `complex_fe_ub` and `inter_fe_ub` use different iteration strategies. Small numerical differences are possible even with identical inputs. Tests A1-A3 use tolerance ~0.1 for ATT comparison. Phase 3b should tighten this after verifying empty-CFE equivalence.

4. **fect_cv gap** (new): `fect_cv` in cv.R does not handle `cfe+nevertreated`. With `se=TRUE, CV=TRUE`, CV routes through `fect_cv` which sends `cfe` to `fect_cfe` (notyettreated). The internal CV loop in `fect_nevertreated` handles the `se=FALSE` path correctly.

---

## Design notes

### gsynth as a predictive routine, not an estimator

**gsynth is NOT an estimator. It is a predictive routine that restricts factor estimation to never-treated units only.** The actual estimator is ife, cfe, etc.

`fect_fe` (ife path) passes the full TT × N matrix to `inter_fe_ub` with treated columns zeroed in II. The C++ solver estimates unit FEs for all N units, including treated units with zero observations — poorly identified.

`fect_nevertreated` (nevertreated path) subsets to control units first, estimates on Nco columns only, then projects counterfactuals. This is the correct approach.

Empirical divergence: ~0.04–0.17 ATT gap between the two paths. Structural, not numerical.

### Future work (Phase 3)

Extend `fect_nevertreated` to accept cfe as the estimator, completing the separation: predictive routine (nevertreated) × estimator (ife/cfe).

---

## Context for new conversation

> I'm working on the fect R package (`~/GitHub/fect`, branch `tianzhu`). Phase 3a (CFE bifurcation in `fect_nevertreated`) is implemented and tested (92/92 pass) but uncommitted. The key addition: `method="cfe"` + `factors.from="nevertreated"` now routes through `fect_nevertreated` using `complex_fe_ub` on never-treated controls, with a three-layer projection and block coordinate descent for treated unit-specific parameters. The draw.error() subsetting bug is fixed. Next step is Phase 3b (merge IFE into CFE) or commit + ship. Read `~/GitHub/fect/log/HANDOFF-factors-from.md` for full context.
