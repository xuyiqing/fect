# Handoff: factors.from Refactoring (REQ-factors-from-001)

## Date
2026-03-15 (updated)

## Status
COMPLETE through Phase 3a + Category I bootstrap tests — all committed and pushed to `cfe`. 125/125 tests pass. Next: Phase 3b (merge IFE into CFE).

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
| 3a | CFE bifurcation in `fect_nevertreated`: method param, `complex_fe_ub` solver, three-layer projection, block coordinate descent, boot.R fixes, 36 new tests | **Done** (`c0cab63`, pushed to `cfe`) |

---

## Target repository
- **Repo**: xuyiqing/fect
- **Branch**: `cfe`
- **Local path**: `~/GitHub/fect`
- **Working branch**: `tianzhu`
- **Prior HEAD**: `c0cab63` (Phase 3a CFE bifurcation)

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

**Test results** (Category I, `5136a03`): 125/125 pass. Added I9-I11 (quantile.CI, em no-op, unbalanced bootstrap). Fixed `est.att.avg` → `est.avg` field name bug in I8 and I9.

---

## ~~Remaining bug: draw.error() CFE subsetting~~ FIXED in Phase 3a

**Status**: **FIXED** — Phase 3a added CFE array subsetting in `draw.error()` (boot.R) for all 5 arrays (X.extra.FE, X.Z, X.Q, X.gamma, X.kappa) before dispatch. Also added `cfe+nevertreated` and `ife+nevertreated` dispatch branches. Tests F1-F4 now cover bootstrap paths.

---

## Open items

### Bootstrap test suite (Category I) — DONE (`5136a03`)

All 11 Category I tests pass (125/125 total). Covers:
- I1-I4: ife/cfe × em=TRUE/FALSE parametric bootstrap (balanced)
- I5-I6: Seed reproducibility (same seed → same SE; different seeds → different SE)
- I7: Different seeds, different SE (ife)
- I8: ATT accuracy + CI coverage (cfe+nevertreated)
- I9: `quantile.CI=TRUE` bias-corrected reflection CI (cfe+nevertreated)
- I10: `em=TRUE` vs `em=FALSE` identical for nevertreated (confirms em is no-op)
- I11: Unbalanced data (5% dropped) forces `_ub`/EM path in `draw.error()` for both IFE and CFE

**Bug fixed**: `est.att.avg` → `est.avg` in I8 and I9. The fect output field for aggregate ATT with CI is `est.avg` (from `fect_boot()` at boot.R line 4653), not `est.att.avg`.

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

**Phase 3a**: Implement CFE bifurcation in `fect_nevertreated` — **DONE** (`c0cab63`, pushed to `cfe`)

#### R/fect_nevertreated.R changes
- **Signature** (line 39): Added 8 params: `method="ife"`, `X.extra.FE`, `X.Z`, `X.Q`, `X.gamma`, `X.kappa`, `Zgamma.id`, `kappaQ.id`
- **CFE array subsetting** (lines 106-162): co/tr splitting for all 5 CFE arrays via `.split_array`. Type A/B extra FE classification with validation (stop if treated levels missing from controls)
- **Method bifurcation** (line 186): Existing IFE wrapped in `if (method == "ife")`, new `else if (method == "cfe")` block at line 641
- **CFE path** (lines 641-1002): `initialFit()` always called (even balanced), CV loop with `complex_fe_ub`, final fit, `est.co.fect` for equivalence test, three-layer projection
- **Layer 2 BCD** (lines ~1065-1285): Block coordinate descent over alpha, kappa, Type-A FE, lambda. Key: when `has_factor && has_alpha`, alpha embedded in augmented `F.hat.aug = [F, 1]` — `lambda_fit` already contains alpha, so separate `alpha_mat` NOT subtracted (`alpha_in_lambda` flag). Convergence: max|change| < 1e-8, typically 2-5 iters
- **r=0 path**: alpha.tr estimation parallel to IFE r=0 path
- **Method string** (line 1788): CFE outputs `"cfe"` not `"gsynth"`
- **Output** (line 1837): Added `gamma`, `kappa` fields
- **7 helper functions** (lines ~1940-2100): `.estimate_kappa_fit`, `.estimate_alpha`, `.estimate_typeA_fit`, `.estimate_lambda_fit`, `.reconstruct_gamma_fit_tr`, `.reconstruct_kappa_fit`, `.extract_and_apply_typeB_fe`

#### R/default.R changes
- **II zeroing** (line 1494): Added `(method == "cfe" && factors.from == "nevertreated")` — II NOT pre-zeroed for cfe+nevertreated
- **SE==FALSE dispatch** (line 2070): `cfe+nevertreated` branch BEFORE existing `cfe` branch → `fect_nevertreated(..., method="cfe")`
- **SE==TRUE dispatch** (line 2266): Added `factors.from = factors.from` to `fect_boot()` call

#### R/boot.R changes
- **Signature** (line 118): Added `factors.from = "notyettreated"` param
- **Initial estimation** (line 254): `cfe+nevertreated` branch → `fect_nevertreated(..., method="cfe")`
- **draw.error() subsetting** (lines 760-786): CFE array subsetting for all 5 arrays before dispatch
- **draw.error() dispatch** (lines 788-840): `ife+nevertreated` → `fect_nevertreated(..., method="ife")`, `cfe+nevertreated` → `fect_nevertreated(..., method="cfe")`, existing `cfe` path uses subsetted arrays
- **Parallel export** (line 912): Added `.reconstruct_gamma_fit_tr`, `.reconstruct_kappa_fit`, `.extract_and_apply_typeB_fe`

#### Key design decisions
1. `complex_fe_ub` returns sigma2/IC/PC/residuals/validX directly (verified from C++ `cfe.cpp` lines 195-199) — no post-hoc computation needed
2. CFE path always calls `initialFit` (matching `fect_cfe` pattern) — ensures Y0.co available
3. Type-B FE extraction: `est.co$fit - fit_no_fe.co` → group means → apply to treated
4. CV loop: same leave-one-out as IFE, kappa NOT subtracted in CV (MSPE primarily sensitive to r)
5. Gamma/kappa reconstruction: map group labels to coefficient rows, compute Z*gamma or Q*kappa

#### Deviations from spec
1. No post-hoc sigma2/IC/PC — `complex_fe_ub` returns them
2. Field name `est.co.best$residuals` not `$e`
3. Helper functions are dot-prefixed at file scope (accessible from boot.R parallel workers)
4. `fect_cv` (cv.R) NOT updated for `cfe+nevertreated` — known gap

#### Post-implementation fixes
1. **`.reconstruct_kappa_fit` dimension mismatch**: 3D array guard — if `length(dim()) == 3`, extract `[,,1]`
2. **Sequential estimation bias → BCD**: Single-pass estimation caused kappa error ~0.22. Replaced with block coordinate descent (converges in 2-5 iters)
3. **Alpha double subtraction**: `F.hat.aug = [F, 1]` means `lambda_fit` contains alpha. Subtracting separate `alpha_mat` caused divergence (~1e17). Fixed with `alpha_in_lambda` flag

#### Tests
- 36 new tests (92 total), all pass
- Categories A-H covering solver equivalence, specification equivalence, accuracy, validation, output, bootstrap, CV, edge cases

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

> I'm working on the fect R package (`~/GitHub/fect`, branch `cfe`). Phase 3a (CFE bifurcation) and Category I bootstrap tests are complete. 125/125 tests pass (`5136a03`, pushed to `cfe`).
>
> **Next step**: Phase 3b — merge IFE into CFE path inside `fect_nevertreated`.
>
> 1. Verify test E0: `complex_fe_ub` with empty CFE arrays ≡ `inter_fe_ub` on same data
> 2. Verify test E4: `ife+nevertreated` ≡ `cfe+nevertreated` (no extras)
> 3. Only then: replace `inter_fe_ub` calls with `complex_fe_ub` in `fect_nevertreated`
> 4. Re-run full test suite to confirm no regressions
> 5. Also open: `fect_cv` gap (cv.R doesn't handle `cfe+nevertreated`)
>
> Read `~/GitHub/fect/log/HANDOFF-factors-from.md` for full context.
