# Change log: merging fect 1.1.10 features into the binary/rolling-CV branch

Date: 2026-09-21
Source of ported behaviour: `../fect` (fect 1.1.10, Licheng Liu's binary revision)
Target: this package (`fect` 2.4.5, "development — binary IFE")
Empirical check script: `../code/binary_empirical.R` (Authoritarian exchange-rate data)

Three features were ported or added. Each section lists the motivation, the
semantics, the files touched, and how it was verified. Short entries were
also added to `NEWS.md`.

---

## 1. Time-varying group indicator for cohort (group) effects

### Motivation
`fect(..., group = "g")` previously stopped with
`A unit in different periods should have the same group index` whenever the
group variable changed within a unit. fect 1.1.10 accepted such indicators
(e.g. regime type by country-year) and reported cohort-specific ATTs.

### Semantics
* **Time-invariant group** (one index per unit): unchanged. Cohort ATT plus
  cohort-specific dynamic effects (`group.output`, `est.group.output`) are
  reported and `plot(show.group = )` works.
* **Time-varying group** (index changes within a unit): the cohort ATT for
  level g is the mean of the treated-cell effects whose group index equals g
  (cell-level definition, identical to fect 1.1.10). Cohort dynamics are not
  defined and are dropped: `group.output` and `est.group.output` are `NULL`,
  and `plot(show.group = )` errors with an explanatory message.
* A message is printed when a time-varying indicator is detected.
* New fit-object fields: `group.time.varying` (logical); `group.att` is now
  named by group level.

### Files
| File | Change |
|---|---|
| `R/default.R` | Replaced the hard stop with detection of a time-varying indicator (`n.group.per.unit > 1`). Missing cells (coded 0) are filled with the unit's index only in the time-invariant case. Level set `g.level` excludes the missing-cell code 0. At output assembly: names `group.att`, stores `group.time.varying`, and for the time-varying case sets `group.output`, `est.group.output`, and `pre.est.group.output` to `NULL`. |
| `R/plot.R` | `show.group` now errors when `group.time.varying` is `TRUE`. |
| `R/print.R` | New "Cohort (group) ATT" section printed whenever `group` was supplied (with or without uncertainty estimates). |
| `man/fect.Rd` | `group` argument documented for both regimes. |
| `NEWS.md` | Entry added. |
| `tests/testthat/test-group-time-varying.R` | New: linear and binary fits, unbalanced panels, jackknife inference, unchanged time-invariant behaviour, plot guard, print output. |

### Verification
* Linear two-way FE, jackknife, `group = "htmil"` (script model 7): cohort ATTs
  0.0241 and 0.1071 with standard errors identical to the stored fect 1.1.10
  results in `../results_emp/fixed.RData`.
* Hand computation `tapply(eff[D==1], G[D==1], mean)` matches `group.att`
  in every test case.

---

## 2. Probit (binary) IFE with treatment reversals

### Motivation
`binary = TRUE` rejected any panel with treatment reversals
(`Binary IFE supports staggered adoption only`). fect 1.1.10 fitted the
Probit model on all untreated cells regardless of reversals. The empirical
script (democratisation and autocratisation episodes) needs this.

### Semantics
* The Probit fitter already works from the untreated-cell mask `II`, which is
  the same mask used by the linear estimators, so no estimation code changed.
* Switch-off effects (`att.off`, `time.off`, `est.att.off`) are now produced
  for binary fits exactly as for linear fits.
* Jackknife and nonparametric bootstrap work with reversals. Parametric
  bootstrap remains unavailable with reversals (pre-existing rule).
* Rolling CV: the "future tail" dropped from a held-out unit's training window
  now covers every later observed cell of that unit, not only later pre-onset
  cells. Under staggered adoption this is identical to the old behaviour;
  with reversals it also removes untreated cells observed after the unit
  switched off, so the training data never contain that unit's future.

### Files
| File | Change |
|---|---|
| `R/default.R` | Removed the `hasRevs && binary` stop. |
| `R/cv_binary.R` | Removed the `hasRevs` stop in `fect_binary_cv()`. |
| `R/cv-helpers.R` | `.build_cv_mask_rolling()`: future tail computed from all `II == 1` cells after the hold-out block. |
| `man/fect.Rd`, `man/r.cv.rolling.Rd`, `R/cv-rolling.R` | Wording updated (reversals allowed). |
| `NEWS.md` | Entry added; the "requires staggered adoption" sentence removed. |
| `tests/testthat/test-binary-rolling.R` | The test that asserted the old error now asserts a successful fit with switch-off output. |
| `tests/testthat/test-binary-reversal.R` | New: fixed-rank fits (r = 0, 1), rolling CV fold integrity with reversals, jackknife inference. |
| `R/cv-helpers.R`, `R/cv_binary.R` | Follow-up (2026-09-23): `cv.eligible.units` now comes from the fold builder, which counts only untreated cells before the first treated period. It used to count all untreated cells, so with reversals it could list units no fold can draw. Regression test in `test-binary-reversal.R`. |

### Verification
* All models of `binary_empirical.R` run: CV over r = 0..3 and the three
  Probit jackknife fits with `group = "htmil"`; gap plots render.
* Overall ATT versus stored fect 1.1.10 results: r = 2: -0.0218 vs -0.0216;
  r = 0: 0.0329 vs 0.0325; r = 1: 0.0129 vs -0.0067.
* Known difference, not introduced here: Probit covariate coefficients and
  their standard errors differ from fect 1.1.10 at every rank, including
  r = 0 (e.g. `logrgdppc` 0.31 vs 0.005). This comes from the branch's
  Probit fitter/initialiser (`BiInitialFit` + `inter_fe_d_*_ub`), which was
  not modified. Worth investigating before relying on the coefficient table.

---

## 3. PC information criterion for the Probit model

### Motivation
fect 1.1.10 reported, in its binary CV table, a Bai-Ng style criterion
(labelled "IC" there) that depends only on the full-sample fit, not on the CV
split. This branch's C++ `IC` slot is a BIC-type penalty with different values.

### Definition
```
PC(r) = r * (N + T) / (N * T) * log(N * T / (N + T)) - 2 * mean log-likelihood
```
where the mean log-likelihood is averaged over the cells used in the fit
(`loglikelihood` slot from C++). Computed in R by `.fect_binary_pc()`; the
C++ `IC` is left unchanged and still reported.

### Files
| File | Change |
|---|---|
| `R/cv_binary.R` | New helper `.fect_binary_pc()`. `CV.out` gains a `PC` column (columns: r, IC, PC, Log-likelihood, MSPE, MSPE.SE, Classification.error). One message line per rank prints PC, log-likelihood and MSPE, with `*` on the selected rank. |
| `R/fe.R` | Fixed-rank binary fits return `$PC`. |
| `R/cv-rolling.R` | `r.cv.rolling(binary = TRUE)` returns `pc` and `loglik` columns in `$mspe`. |
| `man/fect.Rd` | `CV.out` description expanded for the binary case. |
| `NEWS.md` | Entry added. |
| `tests/testthat/test-binary-pc.R` | New: formula check, independence from fold configuration, agreement between `CV.out` and fixed-rank fits, selected-rank `$PC`. |

### Verification
* Penalty term reproduces fect 1.1.10 exactly on the empirical panel
  (0.1342 per factor). PC levels differ only through the log-likelihood
  (see the fitter note in section 2).
* Rank selection still uses probability MSPE with the chosen `cv.rule`.
  Rolling folds are random unless `seed` is passed; PC is unaffected.

---

## Test status (run with `NOT_CRAN=true` under `pkgload::load_all(".")`)

| File | Result |
|---|---|
| `test-group-time-varying.R` | 20 pass |
| `test-binary-reversal.R` | 28 pass (24 before the 2026-09-23 `cv.eligible.units` test) |
| `test-binary-pc.R` | 10 pass |
| `test-binary-rolling.R` | 151 pass |
| `test-fect-basic.R` | pass |
| `test-group-fe.R` | 1 pre-existing error ("legacy index[3:] does NOT enforce nesting"), unrelated |

Note: `test-binary-rolling.R` calls internal helpers without `fect:::`, so it
only passes under a source load, not against an installed copy.

## Installation

```
cd /Users/liulch/Desktop/binary_replication
R CMD INSTALL --no-multiarch fect-codex-gsynth-binary-rolling-cv
```

Installed on 2026-09-21 into the system R 4.4-arm64 library.
