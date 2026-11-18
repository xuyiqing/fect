# fect: a guide for AI agents

This file is for an AI assistant (Claude, ChatGPT, Codex, Gemini and the like) that is helping someone analyze panel
data with the R package **fect**. Read it once, then work with the user as described below. It is written for fect
2.4.9. When it disagrees with the installed package, trust the package: check `packageVersion("fect")`, `args(fect)`
and `?fect`.

- User manual (Quarto book): <https://yiqingxu.org/packages/fect/>
- Source and issues: <https://github.com/xuyiqing/fect>
- Main paper: Liu, Wang and Xu (2024), "A Practical Guide to Counterfactual Estimators for Causal Inference with
  Time-Series Cross-Sectional Data", *American Journal of Political Science* 68(1): 160--76.

## What fect does

fect estimates the effect of a binary treatment in panel data (units observed over time). Every method works in the
same three steps:

1. Fit a model of the untreated outcome Y(0) on the untreated cells (unit-periods without treatment).
2. Predict Y(0) for the treated cells.
3. Average the differences between the observed outcomes and these predictions.

The result is the average treatment effect on the treated (ATT), overall and by period relative to treatment onset.
The methods differ only in the model for Y(0) and in which untreated cells fit it. Because fect never uses treated
cells to fit the model, it avoids the negative-weighting problem of two-way fixed effects regressions with staggered
treatment timing.

fect does not handle continuous treatments, and it assumes no feedback, meaning past outcomes must not drive later treatment.

## Talk with the user before writing code

Before writing any code, ask the user about the design. Ask only what you cannot see in their data.

1. What is the unit, what is the time period, what is the outcome, and what is the treatment?
2. Is the treatment coded 0/1 in every unit-period? Does it switch off for some units (reversal), or does it stay
   on once it starts?
3. Do all treated units start at the same time, or at different times (staggered adoption)?
4. How many units are ever treated, and how many are never treated? How many periods come before treatment?
5. Which covariates matter, and are they affected by the treatment? (Covariates affected by the treatment should
   not be controlled for.)
6. Should observations be weighted, and at what level should uncertainty be clustered?
7. What claim does the user want to make, and for which audience?

Then look at the data together. Plot the treatment pattern with the **panelView** package
(`panelview(Y ~ D, data = df, index = c("unit", "time"), type = "treat")`). Propose a short plan (method, inference,
diagnostics) and agree on it before running long jobs.

## Data requirements

- Long format: one row per unit and period. `index = c("unit_id", "time_id")` names the two columns.
- `D` must be 0/1. Missing `D` or index values stop the fit. Missing `Y` or `X` leave those cells out.
- Give the variables in a formula, `Y ~ D + X1 + X2`, with bare column names only. Terms such as `log(Y)`,
  `factor(g)` or `X1:X2` stop with an error, so create the column first. Covariates must be numeric. Turn a
  factor into dummy columns.
- A unit with fewer than `min.T0` untreated periods is dropped (5 by default for factor methods, 1 for `"fe"`), and
  listed in `fit$remove.id`. Tell the user when this happens.

## Choosing a method

| Data and concern | `method` |
|---|---|
| Many treated units, and parallel trends are plausible after unit and time fixed effects and covariates | `"fe"` (the default) |
| Pre-treatment estimates show trends, which suggests unobserved confounders that vary over time | `"ife"` (latent factors) or `"mc"` (matrix completion) |
| Few treated units (even one), treatment never switches off, many never-treated units | `"gsynth"` |
| The confounding has a known form, such as group-by-period shocks, time-invariant covariates with time-varying effects, unit-specific trends | `"cfe"` |

`"ife"` and `"gsynth"` fit the same model and differ only in `time.component.from`. `"ife"` estimates the factors
from all untreated cells, including treated units before treatment (`"notyettreated"`, the default). `"gsynth"`
estimates them from never-treated units only (`"nevertreated"`). When `r` (the number of factors) is not given,
cross-validation chooses it from 0 to 5. `"mc"` chooses its penalty `lambda` the same way.

A typical sequence:

1. Fit `"fe"` with `se = TRUE`. Look at the pre-treatment estimates in the gap plot and the tests in `fit$test.out`.
2. If the pre-treatment estimates move away from zero, fit `"ife"` and `"mc"`, and repeat the checks.
   `fect_mspe()` compares the prediction errors of fitted models.
3. With few treated units and no reversal, use `"gsynth"` with `vartype = "parametric"` or `vartype = "conformal"`.

## Choosing the uncertainty estimate

Set `se = TRUE` to get standard errors, intervals and p-values. `vartype` chooses how:

- `"bootstrap"` (the default) resamples whole units, or clusters with `cl = "column"`. It needs many treated units.
- `"jackknife"` leaves out one unit at a time. It also relies on the number of treated units.
- `"parametric"` simulates prediction errors from never-treated units. It is meant for few treated units. It needs
  `time.component.from = "nevertreated"` (as with `method = "gsynth"`), no treatment reversal, not `"mc"`, and at
  least two never-treated units.
- `"conformal"` (new in 2.4.9; check `packageVersion("fect")`) compares each treated unit with the never-treated units, each held out as if it
  were treated, and gives an interval for each treated unit and date with a guarantee that holds in finite
  samples. It works with `"gsynth"`, `"ife"` and `"mc"` and a treatment that never switches off. It gives no
  standard errors.

With one or very few treated units, do not report bootstrap or jackknife intervals as if they were reliable:
they resample the treated units, and there are too few to resample. Use `vartype = "parametric"` (with
`"gsynth"`) or `vartype = "conformal"`, and say what each rests on:

- The parametric bootstrap needs the factor model to hold for every unit, enough factors, and errors that are
  independent across units and alike for treated units and controls. Its intervals are approximate.
- Conformal intervals need the treated units to be no harder to predict than the never-treated units, and
  at least `1/alpha - 1` never-treated units (19 for `alpha = 0.05`). The tuning must be fixed in advance:
  give a single `r` (or `lambda`) with `CV = FALSE`, since a value chosen by cross-validation can break the
  guarantee. With several treated units, the default intervals for averages (`est.att`, `est.avg`) have no
  guarantee; add `conformal.block = TRUE` for the average effect. Add `conformal.diag = 4` for pre-treatment
  checks and `conformal.relabel = 999` for a test of no effect.

`nboots` is 200 by default. Use 1000 or more for intervals built from bootstrap quantiles (`ci.method = "basic"`,
`"percentile"`, `"bc"`, `"bca"` in `estimand()`). The manual's chapter on the
bootstrap and jackknife (<https://yiqingxu.org/packages/fect/07-inference>) explains the choices.

## A first analysis

```r
library(fect)
data(sim_base)
data(sim_gsynth)

fit <- fect(Y ~ D + X1 + X2, data = sim_base, index = c("id", "time"),
            method = "fe", force = "two-way", se = TRUE, nboots = 200,
            parallel = FALSE)
print(fit)               # overall ATT and its interval
fit$est.avg              # the same as a matrix
fit$est.att              # ATT by period relative to onset (1 = first treated period)
plot(fit)                # event-study ("gap") plot with intervals
plot(fit, type = "equiv")  # pre-treatment estimates against an equivalence range
fit$test.out             # F test and equivalence tests of the pre-treatment estimates
```

Few treated units:

```r
fit.gs <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id", "time"),
               method = "gsynth", r = c(0, 5), CV = TRUE, se = TRUE,
               vartype = "parametric", nboots = 500, parallel = FALSE)
fit.gs$r.cv                      # number of factors chosen by cross-validation
plot(fit.gs, type = "counterfactual")  # observed and predicted Y(0) for the treated
plot(fit.gs, type = "loading.overlap") # treated loadings against the controls'

fit.cf <- fect(Y ~ D + X1 + X2, data = sim_gsynth, index = c("id", "time"),
               method = "gsynth", force = "two-way", r = 2, CV = FALSE,
               se = TRUE, vartype = "conformal", alpha = 0.1,
               conformal.block = TRUE, conformal.diag = 4, seed = 1234)
print(fit.cf)                    # guarantee level and the block interval for the average
head(fit.cf$conformal$unit)      # an interval for each treated unit and date
head(fit.cf$conformal$diag$unit) # pre-treatment checks; intervals should cover zero
```

`r = 2` is fixed here for the example. In a real analysis, choose it before looking at the outcomes and report it.

## Diagnostics to run and report

- **Pre-treatment estimates.** The gap plot, the F test (`f.p`) and the equivalence tests (`f.equiv.p`,
  `tost.equiv.p`) in `fit$test.out`. Estimates far from zero before treatment are evidence against the model.
- **Placebo test.** `placeboTest = TRUE, placebo.period = c(-2, 0)` hides the last pre-treatment periods and
  estimates an effect there. It should be near zero. It needs a single `r` (or `lambda`), not cross-validation.
- **Leave-one-out pre-trends.** `loo = TRUE` refits once per pre-treatment period. `dloo = TRUE` (method `"fe"`,
  staggered adoption) gives double leave-one-out placebos.
- **Carryover test.** With treatment that switches off, `carryoverTest = TRUE, carryover.period = c(1, 3)` checks
  for effects that persist after the treatment ends.
- **Factor models.** `plot(fit, type = "factors")`, `type = "loadings"`, and `type = "loading.overlap"`. A treated
  unit whose loadings fall outside the controls' range is a warning sign, because its counterfactual is an extrapolation.
- **Sensitivity.** `fect_sens()` on a fit with a placebo test, then `plot(..., type = "sens")`, gives a Rambachan
  and Roth (2023) sensitivity analysis for violations of parallel trends. It needs the **HonestDiDFEct** package.
  The manual's Sensitivity Analysis chapter (<https://yiqingxu.org/packages/fect/10-sens>) shows the setup.

## Where the results are

| Slot | Contents |
|---|---|
| `att.avg`, `est.avg` | Overall ATT. With `se = TRUE`, `est.avg` adds `S.E.`, `CI.lower`, `CI.upper`, `p.value` |
| `est.att` | ATT by period relative to onset, with the same columns plus `count` |
| `eff`, `Y.ct` | Cell-level effects and predicted untreated outcomes (periods by units) |
| `test.out` | F test, equivalence tests, placebo or carryover test |
| `r.cv`, `lambda.cv` | Tuning chosen by cross-validation |
| `est.beta` | Covariate coefficients with inference |
| `remove.id` | Units dropped before fitting |
| `conformal` | With `vartype = "conformal"`: intervals for each treated unit and date (`unit`), the block, diagnostic and relabeling results, and the guarantee level |

`estimand(fit, "att.cumu", "event.time")` and its siblings give cumulative ATT, log-scale ATT and other targets.
The manual's Alternative Estimands chapter (<https://yiqingxu.org/packages/fect/03-estimands>) explains them.

## Common mistakes to avoid

- Do not invent arguments. Check `args(fect)` or `?fect` before using an option you have not seen work.
- `method = "gsynth"` stops when treatment switches off. `vartype = "parametric"` stops with
  `time.component.from = "notyettreated"`. `placeboTest` stops with cross-validation, so give a single `r`.
- For `"ife"` and `"gsynth"`, fect cross-validates when `r` is not given or is a range such as `c(0, 5)`. A
  single `r` fixes it. `"mc"` works the same way with `lambda`.
- `alpha = 0.05` gives 95% intervals in `est.att` and `est.avg`, and `est.att90` holds the 90% intervals used by the
  equivalence plot.
- `vartype = "conformal"` stops with weights, placebo and carryover tests, `loo`, `cm`, binary outcomes,
  `loading.bound` and treatment that switches off. `estimand()`, `effect()`, the equivalence plot and `fect_sens()`
  need bootstrap draws and stop on a conformal fit.
- `fect()` runs in parallel by default. On a shared machine pass `parallel = FALSE` or a small `cores`.
- Report the version (`packageVersion("fect")`), the method, `r` or `lambda`, the variance type and `nboots`.
- An ATT is causal only under the model's assumptions (parallel trends, or the factor structure, plus no
  anticipation and no feedback). Say which assumption carries the result and what the diagnostics show about it.

## Where to read more

| Topic | Manual chapter |
|---|---|
| Install, datasets | <https://yiqingxu.org/packages/fect/01-start> |
| The FE estimator, tests, weights | <https://yiqingxu.org/packages/fect/02-fect> |
| Other estimands | <https://yiqingxu.org/packages/fect/03-estimands> |
| IFE and MC, cross-validation | <https://yiqingxu.org/packages/fect/04-ife-mc> |
| Complex fixed effects | <https://yiqingxu.org/packages/fect/05-cfe> |
| Gsynth (few treated units) | <https://yiqingxu.org/packages/fect/06-gsynth> |
| Bootstrap and jackknife | <https://yiqingxu.org/packages/fect/07-inference> |
| Conformal inference (few treated units) | <https://yiqingxu.org/packages/fect/08-conformal> |
| Sensitivity analysis | <https://yiqingxu.org/packages/fect/10-sens> |
| Plots | <https://yiqingxu.org/packages/fect/11-plots> |
| Cheatsheet: every option and what fect changes or refuses | <https://yiqingxu.org/packages/fect/aa-cheatsheet> |
