<!-- markdownlint-disable MD025 -->
# fect 2.4.7

Development version, not yet on CRAN.

Version 2.4.6 was a development version on GitHub and was not released on
CRAN. Its changes are listed under fect 2.4.6 below and are part of this
release.

## Changes that affect results

Each bullet names the fits whose numbers change.

* `estimand()` now reads each bootstrap or parametric replicate through its
  own units (`fit$colnames.boot`). Its SEs and CIs for `"att"` with
  `by = "overall"`, `"aptt"`, `"log.att"` and the placebo and carryover series
  were wrong (3-6 times too small on `turnout`) unless the treated units were
  the first columns and adopted at one time (fect #141, #150; gsynth #37,
  #60).
* `imputed_outcomes(replicates = TRUE)` now returns each replicate's own
  treated cells: a unit drawn twice appears twice, so the row count varies by
  replicate, and the mean of `eff` within a replicate is its `att.avg.boot`
  draw. Before, every replicate reused the original cells.
* `effect()` now takes the variance type from the fit (`fit$vartype`), not
  from the call, and so does `estimand(fit, "att.cumu", "event.time")`;
  parametric fits whose `vartype` was not typed as a string (gsynth's
  default, or a variable) got CIs centred near 0 (gsynth #37, #75).
* `att.cumu()`, and so `estimand(fit, "att.cumu", "overall")`, now gives
  parametric fits a normal CI and p-value (estimate +/- z * SE); it took the
  quantiles of the parametric draws, which are centred near 0.
* The cumulative ATT (`att.cumu()`, `effect(cumu = TRUE)`,
  `estimand("att.cumu")`) is now the running sum of the per-period ATTs, as
  documented. It was k times the average effect over all treated cells in
  event times 1 to k, which differs when the number of treated units changes
  over event time (on `turnout`, `method = "gsynth"`, `r = 0`: 21.50 instead
  of 8.26 at k = 10). The `weighted` argument of `att.cumu()`, which had no
  effect, now chooses: `FALSE` (the new default) gives the running sum and
  `TRUE` the old estimate (for jackknife fits with a new S.E.; see the next
  bullet). `att.cumu()` and `effect()` now report the same S.E. (gsynth #75),
  and the same interval and p-value in every row (see the next bullet for
  the first row).
* The first row of `att.cumu()` (the window of one event time) now follows
  the rule of the other rows; it was copied from `fit$est.att`. For
  bootstrap fits it now has the percentile interval and p-value of the
  replicates, as `effect()` does (on `simgsynth`, `method = "ife"`, `r = 2`,
  `nboots = 200`, `seed = 1`: [-0.106, 2.921] with p = 0.10, instead of the
  normal interval [-0.322, 2.875] with p = 0.117). For parametric fits with
  `ci.method = "basic"` it now has the normal interval and p-value of the
  other rows (`method = "gsynth"` on the same data: [-0.143, 2.615] with
  p = 0.079, instead of the basic interval [-0.101, 2.535] with p = 0.08).
  Parametric and jackknife fits with the default `ci.method = "normal"` do
  not change (fect #153).
* For jackknife fits, `att.cumu()`, `effect()` and `estimand("att.cumu")`
  now report the jackknife S.E. of the cumulative ATT, computed as for the
  fit's own jackknife SEs, with a normal interval and p-value. `att.cumu()`
  used the spread of the leave-one-out estimates, about sqrt(N - 1) times
  too small for N units (on `turnout`, 47 units, `method = "gsynth"`,
  `r = 0`, at event time 10: 5.29 instead of 35.50; with `weighted = TRUE`,
  5.73 instead of 38.42), and took their quantiles as the interval.
  `effect()` was sqrt(N/(N - 1)) times too large (35.89, about 1% here) and
  used t critical values.
* For parametric fits, `effect()`, and so
  `estimand(fit, "att.cumu", "event.time")`, now gives the normal interval
  and p-value, as `att.cumu()` and the fit's own intervals do. It used a t
  critical value with `nboots - 1` degrees of freedom, so its intervals were
  wider (by 0.6% at the default `nboots = 200`, 2.5% at 50). The SEs do not
  change.
* For fits made with `W` or `W.agg`, `estimand()` and `effect()` now use
  those weights, in the estimates and in every replicate, as `?estimand`
  said and as `fit$att` and `fit$att.avg` do. So
  `estimand(fit, "att", "overall")` now equals `fit$att.avg`, with the SE of
  `fit$est.avg` (on `simgsynth` with unit weights, `method = "gsynth"`,
  `r = 2`: 4.6408 instead of the unweighted 4.6396). Fits made without
  weights, or with `W.est` only, do not change. The `weights` argument of
  `estimand()` must stay `NULL`; other values still stop, now with a message
  that points to `imputed_outcomes()`, which reports each treated cell's
  effect and weight.
* For jackknife fits, `estimand(fit, "att", "overall")` now takes its SE from
  each leave-one-out fit's overall ATT, as `fit$est.avg` does. It used each
  leave-one-out fit's average of the unit-level ATTs, a different quantity
  when units have different numbers of treated periods (on `turnout`: 3.84
  instead of 3.28).
* `fect_mspe()` now refits a fit made by `gsynth()` with `gsynth()`, so it
  scores the fit's own model. It refitted fect's default fixed-effects model
  and reported that model's score (on gsynth's `simdata`, with
  `force = "two-way"` and `seed = 1` in both `gsynth()` and `fect_mspe()`:
  an MSPE of 2.08 instead of the gsynth model's 1.85, both scored as in the
  next bullet), and it stopped on a fit made with a gsynth-only argument
  such as `inference`.
* `fect_mspe()` and `r.cv.rolling()` now score the model's full prediction
  at each held-out cell: the intercept, the unit and time effects, the
  covariates times their coefficients, the factors and the other CFE terms.
  Before, they read the refit's `Y.ct.full` there. That left out the
  covariates (a refit drops a row whose outcome is hidden and sets its
  covariates to 0), and for never-treated fits (`method = "gsynth"`, or
  `"ife"` or `"cfe"` with `time.component.from = "nevertreated"`) also the
  intercept and the fixed effects. Never-treated CFE fits predicted 0, so
  every such model got the same score. On `sim_gsynth` with `seed = 1234`,
  a correctly specified IFE model with `X1` and `X2` now scores 3.16 instead
  of 45.10, and the three never-treated models of the book's comparison of
  gsynth and CFE score 2.36, 2.31 and 8.87 instead of 100.96, 110.33 and
  110.33. The fits themselves do not change; the scores do, and
  `r.cv.rolling()` can choose a different number of factors (fect #163).
  * With the full prediction, a held-out unit with few periods before the
    held-out block adds large errors, because its loadings are fitted to
    those few periods. On `sim_gsynth` (`force = "two-way"`, `seed = 1`),
    `r.cv.rolling()` chooses r = 0 at the default `min.T0 = 5` and r = 2,
    the true number, with `min.T0 = 10`, for both `method = "gsynth"` and
    `method = "ife"`. A larger `min.T0` gives a steadier choice.
  * A held-out cell whose covariates are missing in the data is now left
    out of the score; it was scored with its covariates set to 0.
  * For CFE fits with extra fixed effects in `index`, the refit cannot place
    a held-out cell in its group, so the prediction there uses the effect of
    the group of hidden cells instead of the cell's own group. A
    never-treated CFE fit with `Z`, `Q` or extra fixed effects in `index`
    does not store those terms for treated units, so `fect_mspe()` does not
    score it at held-out cells of treated units (such cells are held out
    only when a not-yet-treated fit comes first in the list).
* `fect_mspe()` on fits made with `Y`, `D` and `X` column names (no
  formula) now looks up the data and the other arguments of the fit's call
  where `fect_mspe()` is called, as it does for formula fits. So it refits
  the fit's own model when the call passes an argument as a variable, such
  as `min.T0 = min.T0` or `W = W`. Before, such a variable took
  `fect_mspe()`'s own argument of the same name, so the refit could use
  another `min.T0` or leave out the weights (on `simgsynth`, a fit with
  `method = "ife"`, `r = 2` and `min.T0 = min.T0`, where `min.T0` is 7,
  scored with `fect_mspe(fit, k = 3, seed = 1)`: an MSPE of 2.48 instead of
  2.24; with `W = W` instead, 2.48 instead of 2.38; all four scored with
  the full prediction of the bullet above). A data frame
  local to a function was not found, and an argument such as
  `min.T0 = k + 3` could take `fect_mspe()`'s own `k`. A formula fit made
  inside a function can now be scored outside it.
* For parametric fits (`vartype = "parametric"`, gsynth's default) with
  `ci.method = "basic"`, the p-values of the effects (`est.att`, `est.avg`,
  the placebo test and the other effect slots) now compare the estimate with
  the draws, which are simulated with no effect: each p-value is twice the
  smaller of the shares of the centered draws at or above the estimate and
  at or below it. They compared the draws with zero, which gave values near
  1 whatever the estimate (on gsynth's `simdata`, `r = 2`: an ATT of 5.54
  with S.E. 0.25 had p = 0.97 and now has 0, as with `"normal"`; the placebo
  test's p-value is 0.10 instead of 0.90). The intervals and the coefficient
  p-values do not change, and neither do bootstrap fits or fits with
  `ci.method = "normal"`. The cohort effects of parametric fits with `group`
  (`est.group.att`) get these p-values under `ci.method = "basic"` (see the
  bullet on `est.group.att` below).
* The cohort effects of fits with `group` (`est.group.att`) now follow
  `ci.method`: the normal interval and p-value under `ci.method = "normal"`
  (the default) and the basic ones under `"basic"`. The two were swapped
  since v1.0.5, for bootstrap and parametric fits (on `simdata` with
  `group = 1 + id %% 2`, `method = "fe"`, `nboots = 200`, `seed = 1`:
  group 1 [6.045, 9.645] under the default, instead of the basic interval
  [6.051, 9.397]) (fect #155).
* For bootstrap fits with `ci.method = "basic"`, every p-value of the
  effects (`est.att`, `est.avg`, `est.avg.unit`, `est.att.off`, the placebo
  and carryover tests, the calendar, balanced, weighted, cohort and subgroup
  slots, and the `dloo` pre-trend estimates) and of the coefficients
  (`est.beta`) is now the p-value that goes with the basic interval: it is
  below alpha exactly when the (1 - alpha) basic interval excludes 0, at
  every alpha. It was the percentile p-value (twice the smaller share of
  the replicates on either side of zero), which goes with the percentile
  interval, so a basic interval could exclude 0 while p was above 0.05. On
  `simgsynth`, `method = "gsynth"`, `r = 2`, `vartype = "bootstrap"`,
  `nboots = 200`, `seed = 1`, the basic interval excludes 0 at event times
  -15, -13, -6 and 2, where p was 0.07, 0.08, 0.16 and 0.09 and is now
  0.040, 0, 0.044 and 0.044; 22 of the 30 p-values in `est.att` change.
  Only p-values change. Parametric fits keep their p-values, and so do fits
  with the default `ci.method = "normal"`, except in one slot: the
  switch-off effects of each subgroup (`est.group.output[[g]]$att.off`) now
  have the normal p-value, as every other slot under `"normal"`; they had
  the percentile p-value (on `simdata` with `group = 1 + id %% 2`,
  `method = "fe"`, `nboots = 200`, `seed = 1`: 0.041 instead of 0.035 at
  event time -21 in group 1) (fect #158).
* `att.avg.unit` (the "Tr units equally weighted" line of `print(fit)`) and
  `est.avg.unit` now average each treated unit's effects over its observed
  treated cells. A treated unit with a missing outcome in any period, before
  or after treatment, was left out, and the value was `NaN` when every
  treated unit had one. All methods change in the same way (on `simgsynth`,
  `method = "fe"`, with the period-5 row of unit 101 removed: 5.092 instead
  of 4.701; with that row removed for units 101 to 105: 5.021 instead of
  `NaN`, and with `nboots = 50`, `seed = 1`, an S.E. of 0.440 instead of
  `NA`) (fect #156).
* With `normalize = TRUE`, `fit$Y.dat`, `fit$est.cm` and the covariate
  array stored in the fit are now on the outcome's scale; they were divided
  by sd(Y). So `imputed_outcomes()`, `estimand("aptt")`,
  `estimand("log.att")`, `fect_iden()`, the raw and counterfactual plots and
  `plot(type = "hte")` of such fits now agree with `normalize = FALSE` (on
  `sim_base` with the outcome shifted by 30, `method = "fe"`, `cm = TRUE`:
  the APTT at event time 1 is 0.0102 instead of 0.0542, the `fect_iden()`
  statistic e0 is 6.52 instead of 5423.6, and the log ATT at event time 1
  is 0.0092 instead of an error). The ATT and the effects did not depend on
  `normalize` and do not change (fect #157).
* Parametric SEs with `normalize = TRUE` are no longer multiplied by sd(Y)
  (on `turnout`, 35.80 instead of 2.56; gsynth #14).
* The bootstrap now resamples groups of one unit correctly: with one treated
  unit every draw contains it (most draws lacked it and were dropped
  silently, e.g. 92 of 100), and with one reversal unit, or two never-treated
  units under `vartype = "parametric"`, the right units are drawn.
* `criterion = "pc"` now selects the `r` with the lowest PC for
  `method = "gsynth"` and for `method = "cfe"` with
  `time.component.from = "nevertreated"`, as it did for `"ife"`; these two
  used the MSPE rule (on `turnout`, r = 2 instead of 4; gsynth #24).
* `W.agg` alone no longer enters the model fit or the cross-validation scores
  of `method = "gsynth"` and never-treated `"ife"`/`"cfe"` fits with
  `CV = TRUE`, which now equal the unweighted fit except for the weighted
  averages (gsynth #101).
* `wgt.implied` now follows Xu (2017), `Lco (Lco'Lco)^-1 Ltr'`, whose columns
  rebuild the factor part of each treated counterfactual (the old formula used
  the treated units' Gram matrix). It is an Nco x Ntr matrix with unit ids as
  dimnames, also under `loading.bound = "simplex"`, whose weights were stored
  transposed (gsynth #17, #82).
* `loading.bound = "simplex"` is now applied with `CV = TRUE`, for
  `method = "gsynth"` without SEs, and in the parametric bootstrap replicates;
  it was dropped there silently.
* The `loo = TRUE` pre-trend refits now use the fit's `loading.bound` (with its
  `gamma.loading`), `time.component.from` and `para.error`; they used
  unbounded loadings, the not-yet-treated estimator for never-treated `"cfe"`
  fits, and `para.error = "auto"`.
* `method = "ife"` with `time.component.from = "nevertreated"` and `se = TRUE`
  now fits the never-treated model, as `se = FALSE` does, and reports
  `method = "gsynth"`; with `CV = FALSE` it fitted the not-yet-treated model
  (ATT 5.579 instead of 5.543 on `simgsynth`, `r = 2`).
* `method = "cfe"` with `time.component.from = "nevertreated"` now computes
  its bootstrap and jackknife SEs, and those of its `loo = TRUE` placebos,
  with the never-treated model, as it does the point estimate. They came
  from the not-yet-treated model, so SEs, CIs and p-values were off, mostly
  in the pre-treatment periods (on `simgsynth`, `r = 2`, jackknife: S.E. 0.46
  instead of 0.63 at event time -7).
* With `parallel = FALSE`, `se = TRUE` (or `permute = TRUE`) and `CV = TRUE`,
  `seed` now gives the same cross-validation folds, and so the same `r`, as
  `se = FALSE`; it used the folds of `seed + 1`.
* A factor time index is now used in its level order: levels that are
  increasing numbers give the same fit as the numbers, and levels that are
  numbers out of numeric order (as in `factor(as.character(1:30))`) are used
  in level order with a warning. It was sorted as text ("1", "10", "2", ...),
  which scrambled event time and made gsynth report treatment reversals
  (gsynth #13).
* A character time index that holds numbers is now read as numbers; it was
  also sorted as text.
* Covariates that are exact linear combinations of others, absorbed by the
  fixed effects, or constant on the cells used to fit the model are now
  dropped with one warning and an `NA` coefficient, so the estimates equal
  those of the fit without them. Collinear covariates changed the estimates
  silently (ATT 5.543 -> 6.510 on `simgsynth`) or crashed (gsynth #40, #80,
  #83).

### Calls that now stop

These calls ran, often with wrong numbers. They now stop with a message that
says what to do.

* `fect()` and `interFE()` formulas take bare column names only: `log(Y)`,
  `factor(g)`, `X1:X2`, `X1 * X2` and `I(X^2)` stop, naming the term (they
  were read as their bare variables, so `log(Y + 20) ~ D` fitted `Y ~ D`).
  Intercept terms such as `+ 0` or `- 1` are still accepted and ignored
  (gsynth #13, #23, #62).
* `X` given together with a formula stops (it was ignored), and so do
  duplicated names in `X` and an `X` that names the outcome or treatment.
* Non-numeric covariates (factor, character, Date) stop and ask for numeric
  dummy columns; logical covariates are used as 0/1, as before.
* A character time index must hold numbers; values such as `"t01"` stop.
* `cl` must be one column name, and for the cluster bootstrap (`se = TRUE`,
  `vartype = "bootstrap"`) it must be constant within each unit and have at
  least two clusters; a varying or single cluster gave NA or zero SEs
  (gsynth #41, #86).
* `method = "ife"` with `time.component.from = "nevertreated"` and `se = TRUE`
  stops where the never-treated model cannot be estimated, as `se = FALSE`
  does; it fitted the not-yet-treated model instead.
* `cm = TRUE` with `time.component.from = "nevertreated"` stops. It ran
  without `est.cm`, which that model does not estimate.
* `interFE()` stops, naming the covariate, when a covariate does not vary
  and the model has no fixed effects (`force = "none"`), and when a
  covariate is the sum of a unit-level and a period-level variable under
  `force = "two-way"`. Such a covariate got the coefficient `NaN` (about
  1e-16 on unbalanced panels) without a message.
* `binary = TRUE` in `fect()` and `interFE()` stops at once with a message
  that binary outcome (probit) models are not supported in this version. It
  failed later with unrelated errors such as "no applicable method for
  'predict' applied to an object of class "NULL"" or "non-numeric argument
  to mathematical function" (fect #160).
* `estimand(fit, "att.cumu", ..., cells = ...)` stops and points to
  `window`. With `by = "overall"` it ignored `cells` and returned the value
  over the full window (on `simgsynth` with the outcome shifted by 20,
  `method = "fe"`: 50.85 for any `cells`). The cumulative ATT sums the
  per-period ATTs over a range of event times, so use
  `window = c(L, R)` to choose the range (fect #161).
* `plot(type = "gap")` with an `id` that is a control unit or not in the
  data stops with a message naming the unit. So does `id` together with
  `loo`, `dloo` or `show.group`. Before, the gap plot ignored `id`, except
  that an id not in the data stopped with "Some specified units are not in
  the data." (fect #162).

* `method = "cfe"` now takes each unit's `Z` and `kappa` label from that
  unit's own rows, and each period's `Q` and `gamma` label from that
  period's own rows. It read `Z` and `kappa` from the first period and `Q`
  and `gamma` from the first unit of the filled panel, where an absent row
  holds 0, so one missing row in the first unit or the first period changed
  the model: with `Q.type = "linear"`, the trend was 0 in every period the
  first unit lacked; a unit without a first-period row had `Z = 0`
  throughout; units missing the first period shared one `kappa` group, and
  periods the first unit lacked shared one `gamma` group. Estimates on
  unbalanced panels with `Z`, `Q`, `gamma` or `kappa` change, and so can the
  CFE scores of `fect_mspe()` and `r.cv.rolling()` (on `sim_linear` with
  `Q.type = "linear"` and the period-20 row of unit 1 removed: an ATT of
  1.0062 instead of 1.0959, with 1.0064 on the full data; on `simdata`
  with `Z = "L1"` and the first-period row of unit 101 removed: 3.0740
  instead of 3.1322, with 3.0759 on the full data).
  Relabeling the units no longer changes any estimate. `Z` or `kappa` that
  vary within a unit, or `Q` or `gamma` that vary within a period, now stop
  with a message naming the unit or period; they were read from one row.
  Group labels (`gamma`, `kappa`) are now coded by their sorted unique
  values in R, so character, factor, non-integer and negative labels work
  on unbalanced panels too (they were cast to unsigned integers in C++);
  a missing label, or a non-numeric `Z` or `Q`, stops with a message.
  Balanced panels are unchanged (fect #168).
* With `normalize = TRUE`, the default equivalence-test threshold
  (`fit$tost.threshold`, `0.36 * sqrt(sigma2.fect)`) of `method = "gsynth"`
  fits and of `CV = TRUE` fits was too small by a factor sqrt(sd(Y)):
  `sigma2.fect` was multiplied by sd(Y) once instead of squared. The
  threshold, `test.out$tost.equiv.p` and `plot(type = "equiv")` change (on
  `sim_gsynth`, `method = "gsynth"`, `r = 2`: a threshold of 0.489 instead
  of 0.202, and an equivalence p-value of 0.954 instead of 0.998). And a
  `lambda` given to `method = "mc"` now means the same model with or
  without `normalize`: it is divided by sd(Y) before the fit, and
  `lambda.cv`, `lambda.seq` and `eigen.all` are reported on the outcome's
  scale (on `sim_gsynth` with `lambda = 0.01`: an ATT of 5.118 instead of
  5.085, the two-way fixed effects ATT, because the penalty removed the
  whole low-rank part). `fit$data.long`, which `panelview(fit)` draws, is
  now on the outcome's scale for every method, as are `est$VNT`, gsynth's
  `IC`, and the cross-validation tables and messages (each score by its own
  power of sd(Y)); those slots only report, so nothing else changes
  (fect #166).
* For parametric fits (`vartype = "parametric"`, gsynth's default) with
  `ci.method = "basic"`, every p-value of the effects and of the
  coefficients (`est.beta`) is now the p-value that goes with the basic
  interval, as for bootstrap fits since #158: below alpha exactly when the
  (1 - alpha) interval excludes 0. It was a counting rule that agreed with
  the interval only up to the resolution of the draws (2 / `nboots`), so a
  row could print a 95% interval that excludes 0 beside p = 0.05, and for
  the coefficients the percentile rule. P-values move by at most about one
  step (on `sim_gsynth`, `method = "gsynth"`, `r = 2`, `nboots = 200`,
  `seed = 11`: 0.042 instead of 0.05 at event time 1, where the 95%
  interval is [0.055, 2.709]). Intervals do not change (fect #169).

## New features

* `plot(fit, type = "gap", id = ...)` now draws the gaps of the chosen
  treated units, each unit's observed outcome minus its predicted untreated
  outcome (`fit$eff`). Before, it ignored `id` and drew the average over all
  treated units. With one unit, the plot shows the unit's gap in each period
  by time relative to its treatment onset, as gsynth 1.2.x did. With several
  units, it shows their unweighted average at each relative time. A control
  unit, or an id not in the data, stops with a message naming it. On a fit
  with parametric bootstrap draws (`vartype = "parametric"` and
  `keep.sims = TRUE`), the plot draws a band around the gaps, formed from
  their draws with the fit's `ci.method` and `alpha`, as gsynth 1.2.x did
  for one unit; fits with another variance type show point estimates only,
  with a message saying why. `?plot.fect` and the user manual explain that
  one unit's gap in one period is a noisy estimate of its effect, which is
  not identified. For the band, fits now record `ci.method` and `ci.alpha`.
  The plot shows no test statistics, and `return.test = TRUE` returns none
  for it. On `simgsynth` (`Y ~ D + X1 + X2`, `method = "ife"`, `r = 2`,
  `CV = FALSE`), the plot with `id = 101` now shows unit 101's gap in its
  first treated period, 0.338; before, it showed the average over the five
  treated units, 1.277 (fect #162; gsynth #106).

## Bug fixes

* Weights (`W`, `W.est`, `W.agg`) now work with `vartype = "parametric"`
  (fect #73, #150; gsynth #101).
* The cluster bootstrap now handles clusters of unequal size: `keep.sims =
  TRUE` no longer stops (the replicate arrays are padded with `NA`; see
  `fit$colnames.boot`), and no replicate is skipped in the counterfactual
  bands.
* `vartype = "parametric"` now warns that `cl` is ignored, and `print()` no
  longer shows "Cluster SE" for parametric and jackknife fits, which ignore
  `cl`.
* Bootstrap replicates that cannot be estimated are now dropped and counted in
  a message instead of stopping the run; a main fit whose factors are
  collinear stops with a clear message.
* `vartype = "parametric"` with fewer than two never-treated units stops with
  a clear message.
* Cross-validation of `r` for never-treated fits now searches at most
  `r = Nco - 1`, with a message, instead of crashing with
  `Mat::head_cols(): size out of bounds` (gsynth #32, #97).
* When cross-validation is skipped because only `r = 0` was given, or there is
  one never-treated unit, the message now says so instead of blaming too few
  pre-treatment records.
* fect's IFE cross-validation table under `criterion = "pc"` now labels its
  columns correctly.
* `effect()` now works for one unit (`id`) and for fits with one treated unit
  (gsynth #45, #53).
* `method = "cfe"` with `time.component.from = "nevertreated"` now runs with
  a single treated unit; it stopped with "'x' must be an array of at least
  two dimensions" (gsynth #45).
* `att.cumu()` now works on fits without standard errors; it stopped with
  "non-numeric argument to mathematical function".
* A `cells` formula in `estimand()` and `imputed_outcomes()` can now use
  variables of the function where it was written.
* `imputed_outcomes()` now reports the aggregation weights of fits made with
  `W` or `W.agg`, one row per treated cell; it listed every treated cell
  twice, with a `W.agg` column that was mostly `NA`. The fit stores these
  weights in `fit$W.agg`.
* `plot(type = "counterfactual")` now works without SEs, and
  `plot(type = "factors")` with a Date or character time index
  (gsynth #69, #84, #89).
* When dropping the periods without control observations leaves no treated
  observation, fect stops with a plain message instead of an unrelated error
  (gsynth #57).
* `fit$remove.id` now lists the removed units whenever units are removed, not
  only when the first unit is among them.
* `interFE()` now names the fixed effect that absorbs a covariate (the unit and
  time labels were swapped).
* `fect_mspe()` keeps working on fits made with `Y`, `D` and `X` given as
  column names: its refits pass the covariates only in the formula, so the
  new stop for `X` given together with a formula does not affect them.
* `fect_mspe()` now scores a fit whose data has the name of a function where
  `fect_mspe()` is called, such as `df` or `data`; it stopped with "object
  of type 'closure' is not subsettable". It skips a value of that name that
  is not a data frame, and when it finds no data frame, the message says so
  (fect #152).
* `att.cumu(fit, period = c(k, k))` and
  `estimand(fit, "att.cumu", "overall", window = c(k, k))`, a window of one
  event time, now work on fits with standard errors; they stopped with
  "subscript out of bounds" (fect #153).
* Jackknife fits made with `W` or `W.agg` and `placeboTest = TRUE` or
  `carryoverTest = TRUE` no longer stop when the model has covariates
  ("length of 'dimnames' [2] not equal to array extent"). Their
  `est.placebo` and `est.carryover` now have the 90% bounds, as unweighted
  fits do (7 columns; they had 5). On `simgsynth`, `method = "fe"`, unit
  weights `1 + id %% 3`, `placebo.period = c(-2, 0)`, no covariates: 90%
  bounds [-1.491, 4.875] (fect #154).
* `estimand()` on jackknife fits: with `ci.method` left at its default,
  `"att.cumu"`, `"aptt"` and `"log.att"` now use `"normal"`, the only
  method jackknife fits support, instead of stopping (fect #159).
* `estimand()` calls that stop now say why, without version wording: a `by`
  value that is not available (a column name, `"cohort"` or
  `"calendar.time"`, or `"overall"` for `"aptt"` and `"log.att"`) gets a
  message that names the values that work; `cells` or `window` with
  `"aptt"` or `"log.att"`, and the arguments that
  `estimand(fit, "att", "event.time")` does not take, are named; and
  `estimand(fit, "att", "overall")` without `keep.sims = TRUE` says to
  refit with it or to use `vartype = "none"` (fect #161).
* `cm = TRUE` now stores `est.cm` without standard errors and after
  cross-validation of the number of factors (as with `method = "ife"`
  without `r`), so `fect_iden()` and `plot(type = "hte", cm = TRUE)` work on
  those fits; they stopped, saying that `est.cm` was missing.
* `dloo = TRUE` works only with the default
  `time.component.from = "notyettreated"`; with `"nevertreated"`, `fect()`
  stops at once, with a message that names `time.component.from`.

## Documentation

* `?fect`: `seed` is needed for reproducible parallel bootstrap draws, since
  `set.seed()` does not fix them; `formula`, `X`, `index`, `cl`, `W.agg`,
  `criterion`, `loading.bound` and the returned `wgt.implied` are updated.
* `?plot.fect`: `id` applies to the gap, counterfactual and status plots
  only, and `nfactors` to the loadings plot only. The plot chapter of the
  user manual shows the gap plot of chosen units.
* `?fect` (`ci.method`, `normalize`, `binary`, the `att.avg.unit` and
  `est.group.att` values, and the names of `est.avg` and `est.avg.unit`,
  which it gave as `est.att.avg` and `est.att.avg.unit`), `?estimand` (`by`, `cells`, `window`,
  `direction`, `ci.method`), `?att.cumu`, `?fect_mspe`, `?r.cv.rolling`
  (`method`, `min.T0`) and `?interFE` (`binary`) describe the changes
  above. User manual: the comparison of gsynth and CFE in the gsynth
  chapter is rewritten for the new scores, the inference chapter describes
  the basic p-values and the other slots as they now are, and the
  estimands chapter lists which `by` values and filters each type takes.

# fect 2.4.6

Development version on GitHub, not released on CRAN. Its changes are part of 2.4.7.

Several of these also change results; each bullet says which calls.

* Add `dloo` and `dloo.adjust` flags to `fect()`: double (cohort-wise)
  leave-one-out pre-trend placebos, computed as a closed-form overlay on the
  in-sample fit -- **no re-fitting of the imputation model**. Like `loo`,
  `dloo = TRUE` fills the object's `pre.est.att` / `pre.att.bound` /
  `pre.att.boot` slots (and adds `dloo.test.out`), viewed with
  `plot(fit, dloo = TRUE)` and tested through the existing pre-trend machinery.
  The
  double-LOO removes both the attenuation and the staggered-adoption
  contamination bias of the in-sample placebo, and is algebraically identical
  (fixed effects cancel in the DiD) to re-fitting the imputation model on the
  restricted later-adopter control pool for every (cohort, period) -- see the
  Proposition in Li & Strezhnev. With covariates, `X beta` does not cancel in
  that DiD, so the covariate coefficients are re-estimated for each (cohort,
  period) on the same restricted pool (a closed-form two-way fixed-effects
  regression), and the identity with re-fitting holds with covariates too.
  An earlier development build used the full-sample coefficients, which are
  estimated from all untreated observations, including control units'
  observations after a cohort adopts; there, changing one never-treated
  unit's outcome in the last period moved the pre-treatment placebos.
  `dloo.adjust = TRUE` selects Liu (2025)'s
  pre-treatment-average baseline (equivalently, the double-LOO rescaled by
  `(g-2)/(g-1)` per cohort), which is preferable for benchmarking
  period-to-period *changes* in the parallel-trends violation. The flags
  require a staggered-adoption design (no reversal) and a balanced
  pre-treatment panel (post-treatment missingness is allowed); they error
  otherwise. Standard errors reuse `fect`'s own bootstrap: the linear overlay
  is applied to each case-resampled replicate inside the existing bootstrap
  loop (no separate resampler, and no per-replicate panels are retained), and
  the confidence intervals follow `ci.method` exactly as elsewhere in the
  package, so dloo inference matches the rest of `fect`. Only supported for
  `method = "fe"` (additive two-way fixed effects); other methods error. When
  combined with `group`, subgroup-wise placebo series are filled into
  `pre.est.group.output` (viewable via `plot(fit, dloo = TRUE, show.group = ...)`):
  each subgroup averages the same per-unit placebo over its own treated units
  with shared controls, so the subgroups aggregate back to the pooled series.
* `r` now defaults to `NULL`. For `method = "ife"` (or `"gsynth"`) a missing
  `r` means the number of factors is selected by cross-validation over 0 to 5,
  as `method = "mc"` already does for a missing `lambda`. Previously the
  default was `r = 0`, so `method = "ife"` without `r` silently fit the
  two-way fixed-effects model with no factors. With `CV = FALSE` and no `r`,
  the fit falls back to `r = 0` with a message, mirroring MC. An explicit
  `r = 0` is honoured silently. `fe` and `cfe` are unchanged.
* Calls that name the variables as strings, such as
  `fect(Y = "y", D = "d", X = "x", data = df, index = c("id", "time"))`, now
  use the same defaults as formula calls such as `fect(y ~ d + x, ...)`. Two
  defaults were different. First, `cv.method` was `"all_units"` instead of
  `"rolling"`, the documented default since v2.3.0. So when these calls
  cross-validated, they used block cross-validation and printed a deprecation
  note for `"all_units"`. They now use rolling cross-validation, which can
  select a different number of factors `r` (or `lambda`) than before. Set
  `cv.method = "block"` to reproduce earlier results. Second, `nlambda` was 0
  instead of 10, so any of these calls that cross-validated `lambda`, such as
  `method = "mc"` with no `lambda` given, stopped with
  `"nlambda" option misspecified.` They now run; with no `lambda` given they
  try 10 values, as formula calls do. Formula calls are unchanged.
* `seed` now fixes the cross-validation folds when `se = FALSE`. Before,
  `fect()` applied `seed` only to standard errors (`se = TRUE`) and
  permutation tests (`permute = TRUE`). Otherwise the folds were drawn from
  whatever state R's random number generator was in, and `seed` had no
  effect. So two identical calls, such as
  `fect(Y ~ D + X, data = df, index = c("id", "time"), method = "ife", CV = TRUE, r = c(0, 3), seed = 1)`,
  could report different cross-validation results and select a different
  number of factors `r` (or a different `lambda`). Now `seed = s` draws the
  same folds as calling `set.seed(s)` first. This is also what a call with
  `se = TRUE` does by default, so the same `seed` selects the same `r` with
  and without standard errors. Cross-validation results change for calls
  that set `seed` with `se = FALSE`. To reproduce an earlier result that
  relied on a `set.seed()` call before `fect()`, keep that call and drop
  `seed`. Calls without `seed`, and calls with `se = TRUE` or
  `permute = TRUE`, give the same results as before.
* Cross-validation now uses the settings passed to `fect()` in two cases
  where it silently used defaults. With `se = TRUE`, it ignored `cv.rule`,
  `cv.buffer`, `cv.donut`, `min.T0` and `proportion` and ran with `"1se"`,
  1, 1, 5 and 0, so it could select a different number of factors `r` (or a
  different `lambda`) than the same call with `se = FALSE`. And
  `method = "gsynth"`, like `method = "ife"` with
  `time.component.from = "nevertreated"`, always used `cv.rule = "1se"`,
  with or without `se`. Calls that set any of these arguments can now select
  a different `r`. Calls that leave them at their defaults select the same
  `r` as before, with one exception. With `se = TRUE`, the WMSPE, WGMSPE,
  Moment and GMoment columns of the cross-validation table now leave out
  pre-treatment periods with few treated units, as set by the default
  `proportion = 0.3`, as calls without `se` already did. So with
  `criterion = "moment"` the selected `r` can change.
* `method = "cfe"` with `time.component.from = "nevertreated"` now applies
  `cv.rule`. Its cross-validation kept its own rule whatever `cv.rule` was:
  it moved to a larger `r` only when that lowered the error by more than 1%
  relative to the best smaller `r`. With the default `cv.rule = "1se"` it can
  now select a smaller `r`, as the other methods do. `cv.rule = "1pct"`
  (the smallest `r` within 1% of the lowest error) usually gives the old
  selection. The cross-validation table itself is unchanged.
* `cv.donut` is now checked against `cv.nobs` for block cross-validation
  (`cv.method = "block"`). Block folds hold out `cv.nobs` consecutive periods
  and score only the middle `cv.nobs - 2 * cv.donut` of them, so a setting
  such as `cv.nobs = 3, cv.donut = 2` left nothing to score, and the
  cross-validation stopped with the internal message
  `No residuals to score.` It now stops before cross-validating and says
  which values work. `cv.donut` must also be a non-negative whole number.
  Rolling cross-validation, the default, does not use `cv.donut`.
* `effect(plot = TRUE)` now sets line widths with `linewidth`, so it no
  longer triggers ggplot2's warning that `size` for lines is deprecated
  (since ggplot2 3.4.0).
* Fix the equivalence (TOST) test for the leave-one-out pre-trend estimates
  (`loo = TRUE`): the reported p-value now uses the same pre-treatment window
  as the F test, the periods that pass the `proportion` cutoff. It previously
  ranged over every pre-treatment period, so an early period with a handful
  of treated units and a wide standard error could set the reported maximum.
  The in-sample test was not affected.
* User manual, Factor-Based Methods chapter: the LOO pre-trend example now
  fits the IFE model with `r = 2`. The previous call omitted `r` and ran with
  zero factors, so the figure captioned "IFEct" showed the two-way
  fixed-effects pre-trend. The text under the joint-test figures was
  rewritten to match. Thanks to Siyang Zhu (UPF) for the report.
* User manual, Inference chapter: the migration example under "Parametric
  bootstrap: valid regimes" now passes the data by name,
  `fect(data = data, Y = "Y", ...)`. Its three calls passed the data frame
  first without a name, so `fect()` took it as the formula and each call
  stopped with `argument "data" is missing, with no default`.

# fect 2.4.5

* Add `group.fe` to `fect()` for absorbing coarser fixed effects, such as state FE with county-level data. Closes #139. Clustered SE defaults to `group.fe[1]`; override with `cl = "<column>"`.
* Fix `method = "cfe"` with `force = "time"` or `"unit"`, which previously triggered an `Index out of bounds` error. `force = "two-way"` is byte-equivalent.
* Remove unused `sfe` argument and `R/polynomial.R`.

# fect 2.4.4

- `fect()` now returns `$sample`, a logical matrix (same dims as `$Y.dat`) marking cells used in any part of the estimation procedure (main fit, placebo/carryover/balance tests).

# fect 2.4.3

* Fix `future.globals.maxSize` overrun in parallel bootstrap: `quiet_nonpara`
  wrapper no longer captures `fect_boot()`'s full frame.
* Raise `future.globals.maxSize` to 2 GiB locally inside the parallel block.

# fect 2.4.2

## New: `ci.method` argument on `fect()`; legacy `quantile.CI` soft-deprecated

* `fect()` gains a `ci.method = c("normal", "basic")` argument. Default
  `"normal"` (Wald: `θ̂ ± z · SE`) preserves the v2.4.1 default behaviour
  byte-equivalently. `"basic"` (reflected pivot: `2 · θ̂ − quantile(boot, …)`)
  is the literature-standard "percentile" CI per @davison_hinkley1997 §5.2.1
  and what `boot::boot.ci(type = "basic")` returns. All CIs in fect's
  returned `est.*` slots use the requested method uniformly.
* The legacy `quantile.CI` argument is soft-deprecated. Both legacy values
  still work (`quantile.CI = FALSE` → `ci.method = "normal"`; `quantile.CI = TRUE`
  → `ci.method = "basic"`) but emit a one-time deprecation warning when
  user-supplied. Removal targeted for v2.5.0+.
* `ci.method = "basic"` with `nboots < 1000` emits a tail-CI replicate
  warning at fit time (mirrors the `estimand()` `.check_tail_ci_replicates`
  gate). The 5th / 195th order statistics that `basic` reads are unstable
  at small `B` --- @efron1987 §3 and @diciccio_efron1996 §4 recommend
  `B ≥ 1000` for tail-quantile CIs.
* `ci.method = "bca"`, `"bc"`, or `"percentile"` on `fect()` is rejected
  with a clear error pointing the user to `estimand(fit, type, ci.method)`
  for the full 5-method surface. fect's built-in CI machinery covers the
  routine `att` workflow; the alternative estimands (`att.cumu`, `aptt`,
  `log.att`) where bias-corrected CIs matter live on the `estimand()` path.

## New: alternative-estimand additions in `estimand()`

* New `test = c("none", "placebo", "carryover")` argument evaluates the
  requested estimand at pre-treatment placebo cells or early
  post-reversal carryover cells, producing a per-event-time series for
  credibility checks. Closes issue #131. Auto-pairs `direction = "on"`
  with placebo and `direction = "off"` with carryover; auto-validates
  the fit (placebo requires `placeboTest = TRUE` at fit time; carryover
  requires `carryoverTest = TRUE` + a reversal panel).
* `type = "att.cumu"` rejected with a clear error when `test != "none"` ---
  cumulative semantics are defined relative to treatment onset.
* `ci.method` enum extended from `c("basic", "percentile")` to
  `c("basic", "percentile", "bc", "bca", "normal")`. New methods:
  - `"bc"` --- bias-corrected percentile (Efron 1987 minus acceleration)
  - `"bca"` --- bias-corrected accelerated (Efron 1987 in full); cell-level
    jackknife computes the acceleration with no extra refits
  - `"normal"` --- Wald CI `θ̂ ± z · SE`
* `ci.method` now defaults to `NULL`, which triggers a per-type default:
  - `"att"` → `"normal"`
  - `"att.cumu"` → `"basic"` (reflected pivot CI; matches Davison-Hinkley 1997 §5.2.1 and `boot::boot.ci(type = "basic")`)
  - `"aptt"` → `"bca"`
  - `"log.att"` → `"bca"`
  Existing scripts that pass `ci.method` explicitly are unaffected.

## New: `para.error` argument for `vartype = "parametric"`

* `fect()` gains a `para.error = c("auto", "ar", "empirical", "wild")`
  argument selecting the residual-error model the parametric bootstrap
  draws from. Replaces the implicit panel-shape-driven dispatch.
* `"auto"` (default) resolves at fit time and stores the resolved label
  on `fit$para.error`:
  - `"empirical"` on a fully-observed panel
  - `"ar"` on a panel with missing cells
* `"ar"` --- the v2.4.1 behavior: AR(1) error process estimated from
  control residuals. Works on any panel shape.
* `"empirical"` --- i.i.d. column-resample from the main-fit residual
  pool. Requires a fully-observed panel.
* `"wild"` --- Liu 1988 / Mammen 1993 / Cameron-Gelbach-Miller 2008
  unit-level Rademacher sign-flips over the empirical residual pool.
  Requires a fully-observed panel; preserves within-unit dependence.
* `para.error` is silently ignored when `vartype != "parametric"`.

## Changed: tighter EM convergence defaults

* `tol`: default flipped from `1e-3` to **`1e-5`**.
* `max.iteration`: default flipped from `1000` to **`5000`**.

The pre-v2.4.2 default `tol = 1e-3` halted IFE/CFE EM well before
convergence: on factor-DGP simdata the EM stopped at iteration 116
with `att.avg = 2.87`, while running to `tol = 1e-7` (~2000 iters)
produces `att.avg = 2.43` --- an 18% gap between two valid stopping
points of the same procedure on the same data. CFE was worse
(40% gap). The new defaults stop the EM after it has actually
stabilized.

* **Inference at the old default was already correct.** Coverage
  simulations (K=80, known-truth DGP, true τ=3) show empirical
  coverage of 0.96 at both old and new defaults for IFE; bootstrap
  SE matches empirical SE in both. The fix improves
  *reproducibility* and *point-estimate stability across
  versions/machines*, not coverage validity.
* **What this means for users**: numerical output from prior versions
  remains valid inferentially (CIs still cover correctly), but the
  point-estimate values will shift on rerun under v2.4.2 --- typically
  by a few percent on canonical IFE, up to 40% on factor-heavy CFE.
  The new numbers are closer to the EM's actual converged minimum.
* **Speed cost**: ~2-5x slower main fit and bootstrap on
  factor-DGP IFE/CFE because EM iterates more (994 iters at 1e-5
  vs 116 at 1e-3 on simdata). GSC and MC paths unaffected
  (they were already converging within the old defaults).
* New `warning()` when EM hits `max.iteration` without satisfying
  the tol gate --- alerts users to under-converged fits on hard
  cases (e.g., very large N panels, near-collinear factors).

Set `tol = 1e-3, max.iteration = 1000` explicitly to reproduce
pre-v2.4.2 numerical output exactly.

## Bug fixes

* `vartype = "parametric"` × `ci.method ∈ {"basic", "percentile", "bc",
  "bca"}` produced 0% coverage CIs through v2.4.1 (the bootstrap
  distribution is H₀-centered, but reflection-based CIs assume centering
  at θ̂). `estimand()` now applies a variance-preserving location shift
  for parametric fits. The `"normal"` ci.method is byte-stable; the
  other four now produce nominal coverage.
* `vartype = "jackknife"` was previously rejected by `estimand()` with a
  slot-contract error. The slot contract is relaxed; only `ci.method =
  "normal"` is accepted (the Wald-style CI from the Tukey SE), with
  hard-error guidance pointing at `"bootstrap"` for the full ci.method
  surface.
* `log.att` and `aptt` silently dropped bootstrap replicates with
  `Y0_b ≤ 0` via `colMeans(..., na.rm = TRUE)`, contaminating the
  bootstrap distribution. Both now hard-error with actionable guidance
  (pre-transform Y, filter near-zero cells, or use a different
  estimand). `estimand("log.att", ...)` additionally hard-errors at the
  point-estimate level when any treated cell has `Y_obs ≤ 0` or
  `Y0_hat ≤ 0`.
* `vartype = "parametric"` with default `time.component.from =
  "notyettreated"` now produces a clearer error that names the user's
  literal `method` argument (was: "Parametric bootstrap is not valid
  when ..."; now: "vartype = 'parametric' requires time.component.from
  = 'nevertreated'. Your call: method = 'fe', time.component.from =
  'notyettreated'."). The reversal-check gate continues to fire first
  on reversal panels.
* `R/diagtest.R` "F-test Failed" message → "F-test could not be
  computed" --- the test never "failed" in any standard sense; the
  matrix arithmetic was undefined for the input.
* Parallel-worker package version warnings (e.g. "package 'mvtnorm'
  was built under R version X.Y.Z") suppressed via `clusterEvalQ`
  pre-load + targeted `withCallingHandlers`.

## Other changes

* `estimand()` warns when `ci.method` `c("basic", "percentile", "bc",
  "bca")` is requested on a fit with fewer than 1000 bootstrap
  replicates, recommending refit at `nboots = 1000` for stable tail
  quantiles (Efron 1987 §3; DiCiccio & Efron 1996 §4). The point
  estimate and SE are unaffected; the warning fires on every such
  call so the user can decide whether to suppress, refit, or
  proceed with caveat.
* `vartype = "jackknife"` with `Nco > 1000` emits a fit-time warning
  recommending `vartype = "bootstrap"` for tractability (full
  leave-one-out scales linearly in N and is slow at the v2.4.2 EM
  convergence defaults).
* `complex_fe_ub` and `cfe_iter` C++ entries gain optional `fit_init`
  parameter (NULL default preserves pre-existing cold-start behavior).
  This mirrors the existing warm-start infrastructure on
  `inter_fe_ub` / `inter_fe_mc` / inner EM helpers. Not exposed to
  the public API; reserved for future deferred features.

# fect 2.4.1

## New: parametric variance support in `estimand()`

* `estimand()`'s `vartype` argument now accepts `"parametric"` in
  addition to `"bootstrap"`, `"jackknife"`, and `"none"`. When the
  fit was produced with `fect(..., vartype = "parametric",
  keep.sims = TRUE)`, all four `type` values (`"att"`, `"att.cumu"`,
  `"aptt"`, `"log.att"`) work without modification, sourcing
  replicates from the parametric `fit$eff.boot` surface populated
  by the existing fit-time machinery.
* Byte-equality between `estimand(fit, "att", "event.time")` and
  `fit$est.att` is preserved under parametric (asserted by tests).
* The output `vartype` column reports the variance method actually
  used at fit time (read from `fit$vartype`), which may differ from
  the user-supplied `vartype` argument; the argument is informational
  and does not re-aggregate replicates.
* No changes to fit-time machinery, slot semantics, or the v2.4.0
  contract documented in `statsclaw-workspace/fect/ref/po-estimands-contract.md`.

# fect 2.4.0

## New: post-hoc estimands API

* New `estimand(fit, type, by, ...)` typed dispatcher computes
  alternative estimands directly from any fect imputation fit. Shipped
  types: `"att"` (default per-event-time ATT, byte-identical to
  `fit$est.att`), `"att.cumu"` (cumulative ATT, replaces `effect()`),
  `"aptt"` (average proportional treatment effect on the treated;
  Chen & Roth 2024), `"log.att"` (mean log-scale treatment effect).
  Shipped `by` axes: `"event.time"`, `"overall"`, plus reserved
  canonical values `"cohort"` / `"calendar.time"` for future commits.
  Returns a tidy data frame with consistent columns
  `(<by_key>, estimate, se, ci.lo, ci.hi, n_cells, vartype)` regardless
  of `type`. See `?estimand` and the new "Alternative estimands"
  vignette chapter for worked examples and the design rationale.
* New `imputed_outcomes(fit, cells, replicates, direction)` low-level
  accessor returns the cell-level imputed potential-outcome surface as
  a long-form data frame with documented columns
  `(id, time, event.time, cohort, treated, Y_obs, Y0_hat, eff,
  eff_debias, W.agg, [replicate])`. Use this for custom estimands the
  dispatcher does not ship; pipe to dplyr / data.table for arbitrary
  aggregation.
* `cells = ` filter argument on both `imputed_outcomes()` and
  `estimand()` accepts NULL (default), a logical vector, or a one-sided
  formula evaluated against the long-form data
  (e.g. `~ event.time %in% 1:5 & !id %in% bad_ids`). The
  `window = c(L, R)` argument on `estimand()` is sugar over
  `cells = ~ event.time >= L & event.time <= R`.
* `direction = c("on", "off")` on both functions selects the event-time
  grid for reversal panels.
* New `eff_debias` reserved slot on the fit object (NULL for plain
  imputation estimators; populated by future doubly-robust estimators)
  so DR scores can be added to the surface without breaking the
  long-form schema.

## Soft-deprecation: `effect()` and `att.cumu()`

* `effect()` and `att.cumu()` continue to work byte-identically to
  v2.3.x and emit a one-time-per-session message pointing at the
  unified `estimand()` API. Removal not before v3.0.0. Migration:
  - `effect(fit, cumu = TRUE)` → `estimand(fit, "att.cumu", "event.time")`
  - `effect(fit, cumu = FALSE)` → `estimand(fit, "att", "event.time")`
  - `att.cumu(fit, period = c(L, R))` →
    `estimand(fit, "att.cumu", "overall", window = c(L, R))`
  Numerical equality is asserted by package tests.

# fect 2.3.3

## Bug fixes

* Fix `"incorrect number of dimensions"` crash in `diagtest()` that
  surfaced intermittently with `parallel = TRUE` + small `nboots`.
  Two layers: (a) `R/diagtest.R` now uses `drop = FALSE` when filtering
  bootstrap columns by all-non-NA, so a single surviving column stays
  a matrix; (b) the bootstrap parallel path in `R/boot.R` now builds
  the PSOCK cluster via `parallelly::makeClusterPSOCK(rscript_libs =
  .libPaths())` (the same robust pattern used by the CV path),
  wrapped in a 3-attempt retry-with-backoff. If `doParallel` cluster
  init exhausts retries the bootstrap now degrades to sequential
  rather than crashing.
* `carryover.rm` is now stored on the fit object. `plot.fect()` no
  longer reads it from `as.list(x$call)$carryover.rm` (which silently
  gave the wrong K under `do.call()`, programmatic wrappers, or any
  call-rewriting code path). Behaviorally a no-op for fits built via
  named-argument `fect(...)`; correct under any other construction.

## Documentation

* Chapter 2 §Other estimands gains a worked example for post-hoc
  estimands derived from the imputed potential-outcome surface,
  showing APTT (Chen & Roth 2024) with bootstrap CIs from the
  existing fit slots. Issue #126 (ajunquera).

# fect 2.3.2

## Modern visual defaults for `plot.fect()` (visual breaking change)

* Default visual overhaul across all 14 plot types under
  `theme.bw = TRUE`: white panel, plain left-aligned title, thin grey
  reference lines, dashed treatment-onset vline, pre/post lightness
  contrast (`grey50` / `grey20`), publication-sized axis text
  (`cex.main = 11`, `cex.lab = 9`, `cex.axis = 8`, `cex.text = 3.0`),
  and compact legends.
* Placebo / carryover plots render highlighted periods as a single
  accent glyph (orange triangle for placebo, blue diamond for
  carryover, orange triangle for `carryover.rm`) instead of a stacked
  pair of circle + accent. Background rectangle behind each highlight
  period is opt-in via `highlight.fill = TRUE` --- default is glyph
  only, which keeps figures clean for print and grayscale.
* `highlight` argument extended to accept a character subset of
  `c("placebo", "carryover", "carryover.rm")` for selective
  per-test-type highlighting (e.g., `highlight = "placebo"` to
  render carryover periods as plain circles when both tests ran at
  fit time). Backward-compatible: `NULL` / `TRUE` / `FALSE` still
  behave as before.
* Stats annotation block (placebo / carryover / F / equivalence
  p-values) sits at the top-left panel corner with symmetric 2.5%
  inset and `2 × num_stat_lines` top padding so it never grazes the
  leftmost CI. Sized at `cex.text * 1.0` so it does not overpower
  the title. User-supplied `stats.pos` still wins.
* `loadings` (ggpairs) plot: correlation panel reformatted with
  overall + per-group entries; per-group label colors match the
  density-plot fills.
* Migrated off ggplot2 4.0's deprecated `fatten` / `lwd` arguments
  to the `size` / `linewidth` aesthetics. Clears the per-plot
  deprecation warnings.

## `legacy.style = TRUE` escape hatch

New `legacy.style` argument (default `FALSE`). Pass
`legacy.style = TRUE` for byte-identical reproduction of pre-2.3.1
figures (bold centered title, larger axis sizes, solid vline, blue
placebo triangles, no peach rectangle), regardless of `theme.bw`.

## `theme.bw = FALSE` soft-deprecated

Setting `theme.bw = FALSE` now emits a one-time per-session message
flagging removal in v2.5.0. Users who want the gray-panel look
should pass `legacy.style = TRUE` (which honors `theme.bw = FALSE`
exactly), or apply `+ ggplot2::theme_gray()` to the returned plot.

# fect 2.3.1

## New: `W.est` and `W.agg` arguments distinguish survey weights from IPW / balancing weights

`fect()` and `fect.formula()` gain two new arguments that control where
the weight column enters the estimator:

* `W.est` --- weight column for the outcome-model fit (the weighted least
  squares applied inside the IFE / MC / CFE solver).
* `W.agg` --- weight column for the across-treated-obs aggregation
  (`att.on`, `est.avg`, `est.att`).

Both default to `NULL` and fall back to the existing `W` argument when
left unset, so callers who pass only `W = "col"` get the same behavior as
v2.3.0 + the consistency fix below (W enters both fit and aggregation).
Pass the per-role arguments to specify finer behavior:

* Survey / sample weights: `W = "ws"` (or equivalently
  `W.est = W.agg = "ws"`). W enters both fit and aggregation.
* Robust-regression / heteroskedasticity / GLS weights: `W.est = "wr"`
  alone. The fit is weighted, the aggregation is unweighted.
* Inverse-probability / balancing / post-stratification weights:
  `W.agg = "ipw"` alone. The outcome model is fit unweighted (preserving
  doubly-robust properties), and the aggregation is weighted by IPW.

In v2.3.1, `W.est` and `W.agg` (when both supplied) must point to the
same column. Truly distinct columns for fit vs. aggregation (e.g. a
combined survey x IPW design where the outcome model uses survey weights
and the aggregation uses survey x IPW) are scheduled for v2.4.0; the
v2.3.1 design errors with an instructive message if requested.

**Caveat for IPW users.** `W.agg = "ipw"` fits the outcome model
unweighted and applies IPW only at the across-treated-obs aggregation.
This is closer to a doubly-robust estimator than v2.3.0's silent
everywhere-weighting --- but it is not a fully cross-fit doubly-robust
estimator. Residuals on never-treated controls used in any de-bias term
inherit in-sample shrinkage from the outcome fit, which DR theory
requires cross-fitting to eliminate. A fully cross-fit DR path is
scheduled for v3.0.

## Breaking change: weighted fits now have a single, consistent ATT surface

When `W = "<column>"` is supplied to `fect()`, every reported quantity on
the returned fit object now reflects those weights. Prior versions
maintained two parallel ATT pipelines on weighted fits: an unweighted one
(populated into `est.att`, `est.avg`, `att.boot`, `att.vcov`, etc.) and a
W-weighted one (populated into the parallel `est.att.W`, `est.avg.W`,
`att.W.boot`, `att.W.vcov` slots). `plot(fit)` silently substituted the
W-weighted pipeline for rendering while `print(fit)` and `fit$est.att`
returned the unweighted pipeline --- so the same fit object reported
different per-period ATTs and aggregate CIs depending on which surface
the user looked at.

As of 2.3.1, when `W` is non-NULL:

* `fit$est.att`, `fit$est.avg`, `fit$est.att90`, `fit$att.bound`,
  `fit$att.boot`, `fit$att.vcov`, `fit$est.placebo`, `fit$est.carryover`,
  and the `*.off` reverse-treatment counterparts all carry the
  W-weighted aggregations.
* `fit$att`, `fit$time`, `fit$count`, `fit$att.avg`, `fit$att.off`,
  `fit$time.off`, `fit$count.off`, `fit$att.placebo`, `fit$att.carryover`
  similarly carry the W-weighted aggregations from the per-method
  estimator.
* `print(fit)` now labels the obs-level row as `Tr obs sample-weighted (W)`
  (instead of `Tr obs equally weighted`) when W was supplied.
* The redundant `*.W` slots (`est.att.W`, `est.avg.W`, `att.W.boot`,
  `att.W.vcov`, `att.W.bound`, `att.on.W`, `time.on.W`, `count.on.W`,
  `att.avg.W`, `att.on.sum.W`, `W.on.sum`, `att.off.W`, `time.off.W`,
  `count.off.W`, `att.off.sum.W`, `W.off.sum`, `att.placebo.W`,
  `att.carryover.W`, `est.placebo.W`, `est.carryover.W`, `est.att.off.W`,
  `att.off.W.bound`, `att.off.W.vcov`) are no longer present on the
  returned fit object.

If you want the unweighted view of the same fit, refit with `W = NULL`.

The `weight` argument to `plot.fect()` is now a no-op (deprecated),
slated for removal in v2.5.0; passing it emits a deprecation warning.
Internally `plot.fect()` no longer auto-flips between two pipelines ---
it reads the canonical slots, which are already W-weighted when W was
supplied at fit time.

The C++ matrix-completion / IFE / CFE solvers (`inter_fe_mc`,
`inter_fe_ub`, `inter_fe_cfe`) already used W as a fit-time
weighted-least-squares weight in 2.3.0; this release does not change the
fit, only the result-object surface.

# fect 2.3.0

## Rolling-window cross-validation (standard ML design)

* New exported function `r.cv.rolling()`: a standalone user-facing helper
  for picking the number of factors `r` via standard rolling-window CV.
  For each of `k` folds, a fraction `cv.prop` of eligible units (controls
  plus treated pre-treatment) is sampled; only sampled units carry a
  mask in that fold. For each sampled unit, a random anchor time `t*` is
  drawn and the fold's training set excludes:
    1. `cv.nobs` observations starting at `t*` (the held-out, scored block);
    2. `cv.buffer` observations immediately before `t*` (gap buffer, attenuates
       AR-leakage at the past-side train/test boundary --- analogous to
       `cv.donut` for the existing CV strategies, but only on the past side
       since the future side is dropped by construction);
    3. all observations from `t* + cv.nobs` through the unit's
       end-of-eligible (the rolling-window step --- training cannot see
       the future of the held-out block; for treated units, end-of-eligible
       is the cell strictly before treatment onset, so post-treatment
       cells are never masked).
  MSPE is scored at the held-out block only and averaged across folds.
  This is the standard time-series CV design (cf. `forecast::tsCV`,
  `tidymodels::sliding_window`, `caret::createTimeSlices`) adapted to
  panel data.

* Per-fold unit sampling is required: masking every eligible unit at the
  same time would leave no donor data at the masked time points and break
  factor identification at the masked tails. With `cv.prop = 0.2`, every
  eligible unit lands in the holdout roughly `k * cv.prop = 2` times in
  expectation across `k` folds; unsampled units stay fully observed and
  contribute training data at every period.

* New parameters: `cv.buffer` (default 1, past-side buffer length ---
  analogous to `cv.donut` for the existing CV strategies, but only on
  the past side because the future side is dropped by construction);
  `k` (default 10 folds, matching the default for the existing CV
  strategies); `cv.prop` (default 0.2, fraction of eligible units
  sampled per fold; raised from 0.1 after small-panel stability tuning);
  `seed` (optional integer base seed for reproducible per-fold sampling
  and anchor selection).

* Closes the forward-leakage channel that the existing
  `cv.method = "all_units"` / `"treated_units"` (random contiguous-block
  masking) leaves open at `cv.donut = 0 / 1` under serially correlated
  residuals: under rolling-window CV, the train/test boundary on the
  future side is closed by construction.

* Workflow: call
  `r.cv.rolling(formula, data, index, method = "ife", cv.buffer = 1, k = 10, cv.prop = 0.2)`
  to get the chosen `r.cv`, then pass that to
  `fect(..., CV = FALSE, r = r.cv, se = TRUE)` for the inferential fit.

* Identical CV behavior across `method = "ife"` (IFE-EM, internal
  `time.component.from = "notyettreated"`) and `method = "gsynth"` (GSC,
  internal `time.component.from = "nevertreated"`). Both paths populate
  `Y.ct.full` at masked positions: the IFE-EM path via EM imputation,
  the GSC path via the model-implied factor product `F * t(lambda_co)`
  (see "GSC: Y.ct.full populated at control positions" below). Other
  methods (e.g. `"mc"`) are not yet supported.

* Return value is a list with `r.cv`, `cv.rule`, `mspe` (data.frame of
  per-r MSPE averaged across folds, plus fold-SE and held-out cell
  counts), `mspe.per.fold` (r-by-k matrix of per-fold MSPE), and the
  chosen `k`, `cv.nobs`, `cv.buffer`, `cv.prop`. Promoting the design
  to a `cv.method = "rolling"` option inside the main `fect()` CV
  dispatcher is deferred to a future release.

* Empirical motivation: on the Eibl & Hertog (2023) oil-rich panels with
  residual AR(1) of 0.56--0.93, fect's default CV (random anchors,
  `cv.donut` 0 or 1) pegs `r.cv = 5` on every (cell, estimator, rule)
  combination, even at widened `cv.nobs = 6` --- the forward-leakage hides
  the rank overfit. `r.cv.rolling()` recovers ranks consistent with the
  placebo-passing preferred rank for each outcome.

* **Behavior change vs the v2.3.0 development tip's tail-only design**:
  the prior implementation deterministically masked the LAST `cv.nobs`
  observations of each control unit (no folds, no random anchors).
  This was simpler but gave a single MSPE estimate per `r` with no
  fold-to-fold SE. Existing callers' `r.cv.rolling()` invocations will
  produce different numerical results under the new design; the
  selected `r.cv` is typically similar but no longer deterministic for
  a fixed dataset (set `seed` for reproducibility).

## CV API: cv.method = "rolling" in the main dispatcher (additive)

* `fect(CV = TRUE, cv.method = "rolling", ...)` now wires rolling-window
  cross-validation into the main `fect()` CV dispatcher for all factor-model
  methods (ife, cfe, gsynth, mc, both). The same masking logic that powers
  the standalone `r.cv.rolling()` is now available through the standard
  CV API: per-fold sampling of `cv.prop` of eligible units (controls plus
  treated pre-treatment), random anchors per sampled unit, `cv.nobs`-cell
  scored holdout, `cv.buffer`-cell past-side buffer, drop-from-anchor-to-
  end-of-eligible (rolling-window step). Treated post-treatment cells are
  never masked.

* New parameter `cv.buffer` (default 1) controls the past-side buffer for
  rolling CV; replaces the role `cv.donut` plays for block CV. `cv.donut`
  is unchanged for block strategies. Defaults are otherwise unchanged ---
  the existing `cv.method = "all_units"` / `"treated_units"` defaults
  still apply, so this is a fully additive change with no breaking
  defaults.

* `fect_mspe(out, cv.method = "rolling", cv.buffer = 1, ...)` for
  rolling-window-CV-based model comparison.

* `r.cv.rolling(method = "cfe", ...)` extends the standalone helper to
  Complex Fixed Effects. CFE-specific args (Z, gamma, Q, Q.type, kappa,
  extra index columns) are forwarded via `...` and held fixed; rolling
  CV picks `r` only.

* Default `cv.method` flips ("rolling" as the package-wide default),
  `cv.donut → cv.buffer` rename, and `(cv.method, cv.units)` API
  decomposition all remain deferred to a future "CV API unification"
  PR per the plan at statsclaw-workspace/fect/ref/cv-unification-plan.md.

## Bounded factor loadings for GSC

* New argument `loading.bound = "simplex"` (default `"none"`): constrains
  treated-unit factor loadings to the convex hull of control loadings via an
  entropy-regularized simplex projection. Solves, per treated unit `i`:

  ```text
  minimize     (1/gamma) * KL(w || uniform) + || u_pre - F_pre %*% t(Lambda_co) %*% w ||^2
  over w in Delta_{Nco}
  lambda.tr_i = t(Lambda_co) %*% w_i
  ```

  By construction, `Y_hat(0) = F %*% lambda.tr` is a convex combination of
  factor-implied control outcomes, so the counterfactual lies pointwise in
  `conv({Y_hat_co_j})` for every time period. Applies only to
  `method = "ife"` with `time.component.from = "nevertreated"` in this
  version.

* New argument `gamma.loading` (default `NULL`): scalar regularization
  strength for the new `"simplex"` projection. `NULL` triggers 5-fold
  cross-validation over a log-grid. When numeric, `gamma.loading` is used
  directly.

* New argument `gamma.loading.grid` (default `NULL`): user-supplied grid for
  `gamma.loading` CV; `NULL` uses `10^seq(-2, 2, length.out = 9)`.

* The unit-FE scalar `alpha.tr` is NOT bounded; under `loading.bound = "simplex"`
  with `force %in% c("unit", "two-way")`, it is computed as the residual-mean
  `mean(U.tr.pre - F.hat.pre %*% t(lambda.tr))` per treated unit.

* Solver: softmax reparameterization with `stats::optim(method = "L-BFGS-B")`
  and an analytic gradient; mirror-descent fallback on ill-conditioned
  inputs. No new R dependencies.

### Diagnostic outputs (under `loading.bound = "simplex"`)

* `loading.bound`: character, records the setting used.
* `gamma.loading`: the value used (CV-selected or user-supplied).
* `loading.proj.resid`: `Ntr`-vector of `|| U.tr.pre - F.hat.pre %*% lambda.tr ||`
  per treated unit. Values substantially above the control-fit RMSE flag
  treated units lying near or outside `conv(Lambda_co)` (simplex constraint
  binds).

### Semantic change to `wgt.implied` (under `loading.bound = "simplex"` only)

* `wgt.implied` becomes the `Ntr x Nco` simplex-weight matrix: each row sums
  to 1 and is non-negative. This is a direct byproduct of the solver and
  replaces the Moore-Penrose pseudo-inverse representation under the bound.
  When `loading.bound = "none"` (default), `wgt.implied` is unchanged.

### Known caveats

* Percentile bootstrap intervals may under-cover when the simplex constraint
  binds (true treated loading on the boundary of `conv(Lambda_co)`); this is
  the Andrews (1999, 2001) non-standard-limit regime for constrained
  estimators. Detect via `loading.proj.resid`. Boundary-corrected inference
  is deferred to a later release.

* v1 does not support the not-yet-treated IFE dispatch, the MC method, or
  the CFE method. `loading.bound = "simplex"` errors cleanly when combined
  with any of these.

## Parallelism cleanup (Phase A bootstrap)

* Phase A's bootstrap error simulation (`R/boot.R::draw.error`) migrated
  from `foreach %dopar%` to `future.apply::future_lapply`. The old
  `%dopar%` inherited whatever backend was registered globally; after any
  prior parallel fect call, `run_dopar_retry`'s `on.exit` left `doFuture`
  registered, so a subsequent call's Phase A inherited a backend that
  shipped heavy closures per iteration --- producing an ~8x slowdown on
  variant (iii) bootstraps in multi-fit sessions (e.g., a forest plot run).

* The `doFuture::registerDoFuture()` re-registration inside
  `run_dopar_retry`'s `on.exit` was removed; it was the source of the
  global state pollution. The function still falls back to `doParallel`
  if the future backend errors; it just no longer leaves a global
  doFuture registration behind.

* New regression test (`tests/testthat/test-phase-a-future-state.R`):
  asserts two consecutive `fect(parallel = TRUE)` calls in the same R
  process have wall-time ratio < 3x.

## GSC: Y.ct.full populated at control positions

* On the GSC path (`method = "ife"` with
  `time.component.from = "nevertreated"`), `Y.ct.full[, co]` is now
  overwritten with the model-implied factor product `F * t(lambda_co)`
  after the shared `Y.ct.full <- Y.ct` assignment (sourced from
  `est.co.best$factor` and `est.co.best$lambda`, gated on dim agreement
  and non-empty rank). Closes a gap that left `NA` at masked control
  positions because the residual recipe `Y.co - residuals` propagates
  `NA`. Enables user-space rolling CV (`r.cv.rolling()`) for
  `method = "gsynth"`.

* No change to ATT, gap, or `est.avg`: those are computed from
  treated-unit positions and do not consume `Y.ct.full[, co]`. Verified
  on the `simgsynth` anchor (set.seed(11), r=2, force="two-way"):
  ATT.avg unchanged at 4.639593 vs unmodified dev; new control-column
  contents match `F * t(lambda.co)` with max abs diff = 0.

# fect 2.2.1

Parametric-bootstrap fixes (`se = TRUE`, `vartype = "parametric"`):

* Parametric bootstrap on unbalanced panels now correctly reflects
  within-unit serial correlation. Prior versions used a diagonal residual
  covariance in the Gaussian draw, under-estimating ATT standard errors on
  serially correlated data by ≈ √((1+ρ)/(1-ρ)) at AR(1) coefficient ρ.
  Balanced panels and `vartype` ∈ {"bootstrap", "jackknife"} are unaffected.
* Fixed `Unsupported bootstrap method: fe` crash when `method = "gsynth"` or
  `"cfe"` and CV selected `r.cv = 0`. Reported against the `gsynth` wrapper,
  which delegates SE to fect.
* Fixed `one.nonpara` dispatcher so `ife+notyettreated`, `cfe+nevertreated`,
  and `cfe+notyettreated` bootstraps route to the correct Loop-2 estimator;
  previously all three produced incorrect SEs. Introduces internal helpers
  `impute_Y0()` and `valid_controls()`.
* Added hard gate erroring on `ife + notyettreated + parametric`; a coverage
  simulation showed ~80% vs 95% nominal. Use `time.component.from = "nevertreated"`
  or a non-parametric `vartype`.

Parallel cross-validation:

* IFE, MC, and CFE cross-validation now run in parallel via `future_lapply`,
  dispatching `(r, fold)` (or `(lambda, fold)`) flat across workers. Auto-engages
  above `Nco * TT > 20000` for IFE/MC, `> 60000` for CFE.
* The `parallel` argument now accepts five forms: `TRUE`, `FALSE`, `"cv"`,
  `"boot"`, `c("cv", "boot")`. Scalar forms are backward-compatible; string
  forms bypass the auto-threshold.
* MC-only tradeoff: `break_check` short-circuits in serial mode only; parallel
  MC evaluates all candidate lambdas. Prefer `parallel = FALSE` for MC if the
  search typically terminates early.
* The future plan is saved and restored on exit (including error paths), so
  caller-set `future::plan()` is unaffected.
* Internal: `(r, fold)` scoring extracted to `R/cv-helpers.R`;
  `fect_nevertreated.R` parallel CV migrated from `foreach %dopar%` (fold-only)
  to `future_lapply` (flat r × k). The migration also resolves a latent
  worker-visibility issue in the prior `foreach` path.
* Fix: vector `parallel = c("cv", "boot")` no longer errors in `fect_boot`,
  `fit_test`, `permutation`, `fect_sens`, or `did_wrapper` — five sites had
  legacy scalar `parallel == TRUE` / `if (parallel)` checks.

# fect 2.2.0

* Added CFE (Complex Fixed Effects) estimator (`method = "cfe"`)
* Added `time.component.from` parameter for latent factor estimation timing
* Added k-fold cross-validation (`cv.sample`) for nevertreated designs
* Improved plot styling: pre/post shading colors, harmonized `type = "esplot"`
* Fixed EM convergence and solver equivalence issues in CFE routines (C++)
* Fixed bootstrap and parallel setup crashes
* Restructured Quarto book with new CFE chapter

# fect 2.1.1

* Added `codetools` to Imports in DESCRIPTION (required by `trim_closure_env()` in `boot.R`)
* Added `importFrom("utils", "tail")` to NAMESPACE (used in `fect_mspe.R`)
* Bumped version to 2.1.1 and updated Date field for CRAN submission

# fect 2.0.4

* Add new plot `type = "hte"`

# fect 2.0.0

* New syntax
* Merged in **gsynth**

# fect 1.0.0

* First CRAN version
* Fixed bugs

# fect 0.6.5

* Replace fastplm with fixest for fixed effects estimation
* Added plots for heterogeneous treatment effects
* Fixed bugs

# fect 0.4.1

* Added a `NEWS.md` file to track changes to the package.
