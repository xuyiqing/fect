# Binary IFE and rolling cross-validation

This development branch ports the working binary initializer and repairs CV
around the current fect implementation. The reference is Licheng Liu's supplied
fect 1.1.10 source, particularly `R/cv_binary.R` and `BiInitialFit` in
`R/support.R`. The branch starts from fect 2.4.5, dev commit
`33fab606cc42d86ae8e1cdbbf9cacfef67e88606`.

## Supported calls

```r
fit <- fect(Y ~ D + X, data = panel, index = c("id", "time"),
            binary = TRUE, QR = TRUE, force = "two-way",
            CV = TRUE, r = c(0, 3), cv.method = "rolling",
            k = 10, cv.nobs = 3, cv.buffer = 1, min.T0 = 5,
            cv.prop = 0.1, cv.rule = "1se", seed = 123,
            se = TRUE, vartype = "jackknife")
```

The outcome must be 0/1. Treatment may follow staggered adoption (common
adoption is included) or switch off again: reversals are supported, and
switch-off effects are reported in `att.off` and `est.att.off`. Observation
weights, `cm = TRUE`, and a never-treated-only fitting route are
unsupported. Inference supports
nonparametric unit/cluster bootstrap and jackknife. Rank is selected once on
the original sample and held fixed in resampling. The bootstrap implementation
can fail when resampled units do not provide enough support; inspect its usual
success diagnostics. Continuous-outcome equivalence diagnostics are not
automatically attached to binary fits.

## What rolling CV holds out

Each fold samples eligible units. Within each sampled unit, it retains earlier
untreated training observations, removes a past buffer, scores an untreated
block, and excludes the remaining future observations from training. Other
units supply contemporaneous untreated observations. The buffer and block
lengths count observed untreated observations before the unit's first treated
period, including on unbalanced panels. With the defaults, a unit needs at
least 5 + 1 + 3 = 9 such observations to be eligible. Untreated observations
after a reversal are never held out and do not count toward eligibility; when
their unit is sampled, they are dropped from training with the rest of its
future. CV eligibility does not redefine the final ATT population.

Every rank uses the same folds. Excluded outcomes are removed before both
initialization and fitting. The primary score is the mean squared error of
predicted probabilities (Brier score), averaged equally across folds.
Classification error at probability 0.5 is a separate diagnostic. Block CV
also scores only its designated cells; donut neighbors are excluded from
training but do not enter the loss denominator.

A rank must complete every fold and the full-sample fit. A failed or unsupported
rank receives infinite selection loss, with reasons retained in `cv.failures`.
If every rank fails, fitting stops. `max.iteration` controls the numerical cap.
Reaching that cap warns on a fixed-rank fit and disqualifies a CV candidate.
The one-SE rule is a tuning heuristic: overlapping folds do not provide an
independent-sample standard error for inference.

Inspect `CV.out`, `cv.loss.per.fold`, `cv.classification.per.fold`, `cv.folds`,
`cv.counts`, `cv.pooled.mspe`, `cv.failures`, `cv.settings`, and
`cv.eligible.units` (rolling CV only). Fold indices address the column-major
outcome matrix `Y.dat`; `cv.id` marks all excluded cells and `est.id` marks
scored cells. `r.cv.rolling(..., binary = TRUE, method = "ife")` uses the same
binary dispatcher and returns fold losses and diagnostics.

## Numerical compatibility

The current C++ numerical algorithms are retained. In particular, the current
QR routine uses the scale acceleration step; the supplied 1.1.10 QR routine
defaults to ordinary EM (`normal = 1`). At the same likelihood-increment
tolerance these can stop at different iterates. Compare like algorithms using
`normal = 0` in the legacy internal routine. The modern IC formula is also
retained, so old and new IC values need not match. IC is diagnostic here and
does not select binary CV rank.

The initializer uses deterministic fixed-effects least squares and PCA.
Changing only the RNG seed does not create different starting values.
Neither successful fits nor CV scores establish consistency or interval
coverage. Paper replication with treatment reversals requires the separately
preserved legacy package.
