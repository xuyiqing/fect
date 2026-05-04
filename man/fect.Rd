\encoding{UTF-8}
\name{fect}
\alias{fect}
\title{Fixed Effects Counterfactual Estimators}
\description{Implements counterfactual estimators in TSCS data analysis and statistical tools to test their identification
assumptions.}
\usage{fect(formula = NULL, data, Y, D, X = NULL,
            W = NULL, W.est = NULL, W.agg = NULL,
            group = NULL,
            na.rm = FALSE,
            index, force = "two-way",
            time.component.from = "notyettreated", em = TRUE,
            r = 0, lambda = NULL, nlambda = 10,
            CV = NULL, k = 20, cv.prop = 0.1, cv.method = "rolling",
            cv.nobs = 3, cv.donut = 1, cv.buffer = 1, criterion = "mspe",
            binary = FALSE, QR = FALSE,
            method = "fe",  se = FALSE, vartype = "bootstrap",
            para.error = "auto", cl = NULL,
            ci.method = "normal", quantile.CI = NULL,
            nboots = 200, alpha = 0.05,
            parallel = TRUE, cores = NULL, tol = 1e-5,
            max.iteration = 5000, seed = NULL,
            min.T0 = NULL, max.missing = NULL,
            proportion = 0.3, pre.periods = NULL,
            f.threshold = 0.5, tost.threshold = NULL,
            knots = NULL, degree = 2,
            sfe = NULL, cfe = NULL,
            Z = NULL, gamma = NULL, Q = NULL, kappa = NULL, 
            Q.type = NULL,
            Q.bspline.degree = NULL,
            Z.param = NULL, Q.param = NULL,
            balance.period = NULL, fill.missing = FALSE,
            placeboTest = FALSE, placebo.period = NULL,
            carryoverTest = FALSE, carryover.period = NULL, carryover.rm = NULL,
            loo = FALSE, permute = FALSE, m = 2,
            normalize = FALSE, keep.sims = FALSE,
            cm = FALSE,
            loading.bound = "none", gamma.loading = NULL,
            gamma.loading.grid = NULL,
            cv.rule = "1se")}
\arguments{
\item{formula}{an object of class "formula": a symbolic description of the model to be fitted, e.g, Y~D+X1+X2}
\item{data}{a data frame, can be a balanced or unbalanced panel data.}
\item{Y}{the outcome indicator.}
\item{D}{the treatment indicator. The treatment should be binary (0 and 1).}
\item{X}{time-varying covariates. Covariates that have perfect collinearity with specified fixed effects are dropped automatically.}
\item{W}{a string giving the column name of a weight variable. Convenience default that populates both \code{W.est} and \code{W.agg} when those are left \code{NULL}. Suitable for survey or sample weights, where the same column applies to both the outcome-model fit and the across-treated-obs aggregation.}

\item{W.est}{a string giving the column name of a weight variable that enters the outcome-model fit (the weighted least squares applied inside the IFE / MC / CFE solver). When \code{NULL}, falls back to \code{W}. Use this (with \code{W.agg = NULL}) when the weight reflects fit-side considerations and the estimand is the unweighted average ATT across treated cells.}

\item{W.agg}{a string giving the column name of a weight variable that enters the across-treated-obs aggregation (\code{att.on}, \code{est.avg}, \code{est.att}). When \code{NULL}, falls back to \code{W}. Use this (with \code{W.est = NULL}) when the user's estimand differs from "ATT for the treated units in the analysis sample" --- common cases include calibration weights to a target population or post-stratification weights that should adjust the summary but not the model fit.

A clean in-package solution for inverse-probability weights for confounding adjustment is under development for fect 3.0 (a cross-fit doubly-robust path); \code{W.agg} is not a substitute for that work and does not deliver the doubly-robust properties an IPW user expects.

In v2.3.1, \code{W.est} and \code{W.agg} (when both supplied) must point to the same column; truly distinct columns for fit vs. aggregation (e.g. combined survey x IPW designs) are scheduled for v2.4.0.}

\item{group}{the group indicator. If specified, the group-wise ATT will be estimated.}
\item{na.rm}{a logical flag indicating whether to list-wise delete missing observations. Default to FALSE. If \code{na.rm = FALSE}, it allows the situation when Y is missing but D is not missing for some observations. If \code{na.rm = TRUE}, it will list-wise delete observations whose Y, D, or X is missing.}
\item{index}{a character vector specifying the unit (first element) and time (second element) indicators. For most methods, must be of length 2. For \code{method = "cfe"}, additional elements (third, fourth, etc.) specify extra fixed-effect grouping variables. Every observation should be uniquely defined by the pair of the unit and time indicator.}
\item{force}{a string indicating whether unit or time or both fixed effects will be imposed. Must be one of the following, "none", "unit", "time", or "two-way". The default is "two-way".}
\item{time.component.from}{Controls which units provide the time-varying model components (time fixed effects, factor structure, temporal dynamics). Options are \code{"notyettreated"} (default) --- all units contribute during their pre-treatment periods, or \code{"nevertreated"} --- only never-treated units estimate the time components, which are then projected onto treated units.}
\item{em}{a logical flag indicating whether to use the EM algorithm for missing data in the estimation sample. Default is \code{TRUE}. Setting \code{em = FALSE} requires a complete estimation sample and is only compatible with \code{time.component.from = "nevertreated"}.}
\item{r}{an integer specifying the number of factors. If \code{CV = TRUE}, the cross validation procedure will select the optimal number of factors from \code{r} to 5.}
\item{lambda}{a single or sequence of positive numbers specifying the hyper-parameter sequence for matrix completion method. If \code{lambda} is a sequence and \code{CV = 1}, cross-validation will be performed.}
\item{nlambda}{an integer specifying the length of hyper-parameter sequence for matrix completion method. Default is \code{nlambda = 10}.}
\item{CV}{a logical flag indicating whether cross-validation will be performed to select the optimal number of factors or hyper-parameter in matrix completion algorithm. If \code{r} is not specified, the procedure will search through \code{r = 0} to \code{5}.}
\item{k}{an integer specifying number of cross-validation rounds. Default is \code{k = 20}.}
\item{cv.prop}{a numerical value specifying the proportion of testing set compared to sample size during the cross-validation procedure.}
\item{cv.method}{a string specifying the cross-validation masking strategy. One of \code{"rolling"} (default; standard time-series rolling-window CV), \code{"block"} (random scattered anchors with contiguous-block masking), or \code{"loo"} (leave-one-out, available for \code{fect_nevertreated}). The legacy aliases \code{"all_units"} (= \code{"block"}) and \code{"treated_units"} (block masking restricted to treated pre-treatment cells) are still accepted but emit a deprecation message; both will be replaced by the unified \code{(cv.method, cv.units)} API in v2.4.0.}
\item{cv.nobs}{an integer specifying the length of continuous observations within a unit in the testing set. Default is \code{cv.nobs = 3}.}
\item{cv.donut}{an integer specifying the length of removed observations at the head and tail of the continuous observations specified by \code{cv.nobs}. Used by block CV (\code{cv.method = "block"}, or the legacy aliases \code{"all_units"}/\code{"treated_units"}). Default is 1 (matches \code{cv.buffer} for rolling CV).}
\item{cv.buffer}{an integer specifying the length of past-side buffer cells masked from training (but not scored) immediately before each rolling-window holdout. Used only by \code{cv.method = "rolling"}; the future side is dropped from training by construction. Analogous to \code{cv.donut} for block CV but applied only on the past side. Default is 1.}
\item{criterion}{criterion used for model selection. Default is "mspe".
  \code{"mspe"} for the mean squared prediction error,
  \code{"gmspe"} for the geometric-mean squared prediction errors,
  \code{"moment"} for period-weighted residuals in test sets,
  \code{"pc"} for an information criterion method.}
\item{binary}{This version doesn't support this option.}
\item{QR}{This version doesn't support this option.}
\item{method}{a string specifying which imputation algorithm will be used.
  \code{"fe"}, \code{"ife"}, \code{"mc"},
  \code{"gsynth"}, or \code{"cfe"}. Default is \code{"fe"}.}
\item{se}{a logical flag indicating whether uncertainty estimates will be produced.}
\item{vartype}{a string specifying the type of variance estimator, e.g.
  \code{"bootstrap"}. Three values are supported:
  \code{"bootstrap"} (nonparametric cluster-bootstrap; the safe default),
  \code{"jackknife"} (leave-one-unit-out), and
  \code{"parametric"} (two-stage pseudo-treated parametric bootstrap).
  The \code{"parametric"} option is restricted to the
  gsynth-style regime: it requires
  \code{time.component.from = "nevertreated"}, no treatment reversal,
  and \code{method} not in \code{c("mc", "both")}.
  These three conditions correspond to Gates A, B, and C in the
  three-gate defense system (see \code{ARCHITECTURE.md}).
  For all other settings, \code{vartype = "bootstrap"} is recommended.}
\item{para.error}{a string specifying the residual-error model used by
  the parametric bootstrap path; sub-option of \code{vartype = "parametric"}
  (silently ignored otherwise). One of \code{"auto"}, \code{"ar"},
  \code{"empirical"}, or \code{"wild"}. Default \code{"auto"} resolves at
  fit time to \code{"empirical"} on a fully-observed panel and \code{"ar"}
  on a panel with missing cells; the resolved label is stored on
  \code{fit$para.error}. \code{"ar"} estimates an AR(1) error process from
  control residuals (works on any panel shape); \code{"empirical"}
  resamples residuals i.i.d. from the main-fit pool (requires fully-
  observed panel); \code{"wild"} applies unit-level Rademacher sign-flips
  over the empirical residual pool (requires fully-observed panel).}
\item{cl}{a string specifying the cluster for cluster bootstrapping.}
\item{ci.method}{a string controlling how the bootstrap distribution
becomes the CIs reported in fect's \code{est.*} slots. Two values:
\code{"normal"} (default; Wald CI: \eqn{\hat\theta \pm z \cdot SE}) or
\code{"basic"} (reflected pivot:
\eqn{[2\hat\theta - q_{1-\alpha/2},\, 2\hat\theta - q_{\alpha/2}]}, the
literature-standard "percentile" CI per Davison & Hinkley (1997, §5.2.1)
and what \code{boot::boot.ci(type = "basic")} returns). For
\code{"percentile"}, \code{"bc"}, or \code{"bca"} CIs (typically wanted
on alternative estimands like \code{aptt} and \code{log.att}), call
\code{\link{estimand}} after fitting.}

\item{quantile.CI}{deprecated as of v2.4.2. Use \code{ci.method} instead.
\code{quantile.CI = FALSE} maps to \code{ci.method = "normal"};
\code{quantile.CI = TRUE} maps to \code{ci.method = "basic"}. Both
mappings still work but emit a one-time deprecation warning.}
\item{nboots}{an integer specifying the number of bootstrap runs. Ignored if \code{se=FALSE}. Default \code{200}, sufficient for the standard error and the normal-CI structure that \code{fect} ships in \code{est.att}/\code{att.avg}. For tail-quantile-based CI methods accessed via \code{estimand()} (\code{ci.method} \code{"basic"}, \code{"percentile"}, \code{"bc"}, \code{"bca"}), bump to \code{nboots = 1000} or higher (Efron 1987 Section 3; DiCiccio & Efron 1996 Section 4); \code{estimand()} emits a warning when called on under-replicated fits.}
\item{alpha}{the significance level for hypothesis tests and confidence intervals. Default \code{0.05}.}
\item{parallel}{controls which operations run in parallel. Accepted values:
\describe{
  \item{\code{TRUE}}{Enable parallel computing for both CV and bootstrap (default in \code{fect()}).}
  \item{\code{FALSE}}{Disable all parallel computing.}
  \item{\code{"cv"}}{Enable parallel CV only; bootstrap runs serially.}
  \item{\code{"boot"}}{Enable parallel bootstrap only; CV runs serially.}
  \item{\code{c("cv","boot")}}{Explicit form of \code{TRUE}: parallel for both.}
}
When \code{parallel = TRUE}, auto-enable thresholds apply: CV parallelism
engages only when \code{Nco * TT} exceeds the per-method threshold
(\code{ife = 20000}, \code{mc = 20000}, \code{cfe = 60000}). Explicit \code{"cv"} overrides
the threshold. Nested parallelism (calling \code{fect()} from within a
\code{future_lapply} or \code{foreach \%dopar\%} block) should use
\code{parallel = FALSE} to avoid deadlock.
When using \code{parallel} with \code{method = "mc"}, parallel CV computes
all candidate lambda values without early stopping (the serial path uses a
\code{break_check} short-circuit to skip lambdas with diminishing MSPE
returns). This guarantees numerical identity between serial and parallel
results but may compute a few extra lambda values compared to the serial
path. Use \code{parallel = FALSE} to preserve the short-circuit behavior.}
\item{cores}{an integer indicating the number of cores for parallel computing.}
\item{tol}{a positive number indicating the relative tolerance for the
EM update check (\code{||fit_new - fit_old||_F / ||fit_old||_F < tol}).
Default tightened from \code{1e-3} to \code{1e-5} in v2.4.3 because the
older default produced under-converged IFE/CFE point estimates that
shifted up to 40\% relative to the converged value (inference was
preserved, but the reported numbers were stopping-point-dependent).}
\item{max.iteration}{the maximal number of EM iterations. Default raised
from \code{1000} to \code{5000} in v2.4.3 to accommodate the tighter tol;
canonical IFE/CFE fits converge in 700-2000 iters at \code{tol = 1e-5}.
A \code{warning()} is emitted when the EM hits this cap without
satisfying the tol gate (under-convergence diagnostic).}
\item{seed}{an integer seed for random number generation.}
\item{min.T0}{an integer specifying the minimum number of pre-treatment periods for each treated unit.}
\item{max.missing}{an integer specifying the maximum number of missing observations allowed per unit.}
\item{proportion}{a numeric value specifying which pre-treatment periods are used for goodness-of-fit tests.}
\item{pre.periods}{a vector specifying the range of pre-treatment periods used for the goodness-of-fit test.}
\item{f.threshold}{a numeric threshold for an F-test in equivalence testing. Default \code{0.5}.}
\item{tost.threshold}{a numeric threshold for two-one-sided t-tests.}
\item{knots}{a numeric vector specifying knots (currently unused; reserved for future use).}
\item{degree}{an integer specifying the degree (currently unused; reserved for future use).}
\item{sfe}{vector specifying other fixed effects for \code{method="cfe"}.}
\item{cfe}{a vector of lists specifying interactive fixed effects for \code{method="cfe"}.}
\item{Z}{a vector specifying the time-invariant covariates for the Z matrix.}
\item{gamma}{a vector specifying the time-varying covariates for the gamma matrix.}
\item{Q}{a vector specifying the time-varying covariates for the Q matrix.}
\item{kappa}{a vector specifying the time-invariant covariates for the kappa matrix.}
\item{Q.type}{a vector specifying the type of Q matrix.}
\item{Q.bspline.degree}{an integer specifying the degree used when \code{Q.type} includes \code{"bspline"} in \code{method="cfe"}. If NULL, a default degree is chosen based on the number of distinct time values.}
\item{Z.param}{a list specifying the parameters for the Z matrix.}
\item{Q.param}{a list specifying the parameters for the Q matrix.}
\item{balance.period}{a length-2 vector specifying a time range for a balanced sample.}
\item{fill.missing}{a logical flag indicating whether to allow missing observations in a balanced sample.}
\item{placeboTest}{a logical flag indicating whether to perform a placebo test.}
\item{placebo.period}{an integer or 2-element numeric vector specifying pseudo-treatment periods.}
\item{carryoverTest}{a logical flag for carryover tests.}
\item{carryover.period}{an integer or 2-element numeric vector specifying pseudo-carryover periods.}
\item{carryover.rm}{an integer specifying the range of post-treatment periods to treat as carryover.}
\item{loo}{a logical flag for leave-one-period-out goodness-of-fit tests.}
\item{permute}{a logical flag indicating whether to run a permutation test.}
\item{m}{an integer specifying the block length for permutation tests. Default \code{2}.}
\item{normalize}{a logical flag indicating whether to scale outcome and covariates.}
\item{keep.sims}{a logical flag indicating whether to save unit-time level bootstrap effects. Default \code{keep.sims = FALSE}. If \code{se = FALSE}, this argument is ignored.}
\item{cm}{a logical flag indicating whether to enable causal moderation analysis. When \code{TRUE}, the estimator decomposes the treatment effect into effect modification and causal moderation components. Currently available for \code{method = "fe"} and \code{method = "ife"}. Default is \code{FALSE}.}
\item{loading.bound}{a string controlling whether treated-unit factor loadings are bounded inside the convex hull of control loadings. \code{"none"} (default) reproduces standard GSC behavior. \code{"simplex"} constrains each treated unit's loading to be a non-negative convex combination of control loadings via an entropy-regularized simplex projection, ensuring the imputed counterfactual lies pointwise in the convex hull of factor-implied control outcomes. Currently applies only to \code{method = "ife"} or \code{method = "gsynth"} (equivalent forms), and requires \code{time.component.from = "nevertreated"}.}
\item{gamma.loading}{scalar regularization strength for the \code{"simplex"} projection. \code{NULL} (default) triggers 5-fold cross-validation over \code{gamma.loading.grid}. A numeric value is used directly. Ignored when \code{loading.bound = "none"}.}
\item{gamma.loading.grid}{a numeric vector of candidate \code{gamma.loading} values for cross-validation. \code{NULL} (default) uses \code{10^seq(-2, 2, length.out = 9)}. Ignored when \code{loading.bound = "none"} or when \code{gamma.loading} is supplied.}
\item{cv.rule}{a string selecting the cross-validation rule for choosing the number of factors \code{r} (or matrix-completion penalty \code{lambda}). One of:
\describe{
  \item{\code{"1se"} (default)}{The 1-SE rule (Breiman, Friedman, Olshen and Stone 1984; Hastie, Tibshirani and Friedman 2009, Section 7.10): pick the smallest \code{r} whose mean CV criterion is within one fold-SE of the minimum-CV-error \code{r}. Biases toward parsimony in a fold-aware way --- when CV is precise, it allows larger \code{r}; when CV is noisy, it gravitates to simpler models.}
  \item{\code{"min"}}{Pick the \code{r} that minimizes the mean CV criterion (no tolerance).}
  \item{\code{"1pct"}}{Legacy pre-2.3.0 heuristic: pick the smallest \code{r} within 1\% relative tolerance of the best mean CV criterion. Use this for byte-identical reproducibility of pre-2.3.0 fits.}
}
Ignored when \code{CV = FALSE}.}
}
\details{
  \code{fect} implements counterfactual estimators for TSCS data. It first imputes counterfactuals by fitting an outcome model using untreated observations, then estimates the individual treatment effect as the difference between observed and predicted outcomes. Finally, it computes average treatment effects on the treated (ATT) and period-specific ATTs. Placebo and equivalence tests help evaluate identification assumptions.
}
\value{
  \item{Y.dat}{T-by-N matrix of the outcome variable.}
  \item{D.dat}{T-by-N matrix of the treatment variable.}
  \item{I.dat}{T-by-N matrix of observation indicators (observed/missing).}
  \item{Y}{name of the outcome variable.}
  \item{D}{name of the treatment variable.}
  \item{X}{name of any time-varying covariates.}
  \item{W}{name of the weight variable.}
  \item{index}{name of the unit and time indicators.}
  \item{force}{specified fixed effects option.}
  \item{T}{number of time periods.}
  \item{N}{number of units.}
  \item{p}{number of time-varying observables.}
  \item{r.cv}{number of factors (selected by cross-validation if needed).}
  \item{lambda.cv}{optimal hyper-parameter for matrix completion, if applicable.}
  \item{beta}{coefficients for any covariates in an interactive fixed effects model.}
  \item{sigma2}{mean squared error.}
  \item{IC}{information criterion.}
  \item{est}{results of the fitted model.}
  \item{MSPE}{mean squared prediction error from cross-validation.}
  \item{CV.out}{results of the cross-validation procedure.}
  \item{niter}{number of iterations.}
  \item{factor}{estimated time-varying factors.}
  \item{lambda}{estimated loadings.}
  \item{lambda.tr}{estimated loadings for treated units.}
  \item{lambda.co}{estimated loadings for control units.}
  \item{mu}{estimated grand mean.}
  \item{xi}{estimated time fixed effects.}
  \item{alpha}{estimated unit fixed effects.}
  \item{alpha.tr}{estimated unit fixed effects for treated units.}
  \item{alpha.co}{estimated unit fixed effects for control units.}
  \item{validX}{logical indicating if valid covariates exist.}
  \item{validF}{logical indicating if factors exist.}
  \item{id}{vector of unit IDs.}
  \item{rawtime}{vector of time periods.}
  \item{obs.missing}{matrix indicating missingness patterns.}
  \item{Y.ct}{T-by-N matrix of predicted outcomes under no treatment.}
  \item{eff}{T-by-N matrix of estimated individual treatment effects.}
  \item{res}{residuals for observed values.}
  \item{eff.pre}{effects for treated units in pre-treatment periods.}
  \item{eff.pre.equiv}{pre-treatment effects under baseline (two-way FE) model.}
  \item{pre.sd}{by-period residual standard deviations for pre-treatment ATT.}
  \item{att.avg}{overall average treatment effect on the treated.}
  \item{att.avg.W}{weighted ATT.}
  \item{att.avg.unit}{by-unit average treatment effect on the treated.}
  \item{time}{time index for switch-on treatment effect.}
  \item{count}{count of observations for each switch-on effect time.}
  \item{att}{switch-on treatment effect.}
  \item{att.on.W}{weighted switch-on effect.}
  \item{time.off}{time index for switch-off treatment effect.}
  \item{att.off}{switch-off treatment effect.}
  \item{att.off.W}{weighted switch-off effect.}
  \item{count.off}{count for each switch-off period.}
  \item{att.placebo}{ATT for placebo periods.}
  \item{att.carryover}{ATT for carryover periods.}
  \item{eff.calendar}{ATT by calendar time.}
  \item{eff.calendar.fit}{loess-fitted ATT by calendar time.}
  \item{N.calandar}{number of treated observations each calendar period.}
  \item{balance.avg.att}{ATT for balanced sample.}
  \item{balance.att}{switch-on ATT for balanced sample.}
  \item{balance.time}{time index for balanced sample.}
  \item{balance.count}{count for each time in balanced sample.}
  \item{balance.att.placebo}{ATT for placebo period in balanced sample.}
  \item{group.att}{ATT for different groups.}
  \item{group.output}{list of switch-on treatment effects by group.}
  \item{est.att.avg}{inference for \code{att.avg}.}
  \item{est.att.avg.unit}{inference for \code{att.avg.unit}.}
  \item{est.att}{inference for \code{att}.}
  \item{est.att.W}{inference for weighted \code{att}.}
  \item{est.att.off}{inference for switch-off.}
  \item{est.att.off.W}{inference for weighted switch-off.}
  \item{est.placebo}{inference for placebo ATT.}
  \item{est.carryover}{inference for carryover ATT.}
  \item{est.eff.calendar}{inference for \code{eff.calendar}.}
  \item{est.eff.calendar.fit}{inference for \code{eff.calendar.fit}.}
  \item{est.balance.att}{inference for balanced sample switch-on.}
  \item{est.balance.avg}{inference for balanced sample average ATT.}
  \item{est.balance.placebo}{inference for balanced sample placebo.}
  \item{est.avg.W}{inference for \code{att.avg.W}.}
  \item{est.beta}{inference for \code{beta}.}
  \item{est.group.att}{inference for group-specific ATT.}
  \item{est.group.output}{inference for group output.}
  \item{att.avg.boot}{bootstrap draws for \code{att.avg}.}
  \item{att.avg.unit.boot}{bootstrap draws for \code{att.avg.unit}.}
  \item{att.count.boot}{bootstrap draws for \code{count}.}
  \item{att.off.boot}{bootstrap draws for \code{att.off}.}
  \item{att.off.count.boot}{bootstrap draws for \code{count.off}.}
  \item{att.placebo.boot}{bootstrap draws for \code{att.placebo}.}
  \item{att.carryover.boot}{bootstrap draws for \code{att.carryover}.}
  \item{balance.att.boot}{bootstrap draws for \code{balance.att}.}
  \item{att.bound}{equivalence confidence interval for pre-trend.}
  \item{att.off.bound}{equivalence confidence interval for switch-off.}
  \item{beta.boot}{bootstrap draws for \code{beta}.}
  \item{test.out}{F-test and equivalence test results for pre-treatment fit.}
  \item{loo.test.out}{leave-one-period-out test results.}
  \item{permute}{permutation test results.}
}
\author{
  Licheng Liu; Ye Wang; Yiqing Xu; Ziyi Liu
}
\references{
  Athey, S., Bayati, M., Doudchenko, N., Imbens, G., and Khosravi, K. (2021).
  Matrix completion methods for causal panel data models.
  \emph{Journal of the American Statistical Association, 116}(536), 1716-1730.

  Bai, J. (2009).
  Panel data models with interactive fixed effects.
  \emph{Econometrica, 77}(4), 1229-1279.

  Liu, L., Wang, Y., and Xu, Y. (2022).
  A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional Data.
  \emph{American Journal of Political Science, 68}(1), 160-176.

  Xu, Y. (2017).
  Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models.
  \emph{Political Analysis, 25}(1), 57-76.
}
\seealso{
  \code{\link{plot.fect}}, \code{\link{print.fect}}
}
\examples{
library(fect)
data(simdata)
out <- fect(Y ~ D + X1 + X2, data = simdata,
            index = c("id","time"), force = "two-way",
            CV = TRUE, r = c(0, 5), se = 0, parallel = FALSE)
}
