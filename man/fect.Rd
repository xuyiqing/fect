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
            r = NULL, lambda = NULL, nlambda = 10,
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
            group.fe = NULL,
            cfe = NULL,
            Z = NULL, gamma = NULL, Q = NULL, kappa = NULL,
            Q.type = NULL,
            Q.bspline.degree = NULL,
            Z.param = NULL, Q.param = NULL,
            balance.period = NULL, fill.missing = FALSE,
            placeboTest = FALSE, placebo.period = NULL,
            carryoverTest = FALSE, carryover.period = NULL, carryover.rm = NULL,
            loo = FALSE, dloo = FALSE, dloo.adjust = FALSE,
            permute = FALSE, m = 2,
            normalize = FALSE, keep.sims = FALSE,
            cm = FALSE,
            loading.bound = "none", gamma.loading = NULL,
            gamma.loading.grid = NULL,
            cv.rule = "min")}
\arguments{
\item{formula}{an object of class "formula" naming the outcome, the treatment and the covariates, e.g. \code{Y ~ D + X1 + X2}: the outcome on the left, the treatment first on the right, then the covariates. Every term must be a bare column name of \code{data}. Other terms stop with a message that names them: transformations such as \code{log(Y)} or \code{I(X1^2)}, \code{factor(g)}, and interactions such as \code{X1:X2} or \code{X1 * X2}. Create such variables as columns first (e.g. \code{data$logY <- log(data$Y)}; dummy columns with \code{model.matrix()}). Intercept terms (\code{+ 0}, \code{- 1}, \code{+ 1}) are accepted and ignored, because \code{fect} always includes its own fixed effects. With a formula, give the covariates in the formula; also passing \code{X} stops.}
\item{data}{a data frame, can be a balanced or unbalanced panel data.}
\item{Y}{the outcome indicator.}
\item{D}{the treatment indicator. The treatment should be binary (0 and 1).}
\item{X}{a character vector of covariate (time-varying) column names, used when no \code{formula} is given. Each name must be a column of \code{data}, given once, and not the outcome or the treatment. Covariates must be numeric: a logical column is used as 0/1, and a factor, character or other column stops (create numeric dummy columns first, e.g. with \code{model.matrix()}). A covariate that cannot be estimated on the cells used to fit the model is dropped with a warning, and its coefficient is reported as \code{NA}: one with no variation on those cells, one absorbed by the fixed effects (e.g. constant over time within each unit, with unit fixed effects), or an exact linear combination of the other covariates after the fixed effects are removed (the first covariate of a collinear set is kept, as in \code{lm()}). The estimates then equal those of the fit without the dropped covariates. The cells used to fit the model are the untreated, observed cells; with \code{time.component.from = "nevertreated"} (and \code{method = "gsynth"}), only those of the never-treated units.}
\item{W}{a string giving the column name of a weight variable. Convenience default that populates both \code{W.est} and \code{W.agg} when those are left \code{NULL}. Suitable for survey or sample weights, where the same column applies to both the outcome-model fit and the across-treated-obs aggregation. Weights work with all three \code{vartype} options.}

\item{W.est}{a string giving the column name of a weight variable that enters the outcome-model fit (the weighted least squares applied inside the IFE / MC / CFE solver). When \code{NULL}, falls back to \code{W}. Use this (with \code{W.agg = NULL}) when the weight reflects fit-side considerations and the estimand is the unweighted average ATT across treated cells.}

\item{W.agg}{a string giving the column name of a weight variable that enters the across-treated-obs aggregation (\code{att.on}, \code{est.avg}, \code{est.att}, and the summaries of \code{\link{estimand}}, \code{\link{effect}} and \code{\link{att.cumu}}). When \code{NULL}, falls back to \code{W}. Use this (with \code{W.est = NULL}) when the user's estimand differs from "ATT for the treated units in the analysis sample" --- common cases include calibration weights to a target population or post-stratification weights that should adjust the summary but not the model fit. Given alone, \code{W.agg} enters only the aggregation: never the model fit, and never the cross-validation scores when \code{CV = TRUE}.

A clean in-package solution for inverse-probability weights for confounding adjustment is under development for fect 3.0 (a cross-fit doubly-robust path); \code{W.agg} is not a substitute for that work and does not deliver the doubly-robust properties an IPW user expects.

When both are given, \code{W.est} and \code{W.agg} must name the same column; distinct columns for the fit and the aggregation (e.g. combined survey x IPW designs) are not supported yet.}

\item{group}{the group indicator. If specified, the group-wise ATT will be estimated.}
\item{na.rm}{a logical flag indicating whether to list-wise delete missing observations. Default to FALSE. If \code{na.rm = FALSE}, it allows the situation when Y is missing but D is not missing for some observations. If \code{na.rm = TRUE}, it will list-wise delete observations whose Y, D, or X is missing.}
\item{index}{a character vector specifying the unit (first element) and time (second element) indicators. For most methods, must be of length 2. For \code{method = "cfe"}, additional elements (third, fourth, etc.) specify extra fixed-effect grouping variables. Each label belongs to a unit-period, so a unit may change groups; on an unbalanced panel, a unit whose observed rows all carry one label keeps it at its absent periods. Every observation should be uniquely defined by the pair of the unit and time indicator. The time indicator can be numeric, a \code{Date}, a factor, or character. A factor is used in its level order: if its levels are increasing numbers, \code{fect} uses those numbers; otherwise it uses the level positions and keeps the levels as labels (in \code{rawtime}, the row names of \code{eff}, and \code{data.long}). If the levels are numbers but not in increasing order (e.g. \code{"1", "10", "2"}, as \code{factor(as.character(time))} gives for 1, 2, ..., 10), \code{fect} warns and still uses the level order; pass a numeric index or reorder the levels if that is not the time order. A character time indicator must hold numbers (e.g. \code{"2001"}), which are used as numbers; other values stop with a message (convert them to numbers, a \code{Date}, or a factor whose levels are in time order).}
\item{force}{a string indicating whether unit or time or both fixed effects will be imposed. Must be one of the following, "none", "unit", "time", or "two-way". The default is "two-way".}
\item{time.component.from}{Controls which units provide the time-varying model components (time fixed effects, factor structure, temporal dynamics). Options are \code{"notyettreated"} (default) --- all units contribute during their pre-treatment periods, or \code{"nevertreated"} --- only never-treated units estimate the time components, which are then projected onto treated units. With \code{"nevertreated"}, \code{method = "ife"} fits the same model as \code{method = "gsynth"}, with or without \code{se}, and the fit reports \code{method = "gsynth"}.}
\item{em}{a logical flag indicating whether to use the EM algorithm for missing data in the estimation sample. Default is \code{TRUE}. Setting \code{em = FALSE} requires a complete estimation sample and is only compatible with \code{time.component.from = "nevertreated"}.}
\item{r}{an integer specifying the number of factors, or a range such as \code{c(0, 5)}. Default \code{NULL}: for \code{method = "ife"} or \code{"gsynth"} the number of factors is then selected by cross-validation over 0 to 5, as \code{method = "mc"} does for \code{lambda}; for \code{method = "fe"} or \code{"cfe"} it means no factors. With \code{CV = FALSE} and no \code{r}, the fit falls back to \code{r = 0} (the FEct model) with a message. A single number with \code{CV = FALSE} (the default when \code{r} is a single number) fixes the number of factors; \code{r = 0} then estimates no factors. With \code{CV = TRUE}, a range \code{c(a, b)} is searched from \code{a} to \code{b}; a single number is the only candidate, so no search takes place. When the factors are estimated from the never-treated units only (\code{method = "gsynth"}, or \code{time.component.from = "nevertreated"}), at most \code{Nco - 1} factors can be estimated, where \code{Nco} is the number of never-treated units; cross-validation then searches no further, with a message.}
\item{lambda}{a single or sequence of positive numbers specifying the hyper-parameter sequence for matrix completion method. If \code{lambda} is a sequence and \code{CV = 1}, cross-validation will be performed.}
\item{nlambda}{an integer specifying the length of hyper-parameter sequence for matrix completion method. Default is \code{nlambda = 10}.}
\item{CV}{a logical flag indicating whether cross-validation will be performed to select the optimal number of factors or hyper-parameter in matrix completion algorithm. If \code{r} is not specified, the procedure searches through \code{r = 0} to \code{5}; for \code{method = "ife"} and \code{"gsynth"} this is the default behaviour.}
\item{k}{an integer specifying number of cross-validation rounds. Default is \code{k = 20}.}
\item{cv.prop}{a number between 0 and 1. With \code{cv.method = "rolling"} (the default), the share of eligible units sampled in each fold; with \code{cv.method = "block"}, the share of observations masked in each fold. Default \code{0.1}.}
\item{cv.method}{a string specifying the cross-validation masking strategy. One of \code{"rolling"} (default; standard time-series rolling-window CV), \code{"block"} (random scattered anchors with contiguous-block masking), or \code{"loo"} (leave-one-out, available only with \code{time.component.from = "nevertreated"}, including \code{method = "gsynth"}). The legacy names \code{"all_units"} (the same as \code{"block"}) and \code{"treated_units"} (block masking with anchors in treated pre-treatment cells) still work and print a note.}
\item{cv.nobs}{an integer specifying the length of continuous observations within a unit in the testing set. Default is \code{cv.nobs = 3}.}
\item{cv.donut}{an integer specifying the length of removed observations at the head and tail of the continuous observations specified by \code{cv.nobs}. Used by block CV (\code{cv.method = "block"}, or the legacy aliases \code{"all_units"}/\code{"treated_units"}). Default is 1 (matches \code{cv.buffer} for rolling CV). When \code{cv.nobs >= 3}, \code{2 * cv.donut} must be smaller than \code{cv.nobs}, so that each held-out run keeps at least one scored period.}
\item{cv.buffer}{an integer specifying the length of past-side buffer cells masked from training (but not scored) immediately before each rolling-window holdout. Used only by \code{cv.method = "rolling"}; the future side is dropped from training by construction. Analogous to \code{cv.donut} for block CV but applied only on the past side. Default is 1.}
\item{criterion}{criterion used for model selection. Default is "mspe".
  \code{"mspe"} for the mean squared prediction error,
  \code{"gmspe"} for the geometric-mean squared prediction errors,
  \code{"moment"} for period-weighted residuals in test sets,
  \code{"pc"} for an information criterion method: the \code{r} with the
  lowest PC value is selected, also for \code{method = "gsynth"} and for
  \code{method = "cfe"} with \code{time.component.from = "nevertreated"}
  (\code{cv.rule} does not apply). Under \code{"pc"}, the IFE
  cross-validation table (\code{CV.out.ife}) has the columns \code{r},
  \code{sigma2}, \code{IC}, \code{PC}, \code{MSPTATT} and \code{MSE}.}
\item{binary}{Not supported in this version: \code{binary = TRUE} stops with a message. Leave it at \code{FALSE} (the default).}
\item{QR}{This version doesn't support this option.}
\item{method}{a string specifying which imputation algorithm will be used.
  \code{"fe"}, \code{"ife"}, \code{"mc"},
  \code{"gsynth"}, \code{"cfe"}, or \code{"both"}. Default is \code{"fe"}.
  \code{"both"} cross-validates \code{"ife"} (over \code{r}) and \code{"mc"}
  (over \code{lambda}) and fits the one with the lower cross-validation
  error. \code{"fe"} is fitted as \code{"ife"} with \code{r = 0}.}
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
  It also needs at least two never-treated units, and it ignores
  \code{cl} (with a warning).
  For all other settings, \code{vartype = "bootstrap"} is recommended.
  The variance type used is stored in the fit as \code{vartype}, which
  \code{\link{effect}} and \code{\link{att.cumu}} read.}
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
\item{cl}{a string specifying the cluster column for cluster bootstrapping:
the name of one column of \code{data}. The cluster bootstrap resamples
whole units by cluster, so with \code{se = TRUE} and
\code{vartype = "bootstrap"} the cluster must be constant within each unit,
and there must be at least two clusters; \code{fect} stops otherwise.
Clusters may differ in size. \code{vartype = "jackknife"} and
\code{vartype = "parametric"} ignore \code{cl}, with a warning.
When \code{group.fe} is set to a single column and \code{cl} is unset,
\code{cl} auto-defaults to \code{group.fe[1]} --- the natural choice when
treatment varies at the group level (Bertrand, Duflo & Mullainathan 2004).
To override --- for example, to cluster at the unit level when group-level
FE is absorbed --- pass \code{cl = index[1]} (e.g., \code{cl = "id"}).
Note that \code{cl = NULL} does NOT disable clustering: the case bootstrap
always resamples units (which is unit-level clustering); the auto-default
only changes the resample unit to a coarser level. \code{cl = FALSE} is
rejected with a guiding error.}
\item{ci.method}{a string controlling how the bootstrap distribution
becomes the CIs reported in fect's \code{est.*} slots. Two values:
\code{"normal"} (default; Wald CI: \eqn{\hat\theta \pm z \cdot SE}) or
\code{"basic"} (reflected pivot:
\eqn{[2\hat\theta - q_{1-\alpha/2},\, 2\hat\theta - q_{\alpha/2}]}, the
literature-standard "percentile" CI per Davison & Hinkley (1997, §5.2.1)
and what \code{boot::boot.ci(type = "basic")} returns). With
\code{"basic"}, the p-values come from the replicates. Each p-value goes
with the basic interval: it is below \eqn{\alpha} exactly when the
\eqn{1 - \alpha} basic interval excludes zero, at every \eqn{\alpha}, for
the effects and for the coefficients (\code{est.beta}). For parametric
fits, whose replicates are simulated with no effect, the replicates are
first recentered at the estimate, the same shift the basic interval uses,
and the p-value is the one dual to that interval. For
\code{"percentile"}, \code{"bc"}, or \code{"bca"} CIs (typically wanted
on alternative estimands like \code{aptt} and \code{log.att}), call
\code{\link{estimand}} after fitting.}

\item{quantile.CI}{deprecated as of v2.4.2. Use \code{ci.method} instead.
\code{quantile.CI = FALSE} maps to \code{ci.method = "normal"};
\code{quantile.CI = TRUE} maps to \code{ci.method = "basic"}. Both
mappings still work but emit a one-time deprecation warning.}
\item{nboots}{an integer specifying the number of bootstrap runs. Ignored if \code{se=FALSE}. Default \code{200}, sufficient for the standard error and the normal-CI structure that \code{fect} ships in \code{est.att}/\code{att.avg}. For tail-quantile-based CI methods accessed via \code{estimand()} (\code{ci.method} \code{"basic"}, \code{"percentile"}, \code{"bc"}, \code{"bca"}), bump to \code{nboots = 1000} or higher (Efron 1987 Section 3; DiCiccio & Efron 1996 Section 4); \code{estimand()} emits a warning when called on under-replicated fits.}
\item{alpha}{the significance level. Default \code{0.05}. The intervals in \code{est.avg}, \code{est.att} and the other \code{est.*} tables are (1 - \code{alpha}) intervals. \code{att.bound}, \code{est.att90}, \code{att.off.bound} and the \code{(90\%)} columns of \code{est.placebo} and \code{est.carryover} are (1 - 2\code{alpha}) intervals, 90\% at the default; the equivalence plot draws them. The p-values and the tests in \code{test.out} do not depend on \code{alpha}, and \code{Y.avg} always holds 95\% and 90\% intervals.}
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
\item{cores}{an integer indicating the number of cores for parallel computing. Default \code{NULL}: the number of available cores minus two, at most 8.}
\item{tol}{a positive number indicating the relative tolerance for the
EM update check (\code{||fit_new - fit_old||_F / ||fit_old||_F < tol}). With \code{CV = TRUE} the cross-validation folds use \code{max(tol, 1e-3)} and the chosen model is then refit at \code{tol}.
Default tightened from \code{1e-3} to \code{1e-5} in v2.4.3 because the
older default produced under-converged IFE/CFE point estimates that
shifted up to 40\% relative to the converged value (inference was
preserved, but the reported numbers were stopping-point-dependent).}
\item{max.iteration}{the maximal number of EM iterations. Default raised
from \code{1000} to \code{5000} in v2.4.3 to accommodate the tighter tol;
canonical IFE/CFE fits converge in 700-2000 iters at \code{tol = 1e-5}.
A \code{warning()} is emitted when the EM hits this cap without
satisfying the tol gate (under-convergence diagnostic).}
\item{seed}{an integer seed for random number generation. When supplied, it fixes the cross-validation folds and the bootstrap and permutation draws, so repeating a call gives the same results; the same \code{seed} also gives the same cross-validation folds, and so the same \code{r}, with and without \code{se}. Default \code{NULL}: the cross-validation folds and the serial bootstrap (\code{parallel = FALSE}) follow the current state of R's random number generator (see \code{\link{set.seed}}), but the parallel bootstrap draws do not, so calling \code{set.seed()} before \code{fect()} does not make them reproducible. For reproducible bootstrap draws, pass \code{seed}.}
\item{min.T0}{an integer: the minimum number of observations under control a unit needs; units with fewer are removed before estimation (see \code{remove.id} in the value). Default \code{NULL}: 1 for \code{method = "fe"} and 5 for the other methods.}
\item{max.missing}{an integer: the maximum number of missing periods a unit may have. A period is missing for a unit when its row is absent or its outcome or a covariate is missing. Units with more missing periods are removed before estimation (see \code{remove.id} in the value). Default \code{NULL}: no limit.}
\item{proportion}{a numeric value specifying which pre-treatment periods are used for goodness-of-fit tests. In cross-validation, the same periods count in the weighted scores (WMSPE, WGMSPE, Moment and GMoment); periods with at most \code{proportion} times the largest number of treated units get zero weight.}
\item{pre.periods}{a vector specifying the range of pre-treatment periods used for the goodness-of-fit test.}
\item{f.threshold}{the equivalence bound of the F test of the pre-treatment ATTs (\code{se = TRUE}). \code{test.out$f.equiv.p} is the probability that an F statistic with noncentrality \code{N_bar * f.threshold} is no larger than the observed one, where \code{N_bar} is the largest number of treated units in the tested event times; a small value supports equivalence. Default \code{0.5}.}
\item{tost.threshold}{a numeric threshold for the two one-sided tests (TOST) of equivalence: the equivalence range is \code{[-tost.threshold, tost.threshold]}. Default \code{NULL}: 0.36 times the residual standard deviation of the fixed-effects model, \code{0.36 * sqrt(sigma2.fect)}.}
\item{knots}{a numeric vector specifying knots (currently unused; reserved for future use).}
\item{degree}{an integer specifying the degree (currently unused; reserved for future use).}
\item{group.fe}{a character vector of column names naming additional simple
additive fixed-effect groupings to absorb (e.g., \code{group.fe = "state"}
when rows are counties and treatment varies at state level). Each entry must
be a column in \code{data}. Each entry must be \emph{nested in}
\code{index[1]} --- i.e., constant within each level of the unit identifier
--- otherwise an error is raised. When \code{group.fe} is set, \code{method
= "fe"} is silently routed to \code{method = "cfe"} (FE is a subset of CFE,
identical result); \code{method = "ife"}, \code{"mc"}, \code{"both"}, and
\code{"gsynth"} hard-error (use \code{method = "cfe"} with \code{r > 0} for
free latent factors with group-level FE). When \code{group.fe} has length 1
and \code{cl} is unset, \code{cl} auto-defaults to \code{group.fe[1]}. To
cluster at the unit level instead, pass \code{cl = index[1]} (e.g.,
\code{cl = "id"}); \code{cl = NULL} does NOT change behavior (the case
bootstrap always resamples units regardless of \code{cl}), and
\code{cl = FALSE} is rejected with a guiding error.}
\item{cfe}{not used: this argument is ignored and is kept only so that older code still runs. Add extra additive fixed effects with \code{group.fe} (or extra columns in \code{index}), and factors with \code{r}.}
\item{Z}{a vector of column names: time-invariant unit covariates whose coefficients vary by period (grouped by \code{gamma}). Each must be the same in every period of a unit; a value that varies within a unit stops with a message.}
\item{gamma}{a column name giving each period's group for the coefficients on \code{Z}; it must be the same for every unit in a period. Default: one group per period.}
\item{Q}{a vector of column names: known time trends whose loadings vary by unit (grouped by \code{kappa}). Each must be the same for every unit in a period; a value that varies within a period stops with a message. See also \code{Q.type}.}
\item{kappa}{a column name giving each unit's group for the loadings on \code{Q}; it must be the same in every period of a unit. Default: one group per unit.}
\item{Q.type}{a vector of known time trends to build as \code{Q} for \code{method = "cfe"}: any of \code{"linear"} (or \code{1}), \code{"quadratic"} (\code{2}), \code{"cubic"} (\code{3}) and \code{"bspline"}. Cannot be combined with \code{Q}.}
\item{Q.bspline.degree}{an integer specifying the degree used when \code{Q.type} includes \code{"bspline"} in \code{method="cfe"}. If NULL, a default degree is chosen based on the number of distinct time values.}
\item{Z.param}{a list specifying the parameters for the Z matrix.}
\item{Q.param}{a list specifying the parameters for the Q matrix.}
\item{balance.period}{a length-2 vector \code{c(a, b)}, with \code{a <= 0} and \code{b >= 1}, of event times (1 is the first treated period, 0 the last one before it). The treated units whose event times cover every period from \code{a} to \code{b} form a balanced sample, and its effects are reported in the \code{balance.*} and \code{est.balance.*} slots. Cannot be combined with \code{group} or \code{carryover.rm}.}
\item{fill.missing}{a logical flag, used only with \code{balance.period}. With \code{FALSE} (the default), a treated unit enters the balanced sample only if its outcome is observed in every period of the window. With \code{TRUE}, periods of the window whose outcome or covariates are missing, or whose rows are absent, do not exclude the unit; the balanced-sample estimates then average its observed cells.}
\item{placeboTest}{a logical flag indicating whether to perform a placebo test: the periods in \code{placebo.period} are held out of the fit and their effects estimated as if they were treated. Needs \code{placebo.period}. Cannot be combined with \code{carryoverTest = TRUE}, \code{CV = TRUE} (give a single \code{r} or \code{lambda}) or \code{permute = TRUE}; \code{loo = TRUE} turns it off.}
\item{placebo.period}{the pre-treatment periods to hold out for the placebo test, in time relative to treatment onset (1 is the first treated period, 0 the last period before it): one period, or a range \code{c(a, b)} with \code{a <= b <= 0}, e.g. \code{c(-2, 0)}. Required when \code{placeboTest = TRUE}.}
\item{carryoverTest}{a logical flag for carryover tests.}
\item{carryover.period}{an integer or 2-element numeric vector specifying pseudo-carryover periods.}
\item{carryover.rm}{a positive integer \code{k}: for each treated unit, the observations whose time value is 1 to \code{k} above a treated period, but not themselves treated, are removed from the fit, so possible carryover effects after treatment ends do not enter the estimation. With a time index that steps by 1, these are the first \code{k} periods after treatment ends. Cannot be combined with \code{balance.period}.}
\item{loo}{a logical flag for leave-one-period-out goodness-of-fit tests. Each pre-treatment period is held out in turn and the model is refitted with the fit's settings, including the selected \code{r} (or \code{lambda}), \code{loading.bound} with its \code{gamma.loading}, \code{time.component.from} and \code{para.error}.}
\item{dloo}{a logical flag for the double (cohort-wise) leave-one-out pre-trend
placebos. When \code{TRUE}, \code{fect} estimates these directly from the in-sample fit -- without re-fitting
the imputation model. With covariates, the covariate coefficients are re-estimated for each cohort and tested
period on that cohort's comparison sample (its other pre-treatment periods, and the units that adopt later or
never, before the cohort adopts), so data from later periods do not affect the placebos. Requires
\code{time.component.from = "notyettreated"} (the default); with \code{"nevertreated"}, \code{fect} stops.}
\item{dloo.adjust}{a logical flag, used only with \code{dloo = TRUE}. When
\code{TRUE}, the pre-trend estimates use the pre-treatment-average baseline of
Liu (2025) instead of the leave-one-out baseline -- equivalently, the
double-LOO rescaled by \eqn{(g-2)/(g-1)} for a cohort adopting at period
\eqn{g}.}
\item{permute}{a logical flag indicating whether to run a permutation test.}
\item{m}{an integer specifying the block length for permutation tests. Default \code{2}.}
\item{normalize}{a logical flag indicating whether to divide the outcome and covariates by the standard deviation of the outcome during estimation. Default \code{FALSE}. Every slot of the fit is reported on the outcome's scale, so the option changes only the speed of estimation, and a given \code{lambda} (\code{method = "mc"}) is on the outcome's scale too. Before version 2.4.7, \code{Y.dat}, \code{est.cm}, the stored covariates, \code{data.long}, \code{sigma2.fect} (and so the default equivalence threshold), \code{est$VNT}, the gsynth \code{IC}, the cross-validation tables and the mc \code{lambda} slots were on the normalized scale or put back with the wrong power of sd(Y).}
\item{keep.sims}{a logical flag indicating whether to save unit-time level bootstrap effects (\code{eff.boot}, with \code{D.boot}, \code{I.boot} and \code{colnames.boot}; see Value). Default \code{keep.sims = FALSE}. If \code{se = FALSE}, this argument is ignored. \code{\link{effect}}, \code{imputed_outcomes(replicates = TRUE)} and most \code{\link{estimand}} calls need these draws.}
\item{cm}{a logical flag indicating whether to enable causal moderation analysis. When \code{TRUE}, the estimator decomposes the treatment effect into effect modification and causal moderation components. Currently available for \code{method = "fe"} and \code{method = "ife"} with \code{time.component.from = "notyettreated"}. The fit then carries \code{est.cm}, the model of the treated outcomes fitted on the treated cells, which \code{\link{fect_iden}} and \code{plot(type = "hte", cm = TRUE)} use; it is stored with or without standard errors and after cross-validation. Default is \code{FALSE}.}
\item{loading.bound}{a string controlling whether treated-unit factor loadings are bounded inside the convex hull of control loadings. \code{"none"} (default) reproduces standard GSC behavior. \code{"simplex"} constrains each treated unit's loading to be a non-negative convex combination of control loadings via an entropy-regularized simplex projection, ensuring the imputed counterfactual lies pointwise in the convex hull of factor-implied control outcomes. Currently applies only to \code{method = "ife"} or \code{method = "gsynth"} (equivalent forms), and requires \code{time.component.from = "nevertreated"}. The bound applies wherever the model is fitted: with or without cross-validation, with or without \code{se}, in the bootstrap and parametric-bootstrap replicates, and in the \code{loo} refits. With a unit effect (\code{force = "unit"} or \code{"two-way"}), each treated unit's intercept is fitted jointly with its constrained loadings (since 2.4.8); the intercept itself is not bounded. Fits with \code{force = "none"} or \code{"time"} have no intercept to fit. If the selected number of factors is 0, a message says that the bound has no effect.}
\item{gamma.loading}{scalar regularization strength for the \code{"simplex"} projection. \code{NULL} (default) triggers 5-fold cross-validation over \code{gamma.loading.grid}. A numeric value is used directly. Ignored when \code{loading.bound = "none"}.}
\item{gamma.loading.grid}{a numeric vector of candidate \code{gamma.loading} values for cross-validation. \code{NULL} (default) uses \code{10^seq(-2, 2, length.out = 9)}. Ignored when \code{loading.bound = "none"} or when \code{gamma.loading} is supplied.}
\item{cv.rule}{a string selecting the cross-validation rule for choosing the number of factors \code{r} (or matrix-completion penalty \code{lambda}). One of:
\describe{
  \item{\code{"min"} (default)}{Pick the \code{r} that minimizes the mean CV criterion (no tolerance). The default since 2.4.7: the fold errors are heavy-tailed, which inflates the fold SE, so the 1-SE rule under-selects \code{r} (in a simulation study across 48 panel designs it picked the true \code{r} 10 to 15 points less often than the min rule and gave a larger ATT error in every design). Fits made with 2.3.0 to 2.4.6 used \code{"1se"}.}
  \item{\code{"1se"}}{The 1-SE rule (Breiman, Friedman, Olshen and Stone 1984; Hastie, Tibshirani and Friedman 2009, Section 7.10): pick the smallest \code{r} whose mean CV criterion is within one fold-SE of the minimum-CV-error \code{r}. Biases toward parsimony in a fold-aware way --- when CV is precise, it allows larger \code{r}; when CV is noisy, it gravitates to simpler models.}
  \item{\code{"1pct"}}{Legacy pre-2.3.0 heuristic: pick the smallest \code{r} within 1\% relative tolerance of the best mean CV criterion. Use this for byte-identical reproducibility of pre-2.3.0 fits.}
}
Ignored when \code{CV = FALSE}.}
}
\details{
  \code{fect} implements counterfactual estimators for TSCS data. It first imputes counterfactuals by fitting an outcome model using untreated observations, then estimates the individual treatment effect as the difference between observed and predicted outcomes. Finally, it computes average treatment effects on the treated (ATT) and period-specific ATTs. Placebo and equivalence tests help evaluate identification assumptions.
}
\value{
  An object of class \code{"fect"}: a list. Which slots it holds depends on the method and the options; a slot that a fit does not need is absent or \code{NULL}. Matrices called T-by-N have one row per period and one column per unit of the fit (units removed before estimation are not among them, except in \code{obs.missing}). "Event time" is time relative to treatment onset: 1 is the first treated period and 0 the last period before it.

  \emph{Point estimates of the effects.}
  \item{att.avg}{the overall ATT: the average of \code{eff} over the treated cells.}
  \item{att.avg.unit}{the average over treated units of each unit's average effect over its observed treated cells ("Tr units equally weighted" in \code{print()}).}
  \item{att, time, count}{the ATT at each event time, the event times, and the number of treated cells at each.}
  \item{att.off, time.off, count.off}{the same by time relative to treatment exit (1 is the first period after treatment ends); only when treatment switches off for some unit.}
  \item{eff}{T-by-N matrix of cell-level effects, the observed outcome minus \code{Y.ct}; \code{NA} where the outcome is missing.}
  \item{Y.ct}{T-by-N matrix of predicted untreated outcomes; \code{NA} where the outcome is missing.}
  \item{eff.calendar, eff.calendar.fit}{two-column matrices by calendar period: the ATT of the treated cells in that period (\code{ATT-calendar}), or its loess fit (\code{ATT-calendar Fitted}), and their number (\code{count}).}
  \item{N.calendar}{the number of treated cells in each calendar period.}
  \item{calendar.enp}{the equivalent number of parameters of the loess fit in \code{eff.calendar.fit}.}
  \item{att.placebo, att.carryover}{with \code{placeboTest} or \code{carryoverTest}: the ATT of the held-out cells.}
  \item{carry.att, carry.time}{with \code{carryover.rm}: the ATT of the cells marked by \code{carryover.rm}, by time after treatment exit.}
  \item{group.att, group.output}{with \code{group}: the ATT of each group, and for each group a list of its \code{att.on}, \code{time.on}, \code{count.on} (by event time) and, with reversals, \code{att.off}, \code{time.off}, \code{count.off}.}
  \item{group, g.level, G}{with \code{group}: a data frame matching each value of the group variable (\code{rawgroup}) to its code (\code{newgroup}); the codes, named by group value; and the T-by-N matrix of codes.}
  \item{balance.att, balance.time, balance.count, balance.avg.att}{with \code{balance.period}: the ATT by event time in the balanced sample, the event times, the numbers of cells, and the overall ATT of the balanced sample. \code{balance.att.placebo} is added with \code{placeboTest}.}

  With \code{W} (or \code{W.agg}), \code{att.avg}, \code{att}, \code{att.off}, \code{att.placebo}, \code{att.carryover} and their inference slots (\code{est.avg}, \code{est.att}, \code{est.att.off}, ...) hold the weighted estimates; a fit has no separate weighted slots.

  \emph{Inference} (\code{se = TRUE}). Unless noted, intervals are (1 - \code{alpha}) intervals built by \code{ci.method}, and p-values test that the effect is 0.
  \item{est.avg}{1-row matrix for \code{att.avg} with columns \code{ATT.avg}, \code{S.E.}, \code{CI.lower}, \code{CI.upper}, \code{p.value}.}
  \item{est.avg.unit}{the same for \code{att.avg.unit} (first column \code{ATT.avg.unit}).}
  \item{est.att}{matrix for \code{att}, one row per event time (the row names): columns \code{ATT}, \code{S.E.}, \code{CI.lower}, \code{CI.upper}, \code{p.value}, \code{count}.}
  \item{est.att90, att.bound}{\code{est.att} with the (1 - 2\code{alpha}) interval (90\% by default) in place of the (1 - \code{alpha}) one, and that interval alone (\code{CI.lower}, \code{CI.upper}). The equivalence plot (\code{plot(type = "equiv")}) draws it.}
  \item{att.vcov}{the covariance matrix of the bootstrap draws of \code{att} across event times (rows and columns in the order of \code{time}); \code{\link{fect_sens}} reads it.}
  \item{est.att.off, att.off.bound, att.off.vcov}{the same three for \code{att.off} (first column \code{ATT.OFF}, last \code{count.off}).}
  \item{est.eff.calendar, est.eff.calendar.fit}{inference for \code{eff.calendar} and \code{eff.calendar.fit}, one row per calendar period, with the columns of \code{est.att} (the first named \code{ATT-calendar} or \code{ATT-calendar Fitted}).}
  \item{est.placebo, est.carryover}{with \code{placeboTest} or \code{carryoverTest}: a 1-row matrix with the ATT of the held-out cells (column \code{ATT.placebo} or \code{ATT.carryover}; named \code{Coef} when the model has covariates), \code{S.E.}, \code{CI.lower}, \code{CI.upper}, \code{p.value}, and the (1 - 2\code{alpha}) interval in \code{CI.lower(90\%)} and \code{CI.upper(90\%)}.}
  \item{est.beta}{inference for \code{beta}: columns \code{Coef}, \code{S.E.}, \code{CI.lower}, \code{CI.upper}, \code{p.value}; \code{NA} rows for dropped covariates.}
  \item{est.group.att, est.group.output}{with \code{group}: a matrix with one row per group (columns as in \code{est.avg}, first \code{ATT}), and for each group a list with \code{att.on} (as \code{est.att}), \code{att.on.bound}, \code{att.on.boot}, \code{att.off}, \code{att.off.bound}, \code{att.off.boot}, \code{att.placebo} and \code{att.carryover}. The interval and p-value follow \code{ci.method}.}
  \item{est.balance.att, est.balance.avg, est.balance.placebo, balance.att.bound, balance.att.vcov}{with \code{balance.period}: \code{est.att}, \code{est.avg}, \code{est.placebo}, \code{att.bound} and \code{att.vcov} for the balanced sample.}
  \item{Y.avg}{with \code{vartype = "bootstrap"} or \code{"jackknife"}: a data frame of the average observed outcome of the treated units (\code{treated}) and their average predicted untreated outcome (\code{counterfactual}) by \code{period}, with 95\% intervals (\code{lower.tr}, \code{upper.tr}, \code{lower.ct}, \code{upper.ct}) and 90\% intervals (\code{lower90.tr}, ...). The levels are fixed, whatever \code{alpha} and \code{ci.method}: basic intervals for the bootstrap, normal intervals for the jackknife. \code{period} is event time when the treated units start treatment at different times (attribute \code{period_type = "event_time"}) and calendar time when they all start at once.}
  \item{vartype}{the variance estimator used: \code{"bootstrap"}, \code{"jackknife"} or \code{"parametric"}.}
  \item{para.error}{for \code{vartype = "parametric"}: the residual model used, with \code{"auto"} resolved to \code{"empirical"} or \code{"ar"}.}
  \item{ci.method, ci.alpha}{the interval rule (\code{"normal"} or \code{"basic"}) and \code{alpha}. The gap plot with \code{id} uses them for the band of the chosen units; see \code{\link{plot.fect}}.}

  \emph{Bootstrap draws} (\code{se = TRUE}). One column per kept replicate; failed replicates are dropped.
  \item{att.avg.boot, att.avg.unit.boot}{draws of \code{att.avg} and \code{att.avg.unit}.}
  \item{att.boot, att.count.boot}{draws of \code{att} and of \code{count}, one row per event time.}
  \item{att.boot.original}{\code{att.boot} before failed replicates are dropped.}
  \item{att.off.boot, att.off.count.boot}{draws of \code{att.off} and \code{count.off}.}
  \item{att.placebo.boot, att.carryover.boot}{draws of \code{att.placebo} and \code{att.carryover}.}
  \item{balance.att.boot, balance.count.boot}{draws of \code{balance.att} and \code{balance.count}.}
  \item{beta.boot}{draws of \code{beta} (\code{NA} rows for dropped covariates).}
  \item{eff.boot}{with \code{keep.sims = TRUE}, the cell-level treatment effects of each replicate, a T-by-W-by-B array (B = number of replicates kept). Slice \code{b} holds the units \code{colnames.boot[[b]]}, in that order: under the case bootstrap a unit drawn twice appears twice and an undrawn unit not at all; under the parametric bootstrap the treated units come first. W = N, except for the cluster bootstrap (\code{cl}), whose replicates differ in size and are padded with \code{NA} to the widest one, and the jackknife (W = N - 1, B = N).}
  \item{D.boot, I.boot}{with \code{keep.sims = TRUE}, the treatment and observation indicators of each replicate, laid out as \code{eff.boot}.}
  \item{colnames.boot}{with \code{keep.sims = TRUE}, a list with one element per replicate: the column positions, in the original panel, of the replicate's units; its length is the replicate's width.}
  \item{eff.off.boot}{with \code{keep.sims = TRUE} and reversals, a list with one \code{eff.off}-like matrix per replicate.}

  \emph{Tests of the identifying assumptions.}
  \item{test.out}{(\code{se = TRUE}) a list. Without a placebo or carryover test, it tests the pre-treatment ATTs in the event times \code{pre.periods}: \code{f.stat}, \code{f.p} (F test that they are all 0), \code{f.threshold}, \code{f.equiv.p} (equivalence F test), \code{df1}, \code{df2}, \code{N_bar} (the largest number of treated units in those event times), \code{tost.equiv.p} (the largest p-value of the two one-sided tests, one per event time) and \code{tost.threshold}. With \code{placeboTest}: \code{placebo.p} and \code{placebo.equiv.p} (TOST of the placebo ATT); with \code{carryoverTest}: \code{carryover.p} and \code{carryover.equiv.p}. An F statistic that cannot be computed (for example a singular covariance of the draws) is \code{NA}, with a message.}
  \item{loo.test.out, dloo.test.out}{with \code{loo} or \code{dloo}: the same tests on the leave-one-out (or double leave-one-out) pre-treatment estimates.}
  \item{pre.est.att, pre.att.bound, pre.att.boot, pre.est.group.output}{when \code{loo} or \code{dloo} is \code{TRUE}, the pre-treatment placebo estimates by event time (columns \code{ATT}, \code{S.E.}, \code{CI.lower}, \code{CI.upper}, \code{p.value}, \code{count.on}), their (1 - 2\code{alpha}) intervals, their draws, and by group with \code{group}. For \code{dloo} these are the double-LOO (or, with \code{dloo.adjust}, pre-treatment-average) estimates.}
  \item{permute}{with \code{permute = TRUE}: a list with \code{permute.att.avg}, the overall ATT of each permutation, and \code{p}, the permutation p-value.}
  \item{eff.pre, eff.pre.equiv}{matrices of the treated units' pre-treatment cells: the effect (\code{eff}), or the effect under the two-way FE model (\code{eff.equiv}), its event time (\code{period}) and the unit's column position (\code{unit}).}
  \item{pre.sd}{by event time before treatment: the standard deviation of \code{eff.pre.equiv} (\code{sd}), the event time (\code{period}) and the number of cells (\code{count}).}
  \item{eff.off, eff.off.equiv, off.sd}{the same three by time after treatment exit; only with reversals.}
  \item{rmse}{the root mean squared \code{eff} over the treated units' untreated cells used in the fit.}
  \item{sigma2.fect}{the residual variance of the two-way FE model; the default \code{tost.threshold} is \code{0.36 * sqrt(sigma2.fect)}.}
  \item{tost.threshold, proportion, pre.periods}{the equivalence range used (with \code{NULL} resolved), \code{proportion}, and the event times the tests use (resolved from \code{proportion} when \code{pre.periods = NULL}).}

  \emph{The fitted model.}
  \item{method}{the method as the fit records it, which can differ from the \code{method} argument: a fit without factors can report \code{"fe"}; a \code{method = "cfe"} fit with not-yet-treated controls reports \code{"fe"} (\code{r = 0}) or \code{"ife"}; a fit whose factors are estimated from the never-treated units reports \code{"gsynth"} (\code{"cfe"} for \code{method = "cfe"}); a \code{"both"} fit reports the method chosen. See the option tables in the cheatsheet chapter of the user manual.}
  \item{r.cv}{the number of factors (selected by cross-validation if needed).}
  \item{lambda.cv, lambda.seq, lambda.norm, eigen.all}{\code{method = "mc"}: the chosen penalty, the candidates, the chosen penalty divided by the largest of \code{eigen.all}, and the singular values (divided by T x N) of the untreated outcome matrix after the fixed effects, from which the candidates are built.}
  \item{CV.out, CV.out.ife, CV.out.mc}{results of the cross-validation procedure: one row per candidate \code{r} (or \code{lambda}), the scores of each; a row the loop did not score holds \code{NA}. An \code{ife} fit also keeps its table as \code{CV.out.ife}, an \code{mc} fit as \code{CV.out.mc}, and a \code{"both"} fit keeps both tables; \code{CV.out} is the one of the method fitted.}
  \item{CV.out.se, CV.out.ife.se, CV.out.mc.se}{the fold standard errors of the scores in \code{CV.out} (same rows, the score columns), the ones \code{cv.rule} reads; with \code{CV.out.ife} and \code{CV.out.mc} come \code{CV.out.ife.se} and \code{CV.out.mc.se}.}
  \item{rmCV, estCV}{cross-validation: for each fold, the positions of the held-out cells and of the cells scored; kept for internal use.}
  \item{beta}{coefficients of the covariates (column \code{Coef}), one row per requested covariate; \code{NA} for a covariate dropped before the fit (see the \code{X} argument).}
  \item{mu, alpha, xi}{the grand mean, the unit effects (N-by-1) and the time effects (T-by-1), as \code{force} includes them. Never-treated fits keep the unit effects as \code{alpha.tr} and \code{alpha.co} only.}
  \item{alpha.tr, alpha.co}{the unit effects of the treated and the control units. With \code{loading.bound = "simplex"}, \code{alpha.tr} comes from the joint fit with the constrained loadings (since 2.4.8).}
  \item{factor, lambda, lambda.tr, lambda.co}{with factors: the factors (T-by-r) and the loadings of all units (not kept by never-treated fits), of the treated units and of the control units.}
  \item{wgt.implied}{for fits whose factors are estimated from the never-treated units (\code{method = "gsynth"}, or \code{"ife"} and \code{"cfe"} with \code{time.component.from = "nevertreated"}) with at least one factor: the implied weights of the control units for each treated unit (Xu 2017), an Nco x Ntr matrix with the control unit IDs as row names and the treated unit IDs as column names. Column \eqn{i} holds \eqn{w_i = \Lambda_{co} (\Lambda_{co}' \Lambda_{co})^{-1} \lambda_i}{w_i = L_co (L_co' L_co)^(-1) lambda_i}, which satisfies \code{t(lambda.co) \%*\% w_i == lambda.tr[i, ]}; so the factor part of treated unit \eqn{i}'s counterfactual is the \eqn{w_i}-weighted sum of the control units' factor parts. The weights can be negative. Under \code{loading.bound = "simplex"}, column \eqn{i} holds the simplex weights of treated unit \eqn{i}: non-negative, summing to one. \code{NULL} when no factor is estimated.}
  \item{loading.bound, gamma.loading, loading.proj.resid}{for never-treated fits with factors: the bound used; under \code{"simplex"} also the \code{gamma.loading} used and, for each treated unit, the norm of the residual of its constrained fit.}
  \item{gamma, kappa}{\code{method = "cfe"} (and never-treated fits): the estimated coefficients on \code{Z} and loadings on \code{Q}; \code{NULL} without them.}
  \item{est}{the fitted model as the estimator returns it: a list with \code{beta}, \code{mu}, \code{fit} (fitted values), \code{niter}, \code{alpha}, \code{xi}, \code{residuals}, \code{sigma2}, \code{validX}, and, depending on the method, \code{IC}, \code{PC}, \code{validF}, \code{factor}, \code{lambda} and \code{VNT}.}
  \item{est.cm}{with \code{cm = TRUE}: the fitted model of the treated outcomes (see \code{cm}).}
  \item{sigma2}{the residual variance of the fitted model.}
  \item{IC, PC}{information criteria of the fitted model (not for \code{"mc"}); \code{criterion = "pc"} selects on \code{PC}.}
  \item{niter}{the number of iterations of the final fit. A warning is given when it reaches \code{max.iteration}.}
  \item{validX, validF}{1 if covariates (0 when every covariate was dropped) or factors are in the model.}
  \item{res, res.full}{T-by-N matrices of residuals: at the cells used to fit the model (\code{NA} elsewhere), and as the estimator returns them, before the other cells are set to \code{NA}.}
  \item{Y.ct.full}{\code{Y.ct} before the cells with a missing outcome are set to \code{NA}.}
  \item{eff.tr}{\code{eff} for the treated units (T-by-Ntr); not kept by every method.}
  \item{Y.tr.cnt, Y.ct.cnt}{never-treated fits: the predicted untreated outcomes of the treated units (T-by-Ntr) and the fitted outcomes of the control units (T-by-Nco).}

  \emph{Data, sample and settings.}
  \item{Y.dat, D.dat, I.dat}{T-by-N matrices of the outcome, the treatment and the observation indicator (1 if the outcome and covariates are observed).}
  \item{Y, D, X}{the names of the outcome, the treatment and the requested covariates (including any dropped before the fit). Each name occurs twice in the list: the second \code{Y} and \code{D} are the estimator's T-by-N outcome and treatment matrices, and the second \code{X} is its covariate array, T-by-N-by-p (\code{x[names(x) == "X"][[2]]}).}
  \item{I, II}{T-by-N indicators: the cell is observed (\code{I}); the cell was used to fit the model, that is, observed, untreated and not held out for a test (\code{II}).}
  \item{T.on, T.off, T.on.balance}{T-by-N matrices of event time, of time relative to treatment exit (with reversals), and of event time within the balanced sample (with \code{balance.period}); \code{NA} where not defined.}
  \item{obs.missing}{a T x N matrix over all units, including those removed, coding each cell: 1 under treatment; 2 not under treatment (control units, and treated units before treatment); 3 missing; 4 in a unit removed from the sample; 5 held out for the placebo or carryover test; 6 removed by \code{carryover.rm}.}
  \item{obs.missing.balance}{with \code{balance.period}: \code{obs.missing} with the cells of the balanced sample coded 7 (after onset) and 8 (before); \code{NA} otherwise.}
  \item{sample}{a T x N logical matrix over all units: \code{TRUE} for the cells coded 1, 2 or 5 in \code{obs.missing}.}
  \item{unit.type}{for each unit of the fit: 1 untreated in all its observed periods, 2 treated in all, 3 treatment status changes.}
  \item{tr, co, Ntr, Nco}{the column positions of the treated (ever-treated) and the never-treated units, and their numbers.}
  \item{T, N, p}{the numbers of periods, units in the fit and covariates.}
  \item{id}{the IDs of the units in the fit.}
  \item{rawtime}{the periods; for a factor time index whose levels are not increasing numbers, the level labels.}
  \item{remove.id}{IDs of the units removed before estimation (too few untreated periods, see \code{min.T0}, or too many missing periods, see \code{max.missing}). Present only when units were removed.}
  \item{hasRevs}{1 if treatment switches off for some unit.}
  \item{data.long}{the input panel (unit, time, outcome and treatment columns) before units were removed; \code{panelView::panelview(fit)} draws it.}
  \item{index}{the unit and time column names (with \code{group.fe} appended).}
  \item{force}{the fixed effects as a code: 0 \code{"none"}, 1 \code{"unit"}, 2 \code{"time"}, 3 \code{"two-way"}.}
  \item{W}{the name of the weight column. A weighted fit has a second element \code{W}, its T-by-N weight matrix.}
  \item{W.agg}{with \code{W} or \code{W.agg}: the aggregation weights as a T-by-N matrix on the fit's panel (like \code{Y.dat}); read by \code{\link{imputed_outcomes}}, \code{\link{estimand}} and \code{\link{effect}}. Absent for unweighted fits and with \code{W.est} alone.}
  \item{W.est.col, W.agg.col, W.in.fit, W.in.agg}{the weight column used in the fit and in the aggregation (\code{NULL} if none), and whether each is weighted.}
  \item{time.component.from, em}{the values used (after the changes \code{method = "gsynth"} makes).}
  \item{placeboTest, placebo.period, carryoverTest, carryover.period, carryover.rm, balance.period, loo, dloo, dloo.adjust, binary, group.fe, cl.label}{the settings, as used (\code{loo = TRUE} turns off \code{placeboTest} and \code{carryoverTest}); \code{cl.label} is the cluster column name.}
  \item{formula, call}{the formula (formula interface only) and the call.}
  \item{eff_debias}{reserved; \code{NULL}.}
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

  Liu, L., Wang, Y., and Xu, Y. (2024).
  A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional Data.
  \emph{American Journal of Political Science, 68}(1), 160-176.

  Xu, Y. (2017).
  Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models.
  \emph{Political Analysis, 25}(1), 57-76.
}
\seealso{
  \code{\link{plot.fect}}, \code{\link{print.fect}}, \code{\link{estimand}},
  \code{\link{esplot}}, \code{\link{get.cohort}}, \code{\link{did_wrapper}},
  \code{\link{fect_sens}}, \code{\link{fect_mspe}}
}
\examples{
library(fect)
data(simdata)
out <- fect(Y ~ D + X1 + X2, data = simdata,
            index = c("id","time"), force = "two-way",
            CV = TRUE, r = c(0, 5), se = 0, parallel = FALSE)
}
