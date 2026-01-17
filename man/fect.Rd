\encoding{UTF-8}
\name{fect}
\alias{fect}
\title{Fixed Effects Counterfactual Estimators}
\description{Implements counterfactual estimators in TSCS data analysis and statistical tools to test their identification
assumptions.}
\usage{fect(formula = NULL, data, Y, D, X = NULL,
            W = NULL, group = NULL,
            na.rm = FALSE,
            index, force = "two-way", r = 0, lambda = NULL, nlambda = 10,
            CV = NULL, k = 10, cv.prop = 0.1, cv.treat = FALSE,
            cv.nobs = 3, cv.donut = 0, criterion = "mspe",
            binary = FALSE, QR = FALSE,
            method = "fe",  se = FALSE, vartype = "bootstrap", cl = NULL,
            quantile.CI = FALSE, nboots = 200, alpha = 0.05,
            parallel = TRUE, cores = NULL, tol = 1e-3,
            max.iteration = 1000, seed = NULL,
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
            loo = FALSE, permute = FALSE, m = 2, normalize = FALSE, keep.sims = FALSE)}
\arguments{
\item{formula}{an object of class "formula": a symbolic description of the model to be fitted, e.g, Y~D+X1+X2}
\item{data}{a data frame, can be a balanced or unbalanced panel data.}
\item{Y}{the outcome indicator.}
\item{D}{the treatment indicator. The treatment should be binary (0 and 1).}
\item{X}{time-varying covariates. Covariates that have perfect collinearity with specified fixed effects are dropped automatically.}
\item{W}{the weight indicator. If specified, the program will fit the data with a weighted outcome model and calculate weighted average and dynamic treatment effect.}
\item{group}{the group indicator. If specified, the group-wise ATT will be estimated.}
\item{na.rm}{a logical flag indicating whether to list-wise delete missing observations. Default to FALSE. If \code{na.rm = FALSE}, it allows the situation when Y is missing but D is not missing for some observations. If \code{na.rm = TRUE}, it will list-wise delete observations whose Y, D, or X is missing.}
\item{index}{a two-element string vector specifying the unit and time indicators. Must be of length 2. Every observation should be uniquely defined by the pair of the unit and time indicator.}
\item{force}{a string indicating whether unit or time or both fixed effects will be imposed. Must be one of the following, "none", "unit", "time", or "two-way". The default is "two-way".}
\item{r}{an integer specifying the number of factors. If \code{CV = TRUE}, the cross validation procedure will select the optimal number of factors from \code{r} to 5.}
\item{lambda}{a single or sequence of positive numbers specifying the hyper-parameter sequence for matrix completion method. If \code{lambda} is a sequence and \code{CV = 1}, cross-validation will be performed.}
\item{nlambda}{an integer specifying the length of hyper-parameter sequence for matrix completion method. Default is \code{nlambda = 10}.}
\item{CV}{a logical flag indicating whether cross-validation will be performed to select the optimal number of factors or hyper-parameter in matrix completion algorithm. If \code{r} is not specified, the procedure will search through \code{r = 0} to \code{5}.}
\item{k}{an integer specifying number of cross-validation rounds. Default is \code{k = 10}.}
\item{cv.prop}{a numerical value specifying the proportion of testing set compared to sample size during the cross-validation procedure.}
\item{cv.treat}{a logical flag speficying whether to only use observations of treated units as testing set.}
\item{cv.nobs}{an integer specifying the length of continuous observations within a unit in the testing set. Default is \code{cv.nobs = 3}.}
\item{cv.donut}{an integer specifying the length of removed observations at the head and tail of the continuous observations specified by \code{cv.nobs}.}
\item{criterion}{criterion used for model selection. Default is "mspe".
  \code{"mspe"} for the mean squared prediction error,
  \code{"gmspe"} for the geometric-mean squared prediction errors,
  \code{"moment"} for period-weighted residuals in test sets,
  \code{"pc"} for an information criterion method.}
\item{binary}{This version doesn't support this option.}
\item{QR}{This version doesn't support this option.}
\item{method}{a string specifying which imputation algorithm will be used.
  \code{"fe"}, \code{"ife"}, \code{"mc"}, \code{"polynomial"}, \code{"bspline"},
  \code{"gsynth"}, or \code{"cfe"}. Default is \code{"fe"}.}
\item{se}{a logical flag indicating whether uncertainty estimates will be produced.}
\item{vartype}{a string specifying the type of variance estimator, e.g. \code{"bootstrap"}.}
\item{cl}{a string specifying the cluster for cluster bootstrapping.}
\item{quantile.CI}{a logical flag indicating whether to use quantile confidence intervals when bootstrapping.}
\item{nboots}{an integer specifying the number of bootstrap runs. Ignored if \code{se=FALSE}.}
\item{alpha}{the significance level for hypothesis tests and confidence intervals. Default \code{0.05}.}
\item{parallel}{a logical flag indicating whether parallel computing will be used in bootstrapping or cross-validation.}
\item{cores}{an integer indicating the number of cores for parallel computing.}
\item{tol}{a positive number indicating the tolerance level for EM updates.}
\item{max.iteration}{the maximal number of iterations for the EM algorithm.}
\item{seed}{an integer seed for random number generation.}
\item{min.T0}{an integer specifying the minimum number of pre-treatment periods for each treated unit.}
\item{max.missing}{an integer specifying the maximum number of missing observations allowed per unit.}
\item{proportion}{a numeric value specifying which pre-treatment periods are used for goodness-of-fit tests.}
\item{pre.periods}{a vector specifying the range of pre-treatment periods used for the goodness-of-fit test.}
\item{f.threshold}{a numeric threshold for an F-test in equivalence testing. Default \code{0.5}.}
\item{tost.threshold}{a numeric threshold for two-one-sided t-tests.}
\item{knots}{a numeric vector specifying knots for the b-spline.}
\item{degree}{an integer specifying the polynomial or b-spline degree.}
\item{sfe}{vector specifying other fixed effects for \code{method="cfe_old"}.}
\item{cfe}{a vector of lists specifying interactive fixed effects for \code{method="cfe_old"}.}
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
data(fect)
out <- fect(Y ~ D + X1 + X2, data = simdata,
            index = c("id","time"), force = "two-way",
            CV = TRUE, r = c(0, 5), se = 0, parallel = FALSE)
}
