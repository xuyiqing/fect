\name{fect}
\alias{fect}
\title{Fixed Effects Counterfactual Estimators}
\description{Implements counterfactual estimators in TSCS data analysis and statistical tools to test their identification 
assumptions.}
\usage{fect(formula = NULL, data, Y, D, X = NULL, group = NULL,
            na.rm = FALSE, 
            index, force = "two-way", r = 0, lambda = NULL, nlambda = 10,
            CV = NULL, k = 10, cv.prop = 0.1, cv.treat = FALSE, 
            cv.nobs = 3, cv.donut = 0, criterion = "mspe",
            binary = FALSE, QR = FALSE,
            method = "fe",  
            se = FALSE, vartype = "bootstrap", nboots = 200, alpha = 0.05,
            parallel = TRUE, cores = NULL, tol = 0.001, seed = NULL, 
            min.T0 = NULL, max.missing = NULL, 
            proportion = 0.3, pre.periods = NULL, 
            f.threshold = 0.5, tost.threshold = NULL,
            knots = NULL, degree = 2, 
            sfe = NULL, cfe = NULL,
            balance.period = NULL, fill.missing = FALSE,
            placeboTest = FALSE, placebo.period = NULL,
            carryoverTest = FALSE, carryover.period = NULL, carryover.rm = NULL,
            loo = FALSE, permute = FALSE, m = 2, normalize = FALSE)  
}
\arguments{
\item{formula}{an object of class "formula": a symbolic description of the model to be fitted, e.g, Y~D+X1+X2}
\item{data}{a data frame, can be a balanced or unbalanced panel data.}
\item{Y}{the outcome indicator.}
\item{D}{the treatment indicator. The treatment should be binary (0 and 1).}
\item{X}{time-varying covariates. Covariates that have perfect collinearity with specified fixed effects are dropped automatically.}
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
\item{cv.donut}{an integer specifying the length of removed observations at the head and tail of the continuous observations specified by \code{cv.nobs}. 
  These removed observations will not be used to fit the data nor be in the validation set for the cross-validation, e.g, if \code{cv.nobs=3} and \code{cv.donut = 1}, the first and the last observation in each triplet will not be included in the test set. 
  Default is \code{cv.donut = 0}.}
\item{criterion}{criterion used for model selection. Default is "mspe". 
  \code{"mspe"} for the mean squared prediction error,
  \code{"gmspe"} for the geometric-mean squared prediction errors,
  if \code{criterion="moment"}, we average the residuals in test sets by their relative periods to treatments and then average the squares of these period-wise deviations weighted by the number of observations at each period, it yields a better pre-trend fitting on test sets rather than a better prediction ability.
  \code{"pc"} for the information criterion of interactive fixed effects or generalized synthetic control model.}
\item{binary}{This version doesn't support this option. a logical flag indicating whether a probit link function will be used.}
\item{QR}{This version doesn't support this option. a logical flag indicating whether QR decomposition will be used for factor analysis in probit model.}
\item{method}{a string specifying which imputation algorithm will be used. 
  \code{"fe"} for fixed effects model, 
  \code{"ife"} for interactive fixed effects model, 
  \code{"mc"} for matrix copletion method, 
  \code{"polynomial"} for polynomial trend terms, 
  \code{"bspline"} for regression splines,
  \code{"gsynth"} for generalized synthetic control method,
  and \code{"cfe"} for complex fixed effects method 
  Default is \code{method = "fe"}.}
\item{se}{a logical flag indicating whether uncertainty estimates will
  be produced.}
\item{vartype}{a string specifying the type of variance estimator. Choose 
  from \code{vartype = c("bootstrap", "jackknife", "parametric")}. Default value is 
  \code{"bootstrap"}.}
\item{nboots}{an integer specifying the number of bootstrap
  runs. Ignored if \code{se = FALSE}.}
\item{alpha}{significant level for hypothesis test and CIs. Default value is 
  \code{alpha = 0.05}.}
\item{parallel}{a logical flag indicating whether parallel computing
  will be used in bootstrapping and/or cross-validation. Ignored if
  \code{se = FALSE}.}
\item{cores}{an integer indicating the number of cores to be used in
  parallel computing. If not specified, the algorithm will use the
  maximum number of logical cores of your computer (warning: this
  could prevent you from multi-tasking on your computer).}
\item{tol}{a positive number indicating the tolerance level.}
\item{seed}{an integer that sets the seed in random number
  generation. Ignored if \code{se = FALSE} and \code{r} is specified.}
\item{min.T0}{an integer specifying the minimum value of observed periods 
  that a unit is under control.}
\item{max.missing}{an integer. Units with number of missing values greater than 
  it will be removed. Ignored if this parameter is set "NULL"(i.e. \code{max.missing = NULL}, the default setting).}
\item{proportion}{a numeric value specifying pre-treatment periods that have 
  observations larger than the proportion of observations at period 0. 
  These pre-treatment periods are used used for goodness-of-fit test.
  Ignore if \code{se = FALSE}. Deafult is \code{proportion = 0.3}.}
\item{pre.periods}{a vector specifying the range of pre-treatment period used for goodness-of-fit test. If left blank, all 
  pre-treatment periods specified by \code{proportion} will be used. Ignore if \code{se = FALSE}.}
\item{f.threshold}{a numeric value specifying the threshold for the F-statistic in the equivalent test.
  Ignore if \code{se = FALSE}.  Deafult is \code{f.threshold = 0.5}.}
\item{tost.threshold}{a numeric value specifying the threshold for the two-one-sided t-test.
  If \code{alpha=0.05}, TOST checks whether the 90% confidence intervals for the estimated period-wise ATT exceeds this threshold.
  The default value is 0.36 times the standard deviation of the outcome variable after two-way fixed effects are partialed out.}
\item{knots}{a numeric vector speicfying the knots for b-spline curve trend term.}
\item{degree}{an integer speifcying the order of either the b-spline or the 
  polynomial trend term.}
\item{sfe}{a vector specifying other fixed effects in addition to unit or time fixed effects that is used when \code{method="cfe"}.}
\item{cfe}{a vector of lists specifying interactive fixed effects when \code{method="cfe"}. 
  For each list, the value of the first element is the name of the group variable for which fixed effects are to be estimated.
  The value of the second element is the name of a regressor (e.g., a time trend).}
\item{balance.period}{a vector of length 2 specifying the range of periods for a balanced sample which has no missing observation in the specified range.}
\item{fill.missing}{a logical flag indicating whether to allow missing observations in this balanced sample. The default is FALSE.}
\item{placeboTest}{a logic flag indicating whether to perform placebo test.}
\item{placebo.period}{an integer or a two-element numeric vector specifying the 
  range of pre-treatment periods that will be assigned as pseudo treatment periods.}
\item{carryoverTest}{a logic flag indicating whether to perform (no) carryover test.}
\item{carryover.period}{an integer or a two-element numeric vector specifying the 
  range of post-treatment periods that will be assigned as pseudo treatment periods.}
\item{carryover.rm}{an integer specifying the range of post-treatment periods that will be assigned as pseudo treatment periods.}
\item{loo}{a logic flag indicating whether to perform the leave-one-period-out goodness-of-fit test, which is very time-consuming.}
\item{permute}{a logic flag indicating whether to perform permutation test.}
\item{m}{an integer specifying the block length in permutation test. Default value is 
  \code{m = 2}.}
\item{normalize}{a logic flag indicating whether to scale outcome and 
  covariates. Useful for accelerating computing speed when magnitude of data is large. The default is \code{normalize=FALSE}.}
}
\details{
  \code{fect} implements counterfactual estimators in TSCS data analysis. These estimators first impute counterfactuals for 
  each treated observation in a TSCS dataset by fitting an outcome model (fixed effects model, interactive fixed effects model, or
  matrix completion) using the untreated observations. They then estimate the individualistic treatment effect for each treated 
  observation by subtracting the predicted counterfactual outcome from its observed outcome. Finally, the average treatment effect
  on the treated (ATT) or period-specific ATTs are calculated. A placebo test and an equivalence test are included to evaluate the
  validity of identification assumptions behind these estimators. Data must be with a dichotomous treatment.
}
\value{
  \item{Y.dat}{a T-by-N matrix storing data of the outcome variable.}
  \item{D.dat}{a T-by-N matrix storing data of the treatment variable.}
  \item{I.dat}{a T-by-N matrix storing data of the indicator for whether is observed or missing.}
  \item{Y}{name of the outcome variable.}
  \item{D}{name of the treatment variable.}
  \item{X}{name of the time-varying control variables.}
  \item{index}{name of the unit and time indicators.}
  \item{force}{user specified \code{force} option.}
  \item{T}{the number of time periods.}
  \item{N}{the total number of units.}
  \item{p}{the number of time-varying observables.}
  \item{r.cv}{the number of factors included in the model -- either
    supplied by users or automatically chosen via cross-validation.}
  \item{lambda.cv}{the optimal hyper-parameter in matrix completion method
  chosen via cross-validation.}
  \item{beta}{coefficients of time-varying observables from the
    interactive fixed effect model.}
  \item{sigma2}{the mean squared error of interactive fixed effect model.}
  \item{IC}{the information criterion.}
  \item{est}{result of the interactive fixed effect model based on
    observed values.}
  \item{MSPE}{mean squared prediction error of the cross-validated model.}
  \item{CV.out}{result of the cross-validation procedure.}
  \item{niter}{the number of iterations in the estimation of the
    interactive fixed effect model.}
  \item{factor}{estimated time-varying factors.}
  \item{lambda}{estimated loadings.}
  \item{lambda.tr}{estimated loadings for treated units.}
  \item{lambda.co}{estimated loadings for control units.}
  \item{mu}{estimated ground mean. }
  \item{xi}{estimated time fixed effects.}
  \item{alpha}{estimated unit fixed effects.}
  \item{alpha.tr}{estimated unit fixed effects for treated units.}
  \item{alpha.co}{estimated unit fixed effects for control units.}
  \item{validX}{a logic value indicating if multicollinearity exists.}
  \item{validF}{a logic value indicating if factor exists.}
  \item{id}{a vector of unit IDs.}
  \item{rawtime}{a vector of time periods.}
  \item{obs.missing}{a matrix stroing status of each unit at each time point.}
  \item{Y.ct}{a T-by-N matrix storing the predicted Y(0).}
  \item{eff}{a T-by-N matrix storing the difference between actual outcome and predicted Y(0).}
  \item{res}{residuals for observed values.}
  \item{eff.pre}{difference between actual outcome and predicted Y(0) for 
    observations of treated units under control.}
  \item{eff.pre.equiv}{difference between actual outcome and predicted Y(0) for 
    observations of treated units under control based on baseline (two-way fixed effects) model.}
  \item{pre.sd}{by period residual standard deviation for estimated pre-treatment 
    average treatment effects.}
  \item{att.avg}{average treatment effect on the treated.}
  \item{att.avg.unit}{by unit average treatment effect on the treated.}
  \item{time}{term for switch-on treatment effect.}
  \item{count}{count of each term for switch-on treatment effect.}
  \item{att}{switch-on treatment effect.}
  \item{time.off}{term for switch-off treatment effect.}
  \item{att.off}{switch-off treatment effect.}
  \item{count.off}{count of each term for switch-off treatment effect.}
  \item{att.placebo}{average treatment effect for placebo period.}
  \item{att.carryover}{average treatment effect for carryover period.}
  \item{eff.calendar}{average treatment effect for each calendar period.}
  \item{eff.calendar.fit}{loess fitted values of average treatment effect for each calendar period.}
  \item{N.calandar}{number of treated observations at each calendar period.}
  \item{balance.avg.att}{average treatment effect for the balance sample.}
  \item{balance.att}{switch-on treatment effect for the balance sample.}
  \item{balance.time}{term of switch-on treatment effect for the balance sample.}
  \item{balance.count}{count of each term for switch-on treatment effect for the balance sample.}
  \item{balance.att.placebo}{average treatment effect for placebo period of the balance sample.}
  \item{group.att}{average treatment effect for different groups.}
  \item{group.output}{a list saving the switch-on treatment effects for different groups.}
  \item{est.att.avg}{inference for \code{att.avg}.}
  \item{est.att.avg.unit}{inference for \code{att.avg.unit}.}
  \item{est.att}{inference for \code{att.on}.}
  \item{est.att.off}{inference for \code{att.off}.}
  \item{est.placebo}{inference for \code{att.placebo}.}
  \item{est.carryover}{inference for \code{att.carryover}.}
  \item{est.eff.calendar}{inference for \code{eff.calendar}.}
  \item{est.eff.calendar.fit}{inference for \code{eff.calendar.fit}.}
  \item{est.balance.att}{inference for \code{balance.att}.}
  \item{est.balance.avg}{inference for \code{balance.avg.att}.}
  \item{est.balance.placebo}{inference for \code{balance.att.placebo}.}
  \item{est.beta}{inference for \code{beta}.}
  \item{est.group.att}{inference for \code{group.att}.}
  \item{est.group.output}{inference for \code{group.output}.}  
  \item{att.avg.boot}{bootstrap results for \code{att.avg}.}
  \item{att.avg.unit.boot}{bootstrap results for \code{att.avg.unit}.}
  \item{att.count.boot}{bootstrap results for \code{count}.}
  \item{att.off.boot}{bootstrap results for \code{att.avg.off}.}
  \item{att.off.count.boot}{bootstrap results for \code{count.off}.}
  \item{att.placebo.boot}{bootstrap results for \code{att.placebo}.}
  \item{att.carryover.boot}{bootstrap results for \code{att.carryover}.}
  \item{balance.att.boot}{bootstrap results for \code{balance.att}.}
  \item{att.bound}{equivalence confidence interval for equivalence test.}
  \item{att.off.bound}{equivalence confidence interval for equivalence test for 
  switch-off effect.}
  \item{beta.boot}{bootstrap results for \code{beta}.}
  \item{test.out}{goodness-of-fit test and equivalent test results for pre-treatment fitting check.}
  \item{loo.test.out}{leave-one-period-out goodness-of-fit test and equivalent test results for pre-treatment fitting check.}
  \item{permute}{permutation test results for sharp null hypothesis.}
}
\author{
  Licheng Liu; Ye Wang; Yiqing Xu; Ziyi Liu
}
\references{
  Jushan Bai. 2009. "Panel Data Models with Interactive Fixed Effects." Econometrica.
  
  Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models." Political 
  Analysis. 
  
  Athey, Susan, et al. 2021 "Matrix completion methods for causal panel data models." Journal of the American Statistical Association. 
  
  Licheng Liu, et al. 2022. "A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional 
  Data." American Journal of Political Science. 

  For more details about the matrix completion method, see \url{https://github.com/susanathey/MCPanel}. 
}
\seealso{
  \code{\link{plot.fect}} and \code{\link{print.fect}}
}
\examples{
library(fect)
data(fect)
out <- fect(Y ~ D + X1 + X2, data = simdata1, 
            index = c("id","time"), force = "two-way",
            CV = TRUE, r = c(0, 5), se = 0, parallel = FALSE) 
}


