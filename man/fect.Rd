\name{fect}
\alias{fect}
\title{Fixed Effects Counterfactual Estimators}
\description{Implements counterfactual estimators in TSCS data analysis and statistical tools to test their identification 
assumptions.}
\usage{fect(formula = NULL, data, Y, D, X = NULL, group = NULL, na.rm = FALSE, 
       index, force = "unit", cl = "unit", r = 0, lambda = NULL, nlambda = 10,
       CV = TRUE, k = 10, cv.prop = 0.1, cv.treat = TRUE, cv.nobs = 3,
       binary = FALSE, QR = FALSE,
       method = "fe", criterion = "mspe", alpha = 0.05,
       se = FALSE, vartype = "bootstrap", nboots = 200, 
       parallel = FALSE, cores = NULL, tol = 0.001, seed = NULL, 
       min.T0 = 5, max.missing = NULL, pre.period = NULL, off.period = NULL, 
       knots = NULL, degree = 2, wald = FALSE,
       placebo.period = NULL, placeboTest = FALSE, 
       permute = FALSE, m = 2, normalize = FALSE)  
}
\arguments{
\item{formula}{an object of class "formula": a symbolic description of
  the model to be fitted.}
\item{data}{a data frame (must be with a dichotomous treatment but balanced
  is not reqiored).}
\item{Y}{outcome.}
\item{D}{treatment.}
\item{X}{time-varying covariates.}
\item{group}{a string indicator of the group variable for average cohort effect.}
\item{na.rm}{a logical flag indicating whether to list-wise delete
  missing data. The algorithm will report an error if missing data
  exist.}
\item{index}{a two-element string vector specifying the unit (group)
    and time indicators. Must be of length 2.}
\item{force}{a string indicating whether unit or time fixed effects will be
    imposed. Must be one of the following,
    "none", "unit", "time", or "two-way". The default is "unit".}
\item{cl}{a string indicator of the cluster variable. The default value is 
  \code{"unit"}. If \code{cl = "unit"}, bootstrap will be blocked at unit level.}
\item{r}{an integer specifying the number of factors. If \code{CV = TRUE}, 
  the cross validation procedure will select the optimal number of
  factors from \code{r} to 5.}
\item{lambda}{a single or sequence of positive numbers specifying the 
  hyper-parameter sequence for matrix completion method. If \code{lambda} 
  is a sequence and \code{CV = 1}, cross-validation will be performed.}
\item{nlambda}{an integer specifying the length of hyper-parameter sequence 
  for matrix completion method. Default is \code{nlambda = 10}.}
\item{CV}{a logical flag indicating whether cross-validation will be
  performed to select the optimal number of factors or hyper-parameter 
  in matrix completion algorithm. If \code{r} is not specified, the procedure will search through \code{r = 0} to \code{5}.}
\item{k}{an integer specifying number of cross-validation rounds. Default is 
  \code{k = 10}. }
\item{cv.prop}{a numerical value specifying the proportion of testing set compared 
  to sample size during the cross-validation procedure.}
\item{cv.treat}{a logical flag speficying whether to only use observations of treated 
  units as testing set.}
\item{cv.nobs}{an integer specifying the length of continuous observations within a 
  unit in the testing set. Default is \code{cv.nobs = 3}.}
\item{binary}{a logical flag indicating whether a probit link function will be 
  used.}
\item{QR}{a logical flag indicating whether QR decomposition will be used for 
  factor analysis in probit model.}
\item{method}{a string specifying which matrix completion algorithm 
  will be used. \code{"fe"} for two-way fixed effects model, 
  \code{"ife"} for interactive fixed effects model, 
  \code{"mc"} for matrix copletion method, \code{"polynomial"} for 
  polynomial trend terms and \code{"bspline"} for regression splines. 
  Default is \code{method = "ife"}.}
\item{criterion}{criterion used for model selection.}
\item{alpha}{significant level for hypothesis test and CIs. Default value is 
  \code{alpha = 0.05}.}
\item{se}{a logical flag indicating whether uncertainty estimates will
  be produced.}
\item{vartype}{a string specifying the type of variance estimator. Choose 
  from \code{vartype = c("bootstrap", "jackknife")}. Default value is 
  \code{"bootstrap"}.}
\item{nboots}{an integer specifying the number of bootstrap
  runs. Ignored if \code{se = FALSE}.}
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
\item{pre.period}{a two-element numeric vector specifying the range of 
  pre-treatment period used for goodness-of-fit test. If left blank, all 
  pre-treatment periods will be used. Ignore if \code{wald = FALSE}.}
\item{off.period}{a two-element numeric vector specifying the range of 
  pre-treatment period used for goodness-of-fit test in situations of 
  treatment reversals. If left blank, all post-treatment periods will be used. 
  Ignore if \code{wald = FALSE}.}
\item{knots}{a numeric vector speicfying the knots for b-spline curve trend term.}
\item{degree}{an integer speifcying the order of either the b-spline or the 
  polynomial trend term.}
\item{wald}{a logic flag indicating whether to perform wald test for pre-treatment 
  fitting check.}
\item{placebo.period}{an integer or a two-element numeric vector specifying the 
  range of pre-treatment period that will be assigned as "placebo" treatment period.}
\item{placeboTest}{a logic flag indicating whether to perform placebo test.}
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
  \item{Y.dat}{a matrix storing data of the outcome variable.}
  \item{D.dat}{a matrix storing data of the treatment variable.}
  \item{I.dat}{a matrix storing data of the indicator for whether is observed or missing.}
  \item{Y}{name of the outcome variable.}
  \item{D}{name of the treatment variable.}
  \item{X}{name of the time-varying control variables.}
  \item{index}{name of the unit and time indicators.}
  \item{id}{a vector of unit IDs.}
  \item{rawtime}{a vector of time periods.}
  \item{unit.type}{a vector storing type for each unit. \code{1} for always 
  under control, \code{2} for always under treatment, and \code{3} for unit 
  that has reversals.}
  \item{obs.missing}{a matrix stroing status of each unit at each time point. 
    \code{1} for under treatment, \code{2} for under control, \code{3} for missing, \code{4} for removed units.}
  \item{Y.ct}{predicted Y(0).}
  \item{eff}{difference between actual outcome and predicted Y(0).}
  \item{eff.pre}{difference between actual outcome and predicted Y(0) for 
    observations of treated units under control.}
  \item{eff.pre.equiv}{difference between actual outcome and predicted Y(0) for 
    observations of treated units under control based on baseline 
    (two-way fixed effects) model.}
  \item{pre.sd}{by period residual standard deviation for estimated pre-treatment 
    average treatment effects.}
  \item{att.avg}{average treatment effect on the treated.}
  \item{att.avg.unit}{by unit average treatment effect on the treated.}
  \item{att.placebo}{average treatment effect for placebo period.}
  \item{force}{user specified \code{force} option.}
  \item{T}{the number of time periods.}
  \item{N}{the total number of units.}
  \item{p}{the number of time-varying observables.}
  \item{r.cv}{the number of factors included in the model -- either
    supplied by users or automatically chosen via cross-validation.}
  \item{lambda.cv}{the optimal hyper-parameter in matrix completion method
  chosen via cross-validation.}
  \item{res}{residuals for observed values.}
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
  \item{mu}{estimated ground mean. }
  \item{xi}{estimated time fixed effects.}
  \item{alpha}{estimated unit fixed effects.}
  \item{validX}{a logic value indicating if multicollinearity exists.}
  \item{validF}{a logic value indicating if factor exists.}
  \item{time}{term for switch-on treatment effect.}
  \item{att}{switch-on treatment effect.}
  \item{group.att}{switch-on averge cohort effect.}
  \item{est.group.att}{inference for \code{group.att}.}
  \item{count}{count of each term for switch-on treatment effect.}
  \item{time.off}{term for switch-off treatment effect.}
  \item{att.off}{switch-off treatment effect.}
  \item{count.off}{count of each term for switch-off treatment effect.}
  \item{est.att.avg}{inference for \code{att.avg}.}
  \item{est.att.avg.unit}{inference for \code{att.avg.unit}.}
  \item{est.att}{inference for \code{att.on}.}
  \item{est.att.off}{inference for \code{att.off}.}
  \item{est.placebo}{inference for \code{att.placebo}.}
  \item{est.beta}{inference for \code{beta}.}
  \item{att.avg.boot}{bootstrap results for \code{att.avg}.}
  \item{att.avg.unit.boot}{bootstrap results for \code{att.avg.unit}.}
  \item{att.count.boot}{bootstrap results for \code{count}.}
  \item{att.off.boot}{bootstrap results for \code{att.avg.off}.}
  \item{att.off.count.boot}{bootstrap results for \code{count.off}.}
  \item{att.placebo.boot}{bootstrap results for \code{att.placebo}.}
  \item{att.bound}{equivalence confidence interval for equivalence test.}
  \item{att.off.bound}{equivalence confidence interval for equivalence test for 
  switch-off effect.}
  \item{beta.boot}{bootstrap results for \code{beta}.}
  \item{wald}{wald test results for pre-treatment fitting check.}
  \item{permute}{permutation test results for sharp null hypothesis.}
}
\author{
  Licheng Liu; Ye Wang; Yiqing Xu 
}
\references{
  Jushan Bai. 2009. "Panel Data Models with Interactive Fixed
  Effects." Econometrica 77:1229--1279.

  Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference
  with Interactive Fixed Effects Models." Political Analysis, Vol. 25, 
  Iss. 1, January 2017, pp. 57-76. Available at: \url{https://doi.org/10.1017/pan.2016.2}.

  Athey, Susan, et al. 2018 "Matrix completion methods for causal panel data models." arXiv preprint arXiv:1710.10251. Available 
  at: \url{https://https://arxiv.org/abs/1710.10251}. 
  
  Licheng Liu, et al. 2020. "A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional 
  Data." Working paper. Available at: \url{https://polmeth.mit.edu/sites/default/files/documents/Yiqing_Xu.pdf}. 


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
            CV = TRUE, r = c(0, 5), se = FALSE) 
print(out) 
}


