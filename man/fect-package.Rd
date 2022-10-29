\name{fect-package}
\alias{fect-package}
\docType{package}
\title{Fixed Effects Counterfactual Estimators}
\description{
  The package implements counterfactual estimators in TSCS data analysis and statistical tools to test their identification assumptions.}
\details{
  It implements counterfactual estimators in TSCS data analysis. These estimators first impute counterfactuals for each treated 
  observation in a TSCS dataset by fitting an outcome model (fixed effects model, interactive fixed effects model, or matrix 
  completion) using the untreated observations. They then estimate the individualistic treatment effect for each treated observatio
  n by subtracting the predicted counterfactual outcome from its observed outcome. Finally, the average treatment effect on the 
  treated (ATT) or period-specific ATTs are calculated. A placebo test and an equivalence test are included to evaluate the 
  validity of identification assumptions behind these estimators.
  
  See \code{\link{fect}} for details.
}
\author{
  Licheng Liu <liulch@mit.edu>, MIT
  
  Ye Wang <yw1576@nyu.edu>, New York University
  
  Yiqing Xu <yiqingxu@stanford.edu >, Stanford University

  Ziyi Liu <zyliu2020@uchicago.edu>, University of Chicago
}
\references{
  Jushan Bai. 2009. "Panel Data Models with Interactive Fixed Effects." Econometrica.
  
  Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models." Political 
  Analysis. 
  
  Athey, Susan, et al. 2021 "Matrix completion methods for causal panel data models." Journal of the American Statistical Association. 
  
  Licheng Liu, et al. 2022. "A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional 
  Data." American Journal of Political Science. 
}
