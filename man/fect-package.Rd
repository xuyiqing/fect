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
}
\references{
  Jushan Bai. 2009. "Panel Data Models with Interactive Fixed Effects." Econometrica 77:1229--1279.
  
  Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models." Political 
  Analysis, Vol. 25, Iss. 1, January 2017, pp. 57-76. Available at: \url{https://doi.org/10.1017/pan.2016.2}. 
  
  Athey, Susan, et al. 2018 "Matrix completion methods for causal panel data models." arXiv preprint arXiv:1710.10251. Available at
  : \url{https://https://arxiv.org/abs/1710.10251}. 
  
  Licheng Liu, et al. 2020. "A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional 
  Data." Working paper. Available at: \url{https://polmeth.mit.edu/sites/default/files/documents/Yiqing_Xu.pdf}. 
}
