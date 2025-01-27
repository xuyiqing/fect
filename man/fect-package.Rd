\encoding{UTF-8}
\name{fect-package}
\alias{fect-package}
\docType{package}
\title{Fixed Effects Counterfactual Estimators}
\description{
  The \code{fect} package implements counterfactual estimators for Time-Series Cross-Sectional (TSCS) data analysis, along with statistical tools to test their identification assumptions.
}
\details{
  The \code{fect} package provides tools for estimating treatment effects in TSCS datasets using a range of counterfactual estimators. These estimators first impute counterfactuals for treated observations by fitting an outcome model (fixed effects model, interactive fixed effects model, or matrix completion) to untreated observations. Then, the individual treatment effect for each treated observation is calculated as the difference between the observed and predicted counterfactual outcomes.

  The package supports:
  \itemize{
    \item Calculation of the Average Treatment Effect on the Treated (ATT) and period-specific ATTs.
    \item Placebo tests and equivalence tests to evaluate the validity of the identification assumptions.
    \item Robust estimation techniques for unbalanced panel datasets.
  }

  See \code{\link{fect}} for details.
}
\author{
  Licheng Liu (\email{liulch@mit.edu}), Massachusetts Institute of Technology
  
  Ye Wang (\email{yw1576@nyu.edu}), New York University
  
  Yiqing Xu (\email{yiqingxu@stanford.edu}), Stanford University
  
  Ziyi Liu (\email{zyliu2020@uchicago.edu}), University of Chicago
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