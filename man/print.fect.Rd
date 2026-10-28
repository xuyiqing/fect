\encoding{UTF-8}
\name{print.fect}
\alias{print.fect}
\alias{summary.fect}
\alias{print.summary.fect}
\title{Print Results}
\description{Print the results of a \code{\link{fect}} fit: the call, the average treatment effect on the treated with its uncertainty, and the coefficients of covariates. For the object returned by \code{\link{effect}}, prints the cumulative effects.

For \code{method = "ife"} fitted by EM (\code{time.component.from = "notyettreated"}) and \code{method = "fe"}, a line under the estimator reports how the fit stopped: the number of rounds, the final relative change and \code{tol}, or that it hit \code{max.iteration} (see \code{convergence} in \code{\link{fect}}). When the bootstrap draws of the average effect break down (their SD is more than 5 times their IQR/1.349; see \code{boot.diag} in \code{\link{fect}}), the output ends with that check.

\code{summary()} prints the same, and always ends with the bootstrap check when one was computed (since 2.4.9).}
\usage{\method{print}{fect}(x, switch.on = TRUE,
       switch.off = FALSE,time.on.lim = NULL, time.off.lim = NULL, \dots)
\method{summary}{fect}(object, \dots)
\method{print}{summary.fect}(x, \dots)}
\arguments{
  \item{x}{a \code{\link{fect}} object (for \code{print.summary.fect}, the object \code{summary()} returns).}
  \item{object}{a \code{\link{fect}} object.}
  \item{switch.on}{not used in the current version (the period-by-period tables it controlled are no longer printed); see \code{x$est.att}.}
  \item{switch.off}{not used in the current version; see \code{x$est.att.off}.}
  \item{time.on.lim}{not used in the current version.}
  \item{time.off.lim}{not used in the current version.}
  \item{\dots}{other arguments.}
}
\value{\code{print()}: no return value. \code{summary()}: an object of class \code{"summary.fect"} (a list holding the fit as \code{object}), printed by its \code{print} method, which returns it invisibly.}
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
  \code{\link{fect}}, \code{\link{plot.fect}}
}
\examples{
library(fect)
data(simdata)
out <- fect(Y ~ D + X1 + X2, data = simdata, 
            index = c("id","time"), force = "two-way",
            CV = TRUE, r = c(0, 5), se = 0, parallel = FALSE) 
print(out)
summary(out)
}