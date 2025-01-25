\name{interFE}
\alias{interFE}
\title{Interactive Fixed Effects Models}
\description{Estimating interactive fixed effect models.}
\usage{interFE(formula = NULL, data, Y, X = NULL, W = NULL,
         index, r = 0, force = "two-way",
         se = FALSE, nboots = 500, seed = NULL, 
         tol = 1e-3, max_iteration = 500, 
         binary = FALSE, QR = FALSE, normalize = FALSE)
}
\arguments{
  \item{formula}{an object of class "formula": a symbolic description of the model to be fitted.}
  \item{data}{a data frame.}
  \item{Y}{outcome variable.}
  \item{X}{time-varying covariates.}
  \item{W}{weights.}
  \item{index}{a two-element string vector specifying the unit (group) and time indicators.}
  \item{r}{an integer specifying the number of factors.}
<<<<<<< Updated upstream
  \item{force}{a string indicating whether unit or time fixed effects will be
    imposed. Must be one of the following,
    "none", "unit", "time", or "two-way". The default is "two-way".}
  \item{se}{a logical flag indicating whether uncertainty estimates will
    be produced via bootstrapping.}
  \item{nboots}{an integer specifying the number of bootstrap
    runs. Ignored if \code{se = FALSE}.}
  \item{seed}{an integer that sets the seed in random number
    generation. Ignored if  \code{se = FALSE} and \code{r} is specified.}
  \item{tol}{a numeric value that specifies tolerate level.}
  \item{max_iteration}{the maximal number of iterations for the EM algorithm.}
  \item{binary}{a logical flag indicating whether a probit link function will be 
    used.}
  \item{QR}{a logical flag indicating whether QR decomposition will be used for 
    factor analysis in probit model.}
  \item{normalize}{a logic flag indicating whether to scale outcome and 
    covariates. Useful for accelerating computing speed when magnitude of data is large.The default is \code{normalize=FALSE}.}
=======
  \item{force}{a string indicating whether unit/time/both fixed effects will be imposed: "none", "unit", "time", or "two-way".}
  \item{se}{logical; if \code{TRUE}, computes bootstrap-based standard errors.}
  \item{nboots}{number of bootstrap runs (ignored if \code{se=FALSE}).}
  \item{seed}{random seed.}
  \item{tol}{tolerance for EM algorithm.}
  \item{max_iteration}{max number of EM iterations.}
  \item{binary}{logical flag for a probit link (not fully supported here).}
  \item{QR}{logical flag for QR-based factor analysis in probit model (not fully supported).}
  \item{normalize}{logical; if \code{TRUE}, scales outcome/covariates.}
>>>>>>> Stashed changes
}
\details{
  \code{interFE} estimates interactive fixed effect models as in Bai (2009).
}
\value{
  \item{beta}{estimated coefficients.}
  \item{mu}{estimated grand mean.}
  \item{factor}{estimated time-varying factors.}
  \item{lambda}{estimated factor loadings.}
<<<<<<< Updated upstream
  \item{VNT}{a diagonal matrix that consists of the r eigenvalues.}
  \item{niter}{the number of iteration before convergence.}
  \item{alpha}{estimated unit fixed effect (if \code{force} is "unit" or "two-way").}
  \item{xi}{estimated time fixed effect (if \code{force} is "time" or "two-way").}
  \item{residuals}{residuals of the estimated interactive fixed effect model.}
  \item{sigma2}{mean squared error of the residuals.}
  \item{IC}{the information criterion.}
  \item{ValidX}{a logical flag specifying whether there are valid covariates.}
  \item{dat.Y}{a matrix storing data of the outcome variable.}
  \item{dat.X}{an array storing data of the independent variables.}
=======
  \item{VNT}{diagonal matrix of r eigenvalues.}
  \item{niter}{number of iterations before convergence.}
  \item{alpha}{estimated unit fixed effects (if imposed).}
  \item{xi}{estimated time fixed effects (if imposed).}
  \item{residuals}{model residuals.}
  \item{sigma2}{residual mean squared error.}
  \item{IC}{information criterion.}
  \item{ValidX}{logical for whether valid covariates exist.}
  \item{dat.Y}{matrix of outcome data.}
  \item{dat.X}{array of independent variables.}
>>>>>>> Stashed changes
  \item{Y}{name of the outcome variable.}
  \item{X}{name of time-varying control variables.}
  \item{index}{name of unit/time indicators.}
  \item{est.table}{table of final estimates.}
  \item{est.boot}{matrix of bootstrap results.}
}
\author{
  Licheng Liu; Ye Wang; Yiqing Xu
}
\references{Jushan Bai. 2009. "Panel Data Models with Interactive Fixed
  Effects." Econometrica 77:1229--1279.}
\seealso{
  \code{\link{print.interFE}}, \code{\link{fect}}
}
\examples{
library(fect)
data(fect)
d <- simdata1[-(1:150),] # remove the treated units
out <- interFE(Y ~ X1 + X2, data = d, index=c("id","time"),
               r = 2, force = "two-way", nboots = 50)
}


