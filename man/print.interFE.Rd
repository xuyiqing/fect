\name{print.interFE}
\alias{print.interFE}
\title{Print Results}
\description{Print  results of interactive fixed effects estimation.}
\usage{\method{print}{interFE}(x, \dots)} 
\arguments{
  \item{x}{an \code{\link{interFE}} object.}
  \item{\dots}{other argv.}
}
\author{
  Licheng Liu; Ye Wang; Yiqing Xu
}
\value{No return value.}
\references{Jushan Bai. 2009. "Panel Data Models with Interactive Fixed
  Effects." Econometrica 77:1229--1279.}
\seealso{\code{\link{interFE}} and \code{\link{fect}}}
\examples{
library(fect)
data(fect)
d <- simdata1[-(1:150),] # remove the treated units
out <- interFE(Y ~ X1 + X2, data = d, index=c("id","time"),
               r = 2, force = "two-way", nboots = 50)
print(out)
}


