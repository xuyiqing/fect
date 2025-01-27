\encoding{UTF-8}
\name{print.interFE}
\alias{print.interFE}
\title{Print Results}
\description{Print results of interactive fixed effects estimation.}
\usage{\method{print}{interFE}(x, \dots)}
\arguments{
  \item{x}{an \code{\link{interFE}} object.}
  \item{\dots}{other arguments.}
}
\author{
  Licheng Liu; Ye Wang; Yiqing Xu
}
\value{No return value.}
\references{
  Bai, J. (2009). 
  Panel data models with interactive fixed effects. 
  \emph{Econometrica, 77}(4), 1229-1279.
}
\seealso{
  \code{\link{interFE}}, \code{\link{fect}}
}
\examples{
library(fect)
data(fect)
d <- simdata[-(1:150),] # remove the treated units
out <- interFE(Y ~ X1 + X2, data = d, index=c("id","time"),
               r = 2, force = "two-way", nboots = 50)
print(out)
}