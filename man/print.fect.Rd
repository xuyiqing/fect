\name{print.fect}
\alias{print.fect}
\title{Print Results}
\description{Print results of the matrix completion method.}
\usage{\method{print}{fect}(x, switch.on = TRUE,
       switch.off = FALSE,time.on.lim = NULL, time.off.lim = NULL, \dots)} 
\arguments{
  \item{x}{a \code{\link{fect}} object.}
  \item{switch.on}{a logical value that specifies whether to print switch.on effect.}
  \item{switch.off}{a logical value that specifies whether to print switch.off effect.}
  \item{time.on.lim}{a two-element numeric vector specifying the range of term of 
    switch-on treatment effects. e.g. \code{time.on.lim = c(-1,1)}.}
  \item{time.off.lim}{a two-element numeric vector specifying the range of term 
    of switch-off treatment effects. e.g. \code{time.off.lim = c(-1,1)}.}
  \item{\dots}{other argv.}
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

  Athey S, Bayati M, Doudchenko N, et al. Matrix completion methods for causal panel data models[J]. arXiv preprint arXiv:1710.10251, 2017.

  For more details about the matrix completion method, see \url{https://github.com/susanathey/MCPanel}. 
}
\seealso{
  \code{\link{fect}} and \code{\link{plot.fect}}
}


