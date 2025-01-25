\name{print.fect}
\alias{print.fect}
\title{Print Results}
\description{Print results of the matrix completion method.}
\usage{\method{print}{fect}(x, switch.on = TRUE,
       switch.off = FALSE,time.on.lim = NULL, time.off.lim = NULL, \dots)} 
\arguments{
  \item{x}{a \code{\link{fect}} object.}
<<<<<<< Updated upstream
  \item{switch.on}{a logical value that specifies whether to print switch.on effect.}
  \item{switch.off}{a logical value that specifies whether to print switch.off effect.}
  \item{time.on.lim}{a two-element numeric vector specifying the range of term of 
    switch-on treatment effects. e.g. \code{time.on.lim = c(-1,1)}.}
  \item{time.off.lim}{a two-element numeric vector specifying the range of term 
    of switch-off treatment effects. e.g. \code{time.off.lim = c(-1,1)}.}
  \item{\dots}{other argv.}
=======
  \item{switch.on}{logical; if \code{TRUE}, print switch-on effect.}
  \item{switch.off}{logical; if \code{TRUE}, print switch-off effect.}
  \item{time.on.lim}{two-element numeric vector specifying the switch-on effect range.}
  \item{time.off.lim}{two-element numeric vector specifying the switch-off effect range.}
  \item{\dots}{other arguments.}
>>>>>>> Stashed changes
}
\value{No return value.}
\author{
  Licheng Liu; Ye Wang; Yiqing Xu; Ziyi Liu
}
\references{  
<<<<<<< Updated upstream
  Jushan Bai. 2009. "Panel Data Models with Interactive Fixed Effects." Econometrica.
  
  Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models." Political 
  Analysis. 
  
  Athey, Susan, et al. 2021 "Matrix completion methods for causal panel data models." Journal of the American Statistical Association. 
  
  Licheng Liu, et al. 2022. "A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional 
  Data." American Journal of Political Science. 
  
  For more details about the matrix completion method, see \url{https://github.com/susanathey/MCPanel}. 
=======
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
>>>>>>> Stashed changes
}
\seealso{
  \code{\link{fect}}, \code{\link{plot.fect}}
}
\examples{
library(fect)
data(fect)
out <- fect(Y ~ D + X1 + X2, data = simdata1, 
            index = c("id","time"), force = "two-way",
            CV = TRUE, r = c(0, 5), se = 0, parallel = FALSE) 
print(out)
}


