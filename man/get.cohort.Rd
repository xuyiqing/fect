\name{get.cohort}
\alias{get.cohort}
\title{Get the Cohort Index}
\description{Gets the cohort index given a panel data.}
\usage{get.cohort(data, D, index, 
                  varname = NULL, entry.time = NULL)
}
\arguments{
    \item{data}{a data frame, can be a balanced or unbalanced panel data.}
    \item{D}{the treatment indicator. The treatment should be binary (0 and 1).}
    \item{index}{a two-element string vector specifying the unit and time indicators. Must be of length 2. Every observation should be uniquely defined by the pair of the unit and time indicator.}
    \item{varname}{a string specifying the name for the generated cohort index.}
    \item{entry.time}{a list of intervals for first get-treated time.}    
}
\details{
  \code{get.cohort} pre-processes the data and generates the index for different cohorts..
}
\value{
  \item{data}{a new data frame containing the cohort index.}
}
\author{
  Licheng Liu; Ye Wang; Yiqing Xu, Ziyi Liu
}
\references{
  Jushan Bai. 2009. "Panel Data Models with Interactive Fixed Effects." Econometrica.
  
  Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models." Political 
  Analysis. 
  
  Athey, Susan, et al. 2021 "Matrix completion methods for causal panel data models." Journal of the American Statistical Association. 
  
  Licheng Liu, et al. 2022. "A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional 
  Data." American Journal of Political Science. 
  
  For more details about the matrix completion method, see \url{https://github.com/susanathey/MCPanel}. 
}
\seealso{
  \code{\link{fect}} and \code{\link{print.fect}}
}
\examples{
library(fect)
data(fect)
simdata.cohort <- get.cohort(data = simdata,D = 'D',index = c("id","time"))
}