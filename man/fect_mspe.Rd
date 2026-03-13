\name{fect_mspe}
\alias{fect_mspe}
\title{Mean Squared Prediction Error Evaluation for fect Objects}
\description{
Evaluates the prediction accuracy of one or more \code{fect} model fits by hiding a subset of control observations and comparing counterfactual predictions against actual values.
}
\usage{
fect_mspe(
  out.fect,
  hide_mask = NULL,
  hide_n = 20,
  seed = NULL,
  n_rep = 1,
  pre.trend = TRUE,
  pre.trend.n = 2
)
}
\arguments{
  \item{out.fect}{A fitted \code{fect} object or a named list of fitted \code{fect} objects to compare.}
  \item{hide_mask}{Optional matrix (TT x N or N x TT) or vector of length TT*N indicating which observations to hide. If \code{NULL}, positions are randomly selected.}
  \item{hide_n}{Integer; number of control observations to hide per replication. Default is 20.}
  \item{seed}{Optional integer; random seed for reproducibility.}
  \item{n_rep}{Integer; number of replications. Default is 1.}
  \item{pre.trend}{Logical; if \code{TRUE}, candidate hide positions are restricted to pre-treatment periods. Default is \code{TRUE}.}
  \item{pre.trend.n}{Integer; number of pre-treatment periods (counting backwards from treatment onset) from which to draw hide candidates. Default is 2.}
}
\value{
A list containing:
\item{summary}{A data frame with mean RMSE and Bias across replications for each model.}
\item{records}{A data frame with per-replication results including Rep, Model, Hidden_N, RMSE, and Bias.}
\item{hide_mask}{The hide mask used in the last replication.}
\item{fits}{The refitted \code{fect} objects from the last replication.}
}
\examples{
\dontrun{
out <- fect(Y ~ D, data = df, index = c("unit", "time"), method = "ife", r = 2)
mspe <- fect_mspe(out.fect = out, hide_n = 10, seed = 42, n_rep = 5)
mspe$summary
}
}
