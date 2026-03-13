\name{fect_mspe_sim}
\alias{fect_mspe_sim}
\title{Mean Squared Prediction Error Evaluation for fect Objects Using Simulated Counterfactuals}
\description{
Evaluates the prediction accuracy of one or more \code{fect} model fits by hiding observations and comparing counterfactual predictions against known \code{Y0} (true counterfactual) values from the source data.
}
\usage{
fect_mspe_sim(
  out.fect,
  hide_mask = NULL,
  hide_mask_y0 = NULL,
  hide_n = 20,
  seed = NULL,
  n_rep = 1
)
}
\arguments{
  \item{out.fect}{A fitted \code{fect} object or a named list of fitted \code{fect} objects to compare.}
  \item{hide_mask}{A 0/1 matrix (TT x N or N x TT) indicating which observations may be hidden. Required.}
  \item{hide_mask_y0}{Optional matrix of true counterfactual values (\code{Y0}) at hide positions. If \code{NULL}, values are extracted from the \code{Y0} column in the source data.}
  \item{hide_n}{Integer; number of observations to hide per replication. Default is 20.}
  \item{seed}{Optional integer; random seed for reproducibility.}
  \item{n_rep}{Integer; number of replications. Default is 1.}
}
\details{
This function requires a \code{Y0} column in the source data used to fit the original \code{fect} model. It is designed for simulation studies where the true counterfactual outcomes are known.
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
out <- fect(Y ~ D, data = simdata, index = c("unit", "time"), method = "ife", r = 2)
mask <- matrix(0, nrow = length(unique(simdata$time)), ncol = length(unique(simdata$unit)))
## set some positions to 1
mspe <- fect_mspe_sim(out.fect = out, hide_mask = mask, hide_n = 10, seed = 42)
mspe$summary
}
}
