\encoding{UTF-8}
\name{get.cohort}
\alias{get.cohort}
\title{Generate Cohort Information from Panel Data}
\description{
The \code{get.cohort} function processes panel data to generate cohort-related variables, including treatment timing, cohort membership, and treatment dynamics. This function ensures that the data is uniquely identified, handles staggered treatment adoption, and includes options for dropping always-treated units.
}
\usage{
get.cohort(data, 
           D, 
           index, 
           varname = NULL, 
           start0 = FALSE, 
           entry.time = NULL, 
           drop.always.treat = FALSE)
}
\arguments{
  \item{data}{A data frame containing the panel data.}
  \item{D}{A string specifying the name of the treatment indicator variable. The variable must be numeric and binary (0/1).}
  \item{index}{A character vector of length 2 specifying the unit and time variables (e.g., \code{c("id", "time")}). These variables uniquely identify each observation.}
  \item{varname}{A character vector of length 4 specifying the names of the output variables: \code{"FirstTreated"}, \code{"Cohort"}, \code{"Time_to_Treatment"}, and \code{"Time_to_Exit"}. If not provided, default names are used.}
  \item{start0}{A logical flag indicating whether the first post-treatment period should be indexed as 0. Default is \code{FALSE}.}
  \item{entry.time}{An optional list specifying custom entry time ranges for staggered adoption. Each list element should contain a numeric vector of length 2 indicating the start and end times. Default is \code{NULL}.}
  \item{drop.always.treat}{A logical flag indicating whether to remove units that are always treated for all periods. Default is \code{FALSE}.}
}
\details{
This function processes panel data to generate four key variables:
\itemize{
  \item \code{FirstTreated}: The first period when a unit is treated.
  \item \code{Cohort}: Cohort membership based on treatment timing.
  \item \code{Time_to_Treatment}: Time to treatment (negative for pre-treatment, positive for post-treatment).
  \item \code{Time_to_Exit}: Time to treatment cessation (if applicable).
}

The function ensures the data is uniquely identified by the specified \code{index} and handles both balanced and unbalanced panels. Custom cohort definitions can be provided via \code{entry.time}. Always-treated units can be dropped if desired.
}
\value{
A data frame with additional variables:
\item{FirstTreated}{The first period a unit is treated.}
\item{Cohort}{Cohort membership based on treatment timing.}
\item{Time_to_Treatment}{Time to treatment adoption.}
\item{Time_to_Exit}{Time to treatment cessation (if applicable).}
}
\examples{
library(fect)
data(fect)
simdata.cohort <- get.cohort(data = simdata, D = 'D', index = c("id","time"))
print(table(simdata.cohort[,'Cohort']))
}
\seealso{
\code{\link{fect}}, \code{\link{plot.fect}}
}
\author{
Tianzhu Qin; Yiqing Xu; Jinwen Wu
}
