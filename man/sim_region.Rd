\encoding{UTF-8}
\name{sim_region}
\docType{data}
\alias{sim_region}
\title{Simulated unbalanced panel with region-specific time effects}
\description{
A simulated unbalanced panel with 500 units belonging to 5 regions and 20 time periods.
The DGP includes region-specific linear time trends as the confounding source.
Treatment probability and timing depend on region, and units in higher-numbered
regions enter the panel later, creating unbalancedness.
Used to demonstrate additional fixed effects in the CFE estimator.
}
\format{A data frame with columns:
\describe{
  \item{time}{Time period (1 to 20)}
  \item{id}{Unit identifier (1 to 500)}
  \item{region}{Region indicator (1 to 5)}
  \item{Y}{Outcome variable}
  \item{D}{Treatment indicator (0/1)}
  \item{region_time}{Region-by-time interaction factor}
}
}
\keyword{datasets}
