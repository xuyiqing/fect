\encoding{UTF-8}
\name{sim_linear}
\docType{data}
\alias{sim_linear}
\title{Simulated panel data with unit-specific linear time trends (block DID)}
\description{
A simulated balanced panel with 200 units and 50 time periods.
All 80 treated units receive treatment at period 41 (block DID design).
The DGP includes unit-specific slopes on a linear time trend,
where treated units have systematically steeper slopes.
Used to demonstrate unit-specific time trends with \code{Q.type = "linear"} in the CFE estimator.
}
\format{A data frame with 10000 rows and 4 columns:
\describe{
  \item{time}{Time period (1 to 50)}
  \item{id}{Unit identifier (1 to 200)}
  \item{Y}{Outcome variable}
  \item{D}{Treatment indicator (0/1)}
}
}
\keyword{datasets}
