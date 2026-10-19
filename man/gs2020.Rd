\name{gs2020}
\alias{gs2020}
\docType{data}
\title{Minority Candidates and Coethnic Donations in U.S. House Elections}
\description{Replication data from Grumbach and Sahn (2020): U.S. House general elections, 1980 to 2012, at the district-cycle level (495 districts). Used in the user manual to study how the presence of an Asian candidate affects the share of campaign donations from Asian donors.}
\format{A data frame with 7,007 rows and 9 variables:
  \describe{
    \item{cycle}{election year (1980 to 2012, every two years).}
    \item{district_final}{congressional district, e.g. \code{"AK1"}.}
    \item{general_sharetotal_A_all}{share of donations from Asian donors.}
    \item{general_sharetotal_B_all}{share of donations from Black donors.}
    \item{general_sharetotal_H_all}{share of donations from Hispanic donors.}
    \item{cand_A_all}{1 if an Asian candidate ran in the general election, 0 otherwise.}
    \item{cand_B_all}{1 if a Black candidate ran, 0 otherwise.}
    \item{cand_H_all}{1 if a Hispanic candidate ran, 0 otherwise.}
    \item{dfe}{numeric district identifier.}
  }
}
\source{Grumbach, J. M., and Sahn, A. (2020). Race and Representation in Campaign Finance. \emph{American Political Science Review, 114}(1), 206-221.}
\keyword{datasets}
