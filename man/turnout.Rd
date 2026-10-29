\encoding{UTF-8}
\name{turnout}
\docType{data}
\alias{turnout}
\title{EDR and Voter Turnout in the US}
\description{State-level turnout in U.S. presidential elections, 1920 to 2012, with the adoption of Election-Day Registration (EDR), from Xu (2017), Section 5. Used in the Gsynth chapter of the user manual.

The panel is balanced, with 47 states (Alaska, Hawaii and North Dakota are not included) and 24 elections. Nine states adopt EDR: Maine, Minnesota and Wisconsin in 1976; Idaho, New Hampshire and Wyoming in 1996; Iowa and Montana in 2008; Connecticut in 2012. No state repeals any of the three registration policies in the data.}
\format{A data frame with 1,128 rows and 6 columns:
\describe{
  \item{abb}{state abbreviation (character).}
  \item{year}{election year (1920 to 2012, every four years).}
  \item{turnout}{voter turnout, in percent.}
  \item{policy_edr}{1 if the state has Election-Day Registration in that election, 0 otherwise (the treatment).}
  \item{policy_mail_in}{1 if the state has mail-in voter registration, 0 otherwise.}
  \item{policy_motor}{1 if the state has motor-voter registration, 0 otherwise.}
}
}
\keyword{datasets}
\references{
  Springer, M. J. (2014). 
  How the States Shaped the Nation: American Electoral Institutions and Voter Turnout, 1920-2000. 
  University of Chicago Press.
  
  Xu, Y. (2017). 
  Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models. 
  \emph{Political Analysis, 25}(1), 57-76.
}
