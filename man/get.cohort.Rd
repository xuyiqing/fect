\name{get.cohort}
\alias{get.cohort}
\title{Obtain the Cohort Index}
\description{Obtain the cohort index based on the provided panel data.}
\usage{
  get.cohort(data, D, index, 
                  varname = NULL, 
                  start0 = FALSE,
                  entry.time = NULL,
                  drop.always.treat = FALSE)
}
\arguments{
    \item{data}{a data frame, can be a balanced or unbalanced panel data.}
    \item{D}{the treatment indicator. The treatment should be binary (0 and 1).}
    \item{index}{a two-element string vector specifying the unit and time indicators. Must be of length 2. Every observation should be uniquely defined by the pair of the unit and time indicator.}
    \item{varname}{a vector of four names for the new columns, in this order: first treated period, cohort, time relative to treatment, and time relative to treatment exit. Default \code{NULL}: \code{"FirstTreat"}, \code{"Cohort"}, \code{"Time_to_Treatment"} and \code{"Time_to_Exit"}.}
    \item{start0}{a logical flag that indicates whether period 0 is the first post-treatment period. By default, it is set to FALSE, meaning that period 1 is considered the initial post-treatment period.}
    \item{entry.time}{a list of intervals representing the initial calendar time of treatment for each cohort. Units that received treatment within each specified interval are categorized as belonging to a single cohort.}  
    \item{drop.always.treat}{a logical flag that indicates whether to drop all always-treated units.}  
}
\details{
  \code{get.cohort} function preprocesses the data and generates the cohort index for different groups. If not specified in \code{varname}, the function automatically creates four new variables: 'FirstTreat' stores the initial calendar time when a unit is first treated, 'Cohort' stores the cohort time, 'Time_to_Treatment' tracks the period relative to the treatment, and 'Time_to_Exit' tracks the period relative to the treatment exit. Users have the option to specify custom names for these newly generated variables using the \code{varname} parameter. By default, the function sets \code{start0=FALSE}, where the relative period 1 corresponds to the first post-treatment/post-exit-treatment period. If \code{start0=TRUE}, period 0 becomes the first post-treatment/post-exit-treatment period. Additionally, users can categorize units into specific cohorts using the \code{entry.time} parameter. For instance, if \code{entry.time=list(c(1,5),c(6,10))}, units that receive their first treatment within the calendar intervals \code{c(1,5)} and \code{c(6,10)} are grouped into two distinct cohorts. Units that have never been treated are always categorized as the "control cohort".  
}
\value{
  The input data frame with the four new columns. Cohorts are labelled \code{"Cohort:<first treated period>"} (or by the intervals in \code{entry.time}), and never-treated units \code{"Control"}.
}
\author{
  Licheng Liu; Yiqing Xu, Ziyi Liu; Zhongyu Yin
}

