% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{ppp.trim}
\alias{ppp.trim}
\title{Function that removes, from point pattern, species with N lower than a specificed threshold}
\usage{
ppp.trim(comm.ppp, min.abu)
}
\arguments{
\item{comm.ppp}{A marked point pattern of spatstat's 'ppp', where marks identify each
point (individual) to a species.}

\item{min.abu}{The threshold. All species with N >= min.abu will be removed.}
}
\value{
ppp object
}
\description{
Function that removes, from point pattern, species with N lower than a specificed threshold
}
