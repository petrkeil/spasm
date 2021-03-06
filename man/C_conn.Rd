% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures_for_incidence_matrices.R
\name{C_conn}
\alias{C_conn}
\title{Network connectance}
\usage{
C_conn(m, clean = TRUE)
}
\arguments{
\item{m}{Community data matrix with species as rows and sites as columns. Can
contain either incidences (1/0) or abundances (natural numbers).}

\item{clean}{Logical. Do you want to remove zero-sum rows or columns? When the index needs
the "d" fraction (sites where no species occurs), only zero-sum rows are removed.}
}
\value{
A single number, network connectance.
}
\description{
Network connectance, defined as the proportion of all possible links in a network, adopted
from Carsten Dormann's bipartite package.
Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
with zero species are removed from m.
}
\references{
Dormann et al. (2009) Indices, graphs and null models: analyzing
bipartite ecological networks. The Open Ecology Journal, 2: 7-24.
}
