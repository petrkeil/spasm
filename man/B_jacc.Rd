% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures_for_incidence_matrices.R
\name{B_jacc}
\alias{B_jacc}
\title{Pairwise Jaccard beta diversity among sites in a community matrix}
\usage{
B_jacc(m, clean = TRUE)
}
\arguments{
\item{m}{Community data matrix with species as rows and sites as columns. Can
contain either incidences (1/0) or abundances (natural numbers).}

\item{clean}{Logical. Do you want to remove zero-sum rows or columns? When the index needs
the "d" fraction (sites where no species occurs), only zero-sum rows are removed.}
}
\value{
A dist or data.frame objects with the between-sites pairwise beta values.
}
\description{
Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
with zero species are removed from m.
}
