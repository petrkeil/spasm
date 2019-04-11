% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures_for_incidence_matrices.R
\name{C_pears}
\alias{C_pears}
\title{Pearson tetrachoric correaltions for binary data}
\usage{
C_pears(m, lst = FALSE, fun = FALSE)
}
\arguments{
\item{m}{Community data matrix with species as rows and sites as columns. Can
contain either incidences (1/0) or abundances (natural numbers).}

\item{lst}{Should the results be returned as a 'dist' object (FALSE), or
as a 'data.frame' (TRUE)?}

\item{fun}{Additional function (e.g. log-transformation) to be applied
to all pairwise values.}
}
\value{
A dist or data.frame objects with the pairwise association values.
}
\description{
Pearson tetrachoric correaltions for binary data
}
\references{
Hubalek Z. (1982) Coefficients of association and similarity, based
on binary (presence-absence) data: an evaluation. Biol. Rev. 57: 669-689.
}