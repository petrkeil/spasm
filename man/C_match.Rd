% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures_for_incidence_matrices.R
\name{C_match}
\alias{C_match}
\title{Simple matching coefficient of Sokal & Michener (1958)}
\usage{
C_match(m, clean = TRUE)
}
\arguments{
\item{m}{Community data matrix with species as rows and sites as columns. Can
contain either incidences (1/0) or abundances (natural numbers).}

\item{clean}{Logical. Do you want to remove zero-sum rows or columns? When the index needs
the "d" fraction (sites where no species occurs), only zero-sum rows are removed.}
}
\value{
A dist or data.frame objects with the pairwise association values.
}
\description{
Note: Before the calculation, species (rows) with zero occurrences are removed from m.
}
\references{
Sokal R.R. & Michener C.D. (1958) A statistical method for evaluating
systematic relationshps. The University of Kansas Scientific Bulletin 38: 1409-1438.

Hubalek Z. (1982) Coefficients of association and similarity, based
on binary (presence-absence) data: an evaluation. Biol. Rev. 57: 669-689.
}
