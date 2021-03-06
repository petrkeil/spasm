% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures_for_incidence_matrices.R
\name{C_checker}
\alias{C_checker}
\title{Number of checkerboard species pairs}
\usage{
C_checker(m, clean = TRUE)
}
\arguments{
\item{m}{Community data matrix with species as rows and sites as columns. Can
contain either incidences (1/0) or abundances (natural numbers).}

\item{clean}{Logical. Do you want to remove zero-sum rows or columns? When the index needs
the "d" fraction (sites where no species occurs), only zero-sum rows are removed.}
}
\value{
A single number, the number of checkerboard species pairs.
}
\description{
This function follows the code from 'EcoSimR' package by Gotelli, Hard and Ellison,
specifically, it re-uses their 'checker' function.
Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
with zero species are removed from m.
}
