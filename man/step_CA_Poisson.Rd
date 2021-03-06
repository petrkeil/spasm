% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/null_model_simulations.R
\name{step_CA_Poisson}
\alias{step_CA_Poisson}
\title{One realization of the an algorithm that resamples each row according to a
Poisson distribution}
\usage{
step_CA_Poisson(m)
}
\arguments{
\item{m}{Community data matrix with species as rows and sites as columns. Can
contain either incidences (1/0) or abundances (natural numbers).}
}
\description{
First, for each species in the matrix, mean abunance is calculated.
Then N samples are taken for each species from a Poisson distribution, with lambda
equal to the mean abunance.
}
