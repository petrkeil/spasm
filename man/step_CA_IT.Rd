% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/null_model_simulations.R
\name{step_CA_IT}
\alias{step_CA_IT}
\title{One realization of the IT algorithm of Ulrich & Gotelli (2010)}
\usage{
step_CA_IT(m)
}
\arguments{
\item{m}{Community data matrix with species as rows and sites as columns. Can
contain either incidences (1/0) or abundances (natural numbers).}
}
\description{
Uses the IT algorithm by Ulrich & Gotelli (2010). This is the slower
but more transparent (i.e. understandable by non-C++ users) R implementation.
}
