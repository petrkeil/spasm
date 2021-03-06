% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/null_model_simulations.R
\name{Z_score}
\alias{Z_score}
\title{Standardized Z-score for a given null model algorithm and ISA metric (for pairwise metrics)}
\usage{
Z_score(m, algorithm, N.sim, metric, ...)
}
\arguments{
\item{m}{Community data matrix with species as rows and sites as columns. Can
contain either incidences (1/0) or abundances (natural numbers).}

\item{N.sim}{Number of algorithm runs.}

\item{metric}{Character string. Name of spasm's pairwise association metric function.}

\item{...}{Additional arguments to the `metric` function.}

\item{alrgorithm}{Character string. Name of the algorithm to be executed.}
}
\description{
Standardized Z-score for a given null model algorithm and ISA metric (for pairwise metrics)
}
