#' Covariance or correlation matrix calculated from community matrix
#'
#' Calculates covariance or correlation matrix from a community abundance or
#' incidence matrix where species are rows and sites are columns.
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @param m Community data matrix with species as rows and sites as columns. Can
#' contain either incidences (1/0) or abundances (natural numbers).
#' @param transf A transformation function to be applied to the community matrix.
#' The default is no transformation. Other options from 'vegan' function
#' 'decostand' are applicable.
#' @param correlation Should correlations be returned (TRUE) or covariances (FALSE).
#' @param method One of 'peason' (default), 'kendall', or 'spearman'.
#' is the default.
#' @return A dist object with the pairwise association values.
#' @import vegan
#' @export

CA_cov_cor <- function(m,
                       transf = NULL,
                       correlation = TRUE,
                       method = "pearson")
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  if(is.null(transf) == FALSE)
  {
    m <- vegan::decostand(m, method = transf)
  }

  if(correlation)
  {
    D <- as.dist(cor(t(m), method = method))
  }

  else
  {
    D <- as.dist(cov(t(m), method = method))
  }

  return(D)
}

# ------------------------------------------------------------------------------
#' Hellinger distance matrix
#'
#' Calculates Hellinger distances for all pairs of species in a community matrix. In particular,
#' the abundance data are Hellinger-transformed, and then used for calcluation
#' of Euclidean distances.
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @inheritParams CA_cov_cor
#' @return A dist object with the pairwise values.
#' @references Legendre P. & De Caceres M. (2013) Beta diversity as the variance of community
#' data: dissimilarity coefficients and partitioning. Ecology Letters 16: 951-963.
#' @references Legendre P. & Legendre L. (2012) Numerical Ecology. Elsevier. pp 265-335.
#' @import vegan
#' @export

CA_hell <- function(m)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  # Hellinger transformation of abundances
  m <- vegan::decostand(m, method = "hellinger")
  D <- dist(m) # Euclidean distances on the Hellinger-transformed matrix

  return(D)
}

# m <- data.Ulrich[[1]]
# CA_hell(m)

# ------------------------------------------------------------------------------
#' Chi-square distance matrix
#'
#' Calculates chi-square distances for all pairs of species in a community matrix.
#' In particular,
#' the abundance data are chi-square-transformed, and then used for calcluation
#' of Euclidean distances.
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @inheritParams CA_cov_cor
#' @return A dist object with the pairwise values.
#' @references Legendre P. & De Caceres M. (2013) Beta diversity as the variance of community
#' data: dissimilarity coefficients and partitioning. Ecology Letters 16: 951-963.
#' @references Legendre P. & Legendre L. (2012) Numerical Ecology. Elsevier. pp 265-335.
#' @import vegan
#' @export

CA_chi <- function(m)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  # Hellinger transformation of abundances
  m <- vegan::decostand(m, method = "chi.square")
  D <- dist(m) # Euclidean distances on the Chi.shquare-transformed matrix

  return(D)
}

# m <- data.Ulrich[[1]]
# CA_chi(m)

# ------------------------------------------------------------------------------
#' Bray-Curtis dissimilarity matrix
#'
#' Calculates all Bray-Curtis dissimilarities for all pairs of species in a community matrix.
#' This is a wrapper to a specific call to the 'vegdist' function in 'vegan'.
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @inheritParams CA_cov_cor
#' @return A dist object with the pairwise values.
#' @references Legendre P. & De Caceres M. (2013) Beta diversity as the variance of community
#' data: dissimilarity coefficients and partitioning. Ecology Letters 16: 951-963.
#' @import vegan
#' @export

CA_bray <- function(m)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  D <- vegan::vegdist(m,
                      method = "bray")

  return(D)
}

# m <- data.Ulrich[[1]]
# CA_bray(m)

# ------------------------------------------------------------------------------
#' Abundance-based Jaccard (a.k.a. Ruzicka) dissimilarity matrix
#'
#' Calculates abundance-based Jaccard dissimilarities for all species pairs
#' in a community matrix.
#' This is a wrapper to a specific call to the 'vegdist' function in 'vegan'.
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @inheritParams CA_cov_cor
#' @return A dist object with the pairwise values.
#' @references Legendre P. & De Caceres M. (2013) Beta diversity as the variance of community
#' data: dissimilarity coefficients and partitioning. Ecology Letters 16: 951-963.
#' @import vegan
#' @export

CA_ruz <- function(m)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  D <- vegan::vegdist(m,
                      method = "bray")
  D <- 2*D/(1+D)

  return(D)
}


