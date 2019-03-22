#' Covariance or correlation matrix calculated from community matrix
#'
#' Calculates covariance or correlation matrix from a community abundance or
#' incidence matrix where species are rows and sites are columns.
#'
#' @param m Community data matrix with species as rows and sites as columns. Can
#' contain either incidences (1/0) or abundances (natural numbers).
#' @param lst Should the results be returned as a 'dist' object (FALSE), or
#' as a 'data.frame' (TRUE)?
#' @param fun Additional function (e.g. log-transformation) to be applied
#' to all pairwise values.
#' @param transf A transformation function to be applied to the community matrix.
#' The default is Hellinger transformation. Other options from 'vegan' function
#' 'decostand' are applicable.
#' @param correlation Should correlations be returned (TRUE) or covariances (FALSE).
#' @param method One of 'peason' (default), 'kendall', or 'spearman'.
#' @return A dist or data.frame objects with the pairwise association values.
#' @import vegan
#' @export

CA_cov_cor <- function(m,
                       lst = FALSE,
                       fun = FALSE,
                       transf = "hellinger",
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


  if(lst) D <- dist2list(D)

  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)

  return(D)
}


# ------------------------------------------------------------------------------
#' Hellinger distance matrix
#'
#' Calculates Hellinger distance matrix for mixed community matrices. In particular,
#' the abundance data are Hellinger-transformed, and then used for calcluation
#' of Euclidean distances.
#'
#' @inheritParams CA_cov_cor
#' @return A dist or data.frame objects with the pairwise values.
#' @references Legendre P. & De Caceres M. (2013) Beta diversity as the variance of community
#' data: dissimilarity coefficients and partitioning. Ecology Letters 16: 951-963.
#' @references Legendre P. & Legendre L. (2012) Numerical Ecology. Elsevier. pp 265-335.
#' @import vegan
#' @export

CA_hell <- function(m, lst = FALSE, fun = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  # Hellinger transformation of abundances
  m <- vegan::decostand(m, method = "hellinger")
  D <- dist(m) # Euclidean distances on the Hellinger-transformed matrix

  if(lst) D <- dist2list(D)

  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)

  return(D)
}

# ------------------------------------------------------------------------------
#' Chi-square distance matrix
#'
#' Calculates chi-square distance matrix for mixed community matrices. In particular,
#' the abundance data are chi-square-transformed, and then used for calcluation
#' of Euclidean distances.
#'
#' @inheritParams CA_cov_cor
#' @return A dist or data.frame objects with the pairwise values.
#' @references Legendre P. & De Caceres M. (2013) Beta diversity as the variance of community
#' data: dissimilarity coefficients and partitioning. Ecology Letters 16: 951-963.
#' @references Legendre P. & Legendre L. (2012) Numerical Ecology. Elsevier. pp 265-335.
#' @import vegan
#' @export

CA_chi <- function(m, lst = FALSE, fun = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  # Hellinger transformation of abundances
  m <- vegan::decostand(m, method = "chi.square")
  D <- dist(m) # Euclidean distances on the Hellinger-transformed matrix

  if(lst) D <- dist2list(D)

  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)

  return(D)
}


# ------------------------------------------------------------------------------
#' Bray-Curtis dissimilarity matrix
#'
#' Calculates Bray-Curtis dissimilarity matrix for a community matrix.
#' This is a wrapper to a specific call to the 'vegdist' function in 'vegan'.
#'
#' @inheritParams CA_cov_cor
#' @return A dist or data.frame objects with the pairwise values.
#' @references Legendre P. & De Caceres M. (2013) Beta diversity as the variance of community
#' data: dissimilarity coefficients and partitioning. Ecology Letters 16: 951-963.
#' @import vegan
#' @export

CA_bray <- function(m, lst = FALSE, fun = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  D <- vegan::vegdist(m,
                      method = "bray")

  if(lst) D <- dist2list(D)

  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)

  return(D)
}


# ------------------------------------------------------------------------------
#' Abundance-based Jaccard (a.k.a. Ruzicka) dissimilarity matrix
#'
#' Calculates abundance-based Jaccard dissimilarity matrix for a community matrix.
#' This is a wrapper to a specific call to the 'vegdist' function in 'vegan'.
#'
#' @inheritParams CA_cov_cor
#' @return A dist or data.frame objects with the pairwise values.
#' @references Legendre P. & De Caceres M. (2013) Beta diversity as the variance of community
#' data: dissimilarity coefficients and partitioning. Ecology Letters 16: 951-963.
#' @import vegan
#' @export

CA_ruz <- function(m, lst = FALSE, fun = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  D <- vegan::vegdist(m,
                      method = "bray")
  D <- 2*D/(1+D)

  if(lst) D <- dist2list(D)

  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)

  return(D)
}


