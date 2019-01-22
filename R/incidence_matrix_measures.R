#' Forbes' coefficient of association among species in a community matrix
#'
#' @param m Community data matrix with species as rows and sites as columns. Can
#' contain either incidences (1/0) or abundances (natural numbers).
#' @param lst Should the results be returned as a 'dist' object (FALSE), or
#' as a 'data.frame' (TRUE)?
#' @param fun Additional function (e.g. log-transformation) to be applied
#' to all pairwise values.
#'
#' @return A dist or data.frame objects with the pairwise association values.
#' @import vegan
#' @references Forbes S.A. (1907) On the local distribution of certain Illinois
#' fishes: An essay in statistical ecology. Bulletin of the Illinois State
#' Laboratory of Natural History, 7: 273-303.
#' @references Arita H. (2016) Species co-occurrence analysis: pairwise versus
#' matrix-level approaches. Global Ecology and Biogeography, 25: 1397-1400.
#' @export
#'

C_forbes <- function(m, lst = FALSE, fun = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  D <- vegan::designdist(m,
                         method = "a / (((a+b)*(a+c))/(a+b+c+d))",
                         abcd = TRUE,
                         terms = "binary")

  if(lst) D <- dist2list(D)

  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)

  return(D)
}



# ------------------------------------------------------------------------------
#' The classical 'C-score' seggregation metric of Stone & Roberts (1990),
#' and its scaled version
#'
#' @inheritParams C_forbes
#' @param scale Should the raw metric be scaled using the total number of possible
#' site combinations?
#' @return A dist or data.frame objects with the pairwise association values.
#' @import vegan
#' @references Stone L. & Roberts A. (1990) The checkerboard score and
#' species distributions. Oecologia 85: 74-79.
#' @references McNickle G.G. et al. (2018) Checkerboard score-area relationships
#' reveal spatial scales of plant community structure. Oikos 127: 415-426.
#' @export

C_seg <- function(m, lst = FALSE, fun = FALSE, scale = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  D <- vegan::designdist(m,
                         method = "b*c",
                         abcd = TRUE,
                         terms = "binary")
  if(scale)
  {
    # standardization by number of possible site combinations
    # (from Ulrich & Gotelli (2012) and  McNickle et al. (2018)
    n <- ncol(m)
    N <- (n*(n-1))/2
    D <- D/N
  }

  if(lst) D <- dist2list(D)

  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)

  return(D)
}



# ------------------------------------------------------------------------------
#' 'Togetherness' metric of Stone & Roberts (1992), and its scaled version
#'
#' @inheritParams C_forbes
#' @inheritParams C_seg
#' @param scale Should the raw metric be scaled using the total number of possible
#' site combinations?
#'
#' @return A dist or data.frame objects with the pairwise association values.
#' @import vegan
#' @references Stone L. & Roberts A. (1992) Competitive exclusion, or species
#' aggregation? An aid in deciding. Oecologia 91: 419-424.
#' @references Ulrich W. & Gogelli N.J. (2013) Pattern detection in null model
#' analysis. Oikos 122: 2-18.
#' @export

C_tog <- function(m, lst = FALSE, fun = FALSE, scale = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]


  D <- vegan::designdist(m,
                         method = "a*d",
                         abcd = TRUE,
                         terms = "binary")

  if(scale)
  {
    # standardization by number of possible site combinations
    # (from Ulrich & Gotelli (2012) and  McNickle et al. (2018)
    n <- ncol(m)
    N <- (n*(n-1))/2
    D <- D/N
  }

  if(lst) D <- dist2list(D)

  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)

  return(D)
}



# ------------------------------------------------------------------------------
#' Pairwise Jaccard associations among species in a community matrix
#'
#' @inheritParams C_forbes
#' @return A dist or data.frame objects with the pairwise association values.
#' @import vegan
#' @export

C_jacc <- function(m, lst = FALSE, fun = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  D <- vegan::designdist(m,
                         method = "a/(a+b+c)",
                         abcd = TRUE,
                         terms = "binary")

  if(lst) D <- dist2list(D)

  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)

  return(D)
}

# ------------------------------------------------------------------------------
#' Pairwise Sorensen associations among species in a community matrix
#'
#' @inheritParams C_forbes
#' @return A dist or data.frame objects with the pairwise association values.
#' @import vegan
#' @export

C_sor <- function(m, lst = FALSE, fun = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  D <- vegan::designdist(m,
                         method = "(2*a)/(2*a+b+c)",
                         abcd = TRUE,
                         terms = "binary")

  if(lst) D <- dist2list(D)

  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)

  return(D)
}


# ------------------------------------------------------------------------------
#' Pairwise Simpson seggregation among species in a community matrix
#'
#' @inheritParams C_forbes
#' @return A dist or data.frame objects with the pairwise association values.
#' @import vegan
#' @export

C_sim <- function(m, lst = FALSE, fun = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  D <- vegan::designdist(m,
                    method = "pmin(b,c)/(pmin(b,c)+a)",
                    abcd = TRUE,
                    terms = "binary")

  if(lst) D <- dist2list(D)

  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)

  return(D)
}


# ------------------------------------------------------------------------------
#' Association metric of Dale (1999)
#'
#' This is a metric described on page 147 of Dale (1999).
#'
#' @inheritParams C_forbes
#' @return A dist or data.frame objects with the pairwise association values.
#' @import vegan
#' @references Dale M.T.D (1999) Spatial pattern analysis in plant ecology.
#' Cambridge University Press.
#' @export

C_dale <- function(m, lst = FALSE, fun = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]

  D <- vegan::designdist(m,
                         method = "a*d - b*c",
                         abcd = TRUE,
                         terms = "binary")

  if(lst) D <- dist2list(D)

  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)

  return(D)
}

# ------------------------------------------------------------------------------
#' Toth's association
#'
#' This metric is equivalent to calculating the classical Whittaker's beta diversity on
#' a transposed matrix. The idea to use this metric for species associations arose
#' during a discussion with Aniko Toth at IBS 2019 conference in Malaga, and hence
#' the name of the metric.
#'
#' @inheritParams C_forbes
#' @return A single number, which is the ratio of mean occupancy and the number of
#' @export

C_toth <- function(m)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]
  # convert to binary matrix
  m[m > 0] <- 1

  mean.occ <- mean(rowSums(m))
  N.sites <- ncol(m)

  return(N.sites/mean.occ)
}


# ------------------------------------------------------------------------------
#' Variance ratio of Schluter (1984)
#'
#' @inheritParams C_forbes
#' @return A single number, the variance ratio.
#'
#' @references Schluter, D. 1984. A variance test for detecting species associations,
#' with some example applications. Ecology 65: 998-1005.
#' @export

V_ratio <- function (m)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]
  # convert to binary matrix
  m[m > 0] <- 1

  v <- var(colSums(m))/sum(apply(m, 1, var))
  return(v)
}

# ------------------------------------------------------------------------------
#' Network connectance
#'
#' Network connectance, defined as the proportion of all possible links in a network, adopted
#' from Carsten Dormann's bipartite package.
#'
#' @inheritParams C_forbes
#' @return A single number, network connectance.
#' @references Dormann et al. (2009) Indices, graphs and null models: analyzing
#' bipartite ecological networks. The Open Ecology Journal, 2: 7-24.
#' @export

N_connect <- function (m)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]
  # convert to binary matrix
  m[m > 0] <- 1

  res <- sum(m>0)/(nrow(m)*ncol(m))
  return(res)
}


# ------------------------------------------------------------------------------
#' Number of unique species combinations
#'
#' This function follows the code from EcoSimR package by Gotelli, Hard and Ellison,
#' specifically, it re-uses their 'species_combo' function.
#'
#' @inheritParams C_forbes
#' @return A single number, the number of unique species combinations.
#' @export

N_combo <- function (m)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]
  # convert to binary matrix
  m[m > 0] <- 1

  # these two lines are adopted from the 'species_combo' function from EcoSimR:
  res <- ncol(unique(m, MARGIN = 2))
  return(res)
}

# ------------------------------------------------------------------------------
#' Number of checkerboard species pairs
#'
#' This function follows the code from EcoSimR package by Gotelli, Hard and Ellison,
#' specifically, it re-uses their 'checker' function.
#'
#' @inheritParams C_forbes
#' @return A single number, the number of checkerboard species pairs.
#' @export

N_checker <- function (m)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]
  # convert to binary matrix
  m[m > 0] <- 1

  # the following code comes from the 'checker' function from EcoSimR:
  pairwise <- cbind(t(combn(nrow(m), 2)), 0)
  shared <- mat.or.vec(1, nrow(pairwise))
  for (i in 1:nrow(pairwise)) {
    shared[i] <- sum(m[pairwise[i, 1], ] == 1 & m[pairwise[i,
                                                           2], ] == 1)
  }
  res <- sum(shared == 0)
  return(res)
}

















