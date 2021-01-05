#' Forbes' coefficient of association among species in a community matrix
#'
#' Note: Before the calculation, species (rows) with only zero occurrences are removed from m.
#'
#' @param m Community data matrix with species as rows and sites as columns. Can
#' contain either incidences (1/0) or abundances (natural numbers).
#' @param clean Logical. Do you want to remove zero-sum rows or columns? When the index needs
#' the "d" fraction (sites where no species occurs), only zero-sum rows are removed.
#' @return A dist object with the pairwise association values.
#' @import vegan
#' @references Forbes S.A. (1907) On the local distribution of certain Illinois
#' fishes: An essay in statistical ecology. Bulletin of the Illinois State
#' Laboratory of Natural History, 7: 273-303.
#' @references Arita H. (2016) Species co-occurrence analysis: pairwise versus
#' matrix-level approaches. Global Ecology and Biogeography, 25: 1397-1400.
#' @export

C_forbes <- function(m, clean = TRUE)
{
  # eliminate empty rows, but not columns, since the index uses d
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
  }

  D <- vegan::designdist(m,
                         method = "a / (((a+b)*(a+c))/(a+b+c+d))",
                         abcd = TRUE,
                         terms = "binary")
  return(D)
}

# C_forbes(matrix(c(1,1,1,0,0,0,
#                0,0,0,1,1,1), byrow=T, nrow=2))
# C_forbes(matrix(c(1,1,1,1,0,0,
#                0,0,1,1,1,1), byrow=T, nrow=2))
# C_forbes(matrix(c(1,1,1,0,0,0,
#               1,1,1,0,0,0), byrow=T, nrow=2))


# ------------------------------------------------------------------------------
#' Alroy's (2015) modification of the Forbes index
#'
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @param m Community data matrix with species as rows and sites as columns. Can
#' contain either incidences (1/0) or abundances (natural numbers).
#' @return A dist object with the pairwise association values.
#' @import vegan
#' @references Alroy J. (2015) A new twist on a very old binary similarity coefficient.
#' Ecology 96: 575-586.
#' @export
#'

C_alroy <- function(m, clean = TRUE)
{
  # eliminate empty rows and columns
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
    m <- m[, colSums(m) != 0]
  }

  D <- vegan::designdist(m,
                         method = "a*( (a+b+c) + sqrt(a+b+c) )/( (a+b)*(a+c) + a*sqrt(a+b+c) + 0.5*(b*c) )",
                         abcd = TRUE,
                         terms = "binary")
  return(D)
}

#C_alroy(matrix(c(1,1,0,0,1,1,
#               0,0,1,0,0,1), byrow=T, nrow=2))

# ------------------------------------------------------------------------------
#' The classical 'C-score' (C_seg) seggregation metric of Stone & Roberts (1990),
#' and its scaled version
#'
#' The scaled version divides the raw C-score by the number of possible site
#' combinations, which is `n(n-1)/2` according to Ulrich & Gotelli (2013). Also
#' note that their scaling formula also includes `m`, the number of species. Thus,
#' Ulrich & Gotelli's (2013) formula for the scaled C-score is actually
#' `4*(C_seg)/(m*n*(m-1)(n-1))`. For two species (`m=2`) this reduces to  `C_seg / ((n*(n-1))/2)`,
#' which is the formula that I use here for the scaled version.
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @inheritParams C_forbes
#' @param scale Should the raw metric be scaled using the total number of possible
#' site combinations?
#' @return A dist object with the pairwise association values.
#' @import vegan
#' @references Stone L. & Roberts A. (1990) The checkerboard score and
#' species distributions. Oecologia 85: 74-79.
#' @references McNickle G.G. et al. (2018) Checkerboard score-area relationships
#' reveal spatial scales of plant community structure. Oikos 127: 415-426.
#' @references Ulrich W. & Gotelli N.J. (2013) Pattern detection in null model analysis.
#' Oikos 122: 2-18.
#' @export

C_seg <- function(m, scale = TRUE, clean = TRUE)
{
  # eliminate empty rows and columns
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
    m <- m[, colSums(m) != 0]
  }

  D <- vegan::designdist(m,
                         method = "b*c",
                         abcd = TRUE,
                         terms = "binary")
  if(scale)
  {
    # standardization by number of possible site combinations
    # (from Ulrich & Gotelli (2012) and  McNickle et al. (2018)
    n <- ncol(m) # number of sites
    N <- (n*(n-1))/2 # number of possible site combinations
    D <- D/N
  }

  if(length(D) == 0) D <- NA

  return(D)
}



# ------------------------------------------------------------------------------
#' 'Togetherness' (C_tog) metric of Stone & Roberts (1992), and its scaled version
#'
#' The scaled version divides the raw togetherness by the number of possible site
#' combinations, which is `n(n-1)/2` according to Ulrich & Gotelli (2013). Also
#' note that their scaling formula  includes `m`, the number of species. Thus,
#' Ulrich & Gotelli's (2013) formula for the scaled C-score is actually
#' `4*(C_tog)/(m*n*(m-1)(n-1))`. For two species (`m=2`) this reduces to  `C_tog / ((n*(n-1))/2)`,
#' which is the formula that I use here for the scaled version.
#' Note: Before the calculation, species (rows) with zero occurrences are removed from m.
#'
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

C_tog <- function(m, scale = TRUE, clean = TRUE)
{
  # eliminate empty rows (not columns, because we want the d fraction)
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
  }

  # convert to binary
  m[m>1] <- 1

  # according to Stone & Roberts 1992:
  # Tij = Sij(NI + Sij - ri - rj)
  # where
  # Sij = a
  # ri = a + b
  # rj = a + c
  # NI = a + b + c + d
  # which simplifes to a*d
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

 if(length(D) == 0) D <- NA

  return(D)
}

# ------------------------------------------------------------------------------
#' Pairwise Jaccard associations among species in a community matrix
#'
#' This is the positive association measure (as opposed to dissimilarity).
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @inheritParams C_forbes
#' @return A dist or data.frame objects with the pairwise association values.
#' @import vegan
#' @export

C_jacc <- function(m, clean = TRUE)
{
  # eliminate empty rows and columns
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
    m <- m[, colSums(m) != 0]
  }

  D <- vegan::designdist(m,
                         method = "a/(a+b+c)",
                         abcd = TRUE,
                         terms = "binary")
  return(D)
}


# ------------------------------------------------------------------------------
#' Pairwise Jaccard beta diversity among sites in a community matrix
#'
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @inheritParams C_forbes
#' @return A dist or data.frame objects with the between-sites pairwise beta values.
#' @import vegan
#' @export

B_jacc <- function(m, clean = TRUE)
{
  # eliminate empty rows and columns
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
    m <- m[, colSums(m) != 0]
  }

  D <- vegan::designdist(t(m), # the only difference between C_jacc is this transposition
                         method = "a/(a+b+c)",
                         abcd = TRUE,
                         terms = "binary")
  return(D)
}


# ------------------------------------------------------------------------------
#' Pearson tetrachoric correaltions for binary data
#'
#' This is the index A_30 from Hubalek (1982), who considered it to be one of
#' the recommended indices.
#' Note: Before the calculation, species (rows) with zero occurrences are removed from m.
#'
#' @inheritParams C_forbes
#' @return A dist or data.frame objects with the pairwise association values.
#' @references Hubalek Z. (1982) Coefficients of association and similarity, based
#' on binary (presence-absence) data: an evaluation. Biol. Rev. 57: 669-689.
#' @import vegan
#' @export

C_pears <- function(m, clean = TRUE)
{
  # eliminate empty rows (but not empty columns, because we want d)
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
  }

  D <- vegan::designdist(m,
                         method = "(a*d-b*c)/(((a+b)*(c+d)*(a+c)*(b+d))^0.5)",
                         abcd = TRUE,
                         terms = "binary")
  return(D)
}

# C_pears(matrix(c(1,1,1,0,0,0,
#                0,0,0,1,1,1), byrow=T, nrow=2))
#C_pears(matrix(c(1,1,1,1,0,0,
#                0,0,1,1,1,1), byrow=T, nrow=2))
#C_pears(matrix(c(1,1,1,0,0,0,
#                1,1,1,0,0,0), byrow=T, nrow=2))


# ------------------------------------------------------------------------------
#' Simple matching coefficient of Sokal & Michener (1958)
#'
#' Note: Before the calculation, species (rows) with zero occurrences are removed from m.
#'
#' @inheritParams C_forbes
#' @return A dist or data.frame objects with the pairwise association values.
#' @references Sokal R.R. & Michener C.D. (1958) A statistical method for evaluating
#' systematic relationshps. The University of Kansas Scientific Bulletin 38: 1409-1438.
#' @references Hubalek Z. (1982) Coefficients of association and similarity, based
#' on binary (presence-absence) data: an evaluation. Biol. Rev. 57: 669-689.
#' @import vegan
#' @export

C_match <- function(m, clean = TRUE)
{
  # eliminate empty rows (but not empty columns, since we want d)
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
  }

  D <- vegan::designdist(m,
                         method = "(a + d)/(a + b + c + d)",
                         abcd = TRUE,
                         terms = "binary")
  return(D)
}

# C_match(matrix(c(1,1,1,0,0,0,
#                0,0,0,1,1,1), byrow=T, nrow=2))
# C_match(matrix(c(1,1,1,1,0,0,
#                0,0,1,1,1,1), byrow=T, nrow=2))
# C_match(matrix(c(1,1,1,0,0,0,
#                1,1,1,0,0,0), byrow=T, nrow=2))



# ------------------------------------------------------------------------------
#' Pairwise Sorensen associations among species in a community matrix
#'
#' Arita (2017) also quotes this to be the "index of co-incidence" of Dice (1945).
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @inheritParams C_forbes
#' @return A dist or data.frame objects with the pairwise association values.
#' @references Arita H. (2015) Multisite and multispecies measures of overlap,
#' co-occurrence and co-diversity. Ecography 40: 709-718.
#' @references Dice L.R. (1945) Measures of the amount of ecological association
#' between species. Ecology 94: 2403-2414.
#' @import vegan
#' @export

C_sor <- function(m, clean = TRUE)
{
  # eliminate empty rows and columns
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
    m <- m[, colSums(m) != 0]
  }

  D <- vegan::designdist(m,
                         method = "(2*a)/(2*a+b+c)",
                         abcd = TRUE,
                         terms = "binary")
  return(D)
}

# C_sor(matrix(c(1,1,1,0,0,0,
#                0,0,0,1,1,1), byrow=T, nrow=2))
# C_sor(matrix(c(1,1,1,1,0,0,
#                0,0,1,1,1,1), byrow=T, nrow=2))
# C_sor(matrix(c(1,1,1,0,0,0,
#                1,1,1,0,0,0), byrow=T, nrow=2))


# ------------------------------------------------------------------------------
#' Pairwise Simpson association among species in a community matrix
#'
#' This returns the positive association measure (as opposed to dissimilarity).
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @inheritParams C_forbes
#' @return A dist or data.frame objects with the pairwise association values.
#' @import vegan
#' @export

C_sim <- function(m, clean = TRUE)
{
  # eliminate empty rows and columns
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
    m <- m[, colSums(m) != 0]
  }

  D <- vegan::designdist(m,
                    # method = "a/(pmin(b,c)+a)", # the similarity version
                    method = "pmin(b,c)/(pmin(b,c)+a)", # use this for dissimilarity
                    abcd = TRUE,
                    terms = "binary")
  D <- 1-D # convert to similarity (association) measure
  return(D)
}

# C_sim(matrix(c(1,1,1,0,0,0,
#                0,0,0,1,1,1), byrow=T, nrow=2))
# C_sim(matrix(c(1,1,1,1,0,0,
#                0,0,1,1,1,1), byrow=T, nrow=2))
# C_sim(matrix(c(1,1,1,0,0,0,
#                1,1,1,0,0,0), byrow=T, nrow=2))


# ------------------------------------------------------------------------------
#' Whittaker's index of overall association in a community matrix
#'
#' This metric is equivalent to calculating the classical Whittaker's beta diversity on
#' a transposed matrix. The impulse to use this metric for species associations arose
#' during a discussion with Aniko Toth at IBS 2019 conference in Malaga. We then
#' found out that Arita (2015) also used the same idea.
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @inheritParams C_forbes
#' @return A single number, which is the ratio of mean occupancy and the number of
#' @references Arita H. (2015) Multisite and multispecies measures of overlap,
#' co-occurrence and co-diversity. Ecography 40: 709-718.
#' @export

Whittaker <- function(m, clean = TRUE)
{
  # eliminate empty rows and columns
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
    m <- m[, colSums(m) != 0]
  }

  # convert to binary matrix
  m[m > 0] <- 1

  mean.occ <- mean(rowSums(m))
  N.sites <- ncol(m)

  return(N.sites/mean.occ)
}

# m <- data.Atmar[[1]]
# Whittaker(t(m))
# Whittaker(m)

# ------------------------------------------------------------------------------
#' Variance ratio of Schluter (1984)
#'
#' This is the classical variance ratio function; the function is applicable to
#' both incidence and abundance matrices.
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @inheritParams C_forbes
#' @return A single number, the variance ratio.
#'
#' @references Schluter, D. 1984. A variance test for detecting species associations,
#' with some example applications. Ecology 65: 998-1005.
#' @export

C_ratio <- function (m, clean = TRUE)
{
  # eliminate empty rows and columns
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
    m <- m[, colSums(m) != 0]
  }

  # convert to binary matrix
  m[m > 0] <- 1

  v <- var(colSums(m))/sum(apply(m, 1, var))
  return(v)
}

# m <- data.Atmar[[1]]
# C_ratio(m)

# ------------------------------------------------------------------------------
#' Network connectance
#'
#' Network connectance, defined as the proportion of all possible links in a network, adopted
#' from Carsten Dormann's bipartite package.
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @inheritParams C_forbes
#' @return A single number, network connectance.
#' @references Dormann et al. (2009) Indices, graphs and null models: analyzing
#' bipartite ecological networks. The Open Ecology Journal, 2: 7-24.
#' @export

C_conn <- function (m, clean = TRUE)
{
  # eliminate empty rows and columns
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
    m <- m[, colSums(m) != 0]
  }

  # convert to binary matrix
  m[m > 0] <- 1

  res <- sum(m>0) / ((nrow(m)*ncol(m)))
  return(res)
}

# m <- data.Atmar[[1]]
# C_conn(m)
# bipartite::networklevel(m, index="connectance")


# ------------------------------------------------------------------------------
#' Number of unique species combinations
#'
#' This function follows the code from EcoSimR package by Gotelli, Hard and Ellison,
#' specifically, it re-uses their 'species_combo' function.
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @inheritParams C_forbes
#' @return A single number, the number of unique species combinations.
#' @export

C_combo <- function (m, clean = TRUE)
{
  # eliminate empty rows and columns
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
    m <- m[, colSums(m) != 0]
  }

  # these two lines are adopted from the 'species_combo' function from EcoSimR:
  res <- ncol(unique(m, MARGIN = 2))
  return(res)
}

# m <- data.Atmar[[1]]
# C_combo(m)
# EcoSimR::species_combo(m)


# ------------------------------------------------------------------------------
#' Number of checkerboard species pairs
#'
#' This function follows the code from 'EcoSimR' package by Gotelli, Hard and Ellison,
#' specifically, it re-uses their 'checker' function.
#' Note: Before the calculation, species (rows) with zero occurrences, and sites (columns)
#' with zero species are removed from m.
#'
#' @inheritParams C_forbes
#' @return A single number, the number of checkerboard species pairs.
#' @export

C_checker <- function (m, clean = TRUE)
{
  # eliminate empty rows and columns
  if(clean)
  {
    m <- m[rowSums(m) != 0, ]
    m <- m[, colSums(m) != 0]
  }

  # the following code comes from the 'checker' function from EcoSimR:
  pairwise <- cbind(t(combn(nrow(m), 2)), 0)
  shared <- mat.or.vec(1, nrow(pairwise))
  for (i in 1:nrow(pairwise)) {
    shared[i] <- sum(m[pairwise[i, 1], ] == 1 & m[pairwise[i, 2], ] == 1)
  }
  res <- sum(shared == 0)
  return(res)
}

# m <- data.Atmar[[1]]
# C_checker(m)
# EcoSimR::checker(m)
















