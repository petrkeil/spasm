#' Forbes' coefficient of association among species in a community matrix
#' 
#' @param m Community data matrix with species as rows and sites as columns. Can 
#' contain either incidences (1/0) or abundances (natural numbers).
#' @param lst Should the results be returned as a 'dist' object (FALSE), or
#' as a 'data.frame' (TRUE)?
#' @param fun Additional function (e.g. log-transformation) to be applied 
#' to all pairwise values.
#' 
#' @return A list or data.frame objects with the pairwise association values.
#' @import vegan
#' @references Forbes S.A. (1907) On the local distribution of certain Illinois
#' fishes: An essay in statistical ecology. Bulletin of the Illinois State 
#' Laboratory of Natural History, 7: 273-303.
#' @references Arita H. (2016) Species co-occurrence analysis: pairwise versus 
#' matrix-level approaches. Global Ecology and Biogeography, 25: 1397-1400.
#' @export
#' 

C_Forbes <- function(m, lst = FALSE, fun = FALSE)
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
#' 'C-score' metric of Stone & Roberts (1990), and its scaled version
#' 
#' @param m Community data matrix with species as rows and sites as columns. Can 
#' contain either incidences (1/0) or abundances (natural numbers).
#' @param lst Should the results be returned as a 'dist' object (FALSE), or
#' as a 'data.frame' (TRUE)?
#' @param fun Additional function (e.g. log-transformation) to be applied 
#' to all pairwise values.
#' @param scale Should the raw metric be scaled using the total number of possible
#' site combinations?
#' 
#' @return A list or data.frame objects with the pairwise association values.
#' @import vegan
#' @references Stone L. & Roberts A. (1990) The checkerboard score and 
#' species distributions. Oecologia 85: 74-79.
#' @references McNickle G.G. et al. (2018) Checkerboard score-area relationships
#' reveal spatial scales of plant community structure. Oikos 127: 415-426.
#' @export

C_Seg <- function(m, lst = FALSE, fun = FALSE, scale = FALSE)
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
#' @param m Community data matrix with species as rows and sites as columns. Can 
#' contain either incidences (1/0) or abundances (natural numbers).
#' @param lst Should the results be returned as a 'dist' object (FALSE), or
#' as a 'data.frame' (TRUE)?
#' @param fun Additional function (e.g. log-transformation) to be applied 
#' to all pairwise values.
#' @param scale Should the raw metric be scaled using the total number of possible
#' site combinations?
#' 
#' @return A list or data.frame objects with the pairwise association values.
#' @import vegan
#' @references Stone L. & Roberts A. (1992) Competitive exclusion, or species
#' aggregation? An aid in deciding. Oecologia 91: 419-424.
#' @references Ulrich W. & Gogelli N.J. (2013) Pattern detection in null model
#' analysis. Oikos 122: 2-18.
#' @export

C_Tog <- function(m, lst = FALSE, fun = FALSE, scale = FALSE)
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
#' @param m Community data matrix with species as rows and sites as columns. Can 
#' contain either incidences (1/0) or abundances (natural numbers).
#' @param lst Should the results be returned as a 'dist' object (FALSE), or
#' as a 'data.frame' (TRUE)?
#' @param fun Additional function (e.g. log-transformation) to be applied 
#' to all pairwise values.
#' 
#' @return A list or data.frame objects with the pairwise association values.
#' @import vegan
#' @export

C_Jacc <- function(m, lst = FALSE, fun = FALSE)
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
#' Pairwise Simpson seggregation among species in a community matrix
#' 
#' @param m Community data matrix with species as rows and sites as columns. Can 
#' contain either incidences (1/0) or abundances (natural numbers).
#' @param lst Should the results be returned as a 'dist' object (FALSE), or
#' as a 'data.frame' (TRUE)?
#' @param fun Additional function (e.g. log-transformation) to be applied 
#' to all pairwise values.
#' 
#' @return A list or data.frame objects with the pairwise association values.
#' @import vegan
#' @export

C_Sim <- function(m, lst = FALSE, fun = FALSE)
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
#' Variance ratio of Schluter 
#' 
#' @param m Community data matrix with species as rows and sites as columns. Can 
#' contain either incidences (1/0) or abundances (natural numbers).
#' 
#' @return A single number, the variance ratio.
#' 
#' @references Schluter, D. 1984. A variance test for detecting species associations, with some example applications. Ecology 65: 998-1005.
#' 
#' @import vegan
#' @export


V_ratio <- function (m) 
{
  m <- m[which(rowSums(m) > 0), ]
  v <- var(colSums(m))/sum(apply(m, 1, var))
  return(v)
}

























