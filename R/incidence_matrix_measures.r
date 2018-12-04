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
#' Classical 'C-score' metric of Stone -- seggregation among species in a community matrix
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
#' @references Stone L. & Roberts A. (1990) The checkerboard score and 
#' species distributions. Oecologia 85: 74-79.
#' @export

C_Stone <- function(m, lst = FALSE, fun = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]
  

  D <- vegan::designdist(m, 
                         method = "b*c",
                         abcd = TRUE, 
                         terms = "binary")

  if(lst) D <- dist2list(D)
  
  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)  
  
  return(D)
} 

# ------------------------------------------------------------------------------
#' 'C-score' metric of Stone & Roberts, scaled according to Ulrich & Gotelli 2012
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
#' @references Stone L. & Roberts A. (1990) The checkerboard score and 
#' species distributions. Oecologia 85: 74-79.
#' @references McNickle G.G. et al. (2018) Checkerboard score-area relationships
#' reveal spatial scales of plant community structure. Oikos 127: 415-426.
#' @export

C_Ulrich <- function(m, lst = FALSE, fun = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]
  
  
  D <- vegan::designdist(m, 
                         method = "b*c",
                         abcd = TRUE, 
                         terms = "binary")
  
  # standardization by number of sites (from Ulrich & Gotelli 2012)
  n <- ncol(m) 
  N <- (n*(n-1))/2
  
  D <- D/N
  
  if(lst) D <- dist2list(D)
  
  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)  
  
  return(D)
} 


# ------------------------------------------------------------------------------
#' Pairwise Jaccard similarities (proportional overlap) among species community matrix
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
#' Pairwise Simpson dissimilarities among species in a community matrix
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
#' Pairwise pearson correlations among species in community matrix
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

# ONLY USES SPECIES THAT OCCUPY LESS THAN ALL SITES
C_Pearson <- function(m, lst = FALSE, fun = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]
  
  m <- m[rowSums(m) < ncol(m),] # remove species that only have 1s
  
  D <- as.dist(cor(t(m), 
          method = "pearson"))
  if(lst) D <- dist2list(D)
  
  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)  
  
  return(D)
}

#C_Pearson(Atmar[[265]], lst=FALSE, fun=mean)


# ------------------------------------------------------------------------------

#' Conversion of distance matrix to list
#' @export

# function from package 'spaa' by Jinlong Zhang <jinlongzhang01@gmail.com>:
dist2list <- function (dist) 
{
  if (!class(dist) == "dist") {
    stop("the input data must be a dist object.")
  }
  dat <- as.data.frame(as.matrix(dist))
  if (is.null(names(dat))) {
    rownames(dat) <- paste(1:nrow(dat))
  }
  value <- stack(dat)$values
  rnames <- rownames(dat)
  namecol <- expand.grid(rnames, rnames)
  colnames(namecol) <- c("col", "row")
  res <- data.frame(namecol, value)
  return(res)
}

















