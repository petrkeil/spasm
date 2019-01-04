



# ------------------------------------------------------------------------------
CA_Cov_or_Cor <- function(m, lst = FALSE, fun = FALSE, 
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

CA_Gower <- function(m, lst = FALSE, fun = FALSE)
{
  # eliminate empty rows and columns
  m <- m[rowSums(m) != 0, ]
  m <- m[, colSums(m) != 0]
  
  D <- vegan::vegdist(m, 
                      method = "gower")
  
  if(lst) D <- dist2list(D)
  
  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)  
  
  return(D)
}

# ------------------------------------------------------------------------------

CA_Bray <- function(m, lst = FALSE, fun = FALSE)
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

# CA_Bray(t(comm$abund.table))

# ------------------------------------------------------------------------------

CA_Ulrich <- function(m, lst = FALSE, fun = FALSE)
{
  # eliminate empty rows and columns
  n.spec <- nrow(m)
  combs <- combn(1:n.spec, m =2)
  
  D <- matrix(0, nrow = n.spec, ncol=n.spec)
  
  for(i in 1:ncol(combs))
  {
    sub.m <- m[combs[,i], ]
    A <- sub.m[1,]
    B <- sub.m[2,]
    # sites where abundances differ
    AneqB <- 1*((A == B) == FALSE)
    # sites where A has lower abundance
    Amin <- 1*(pmin(A,B) == A) * AneqB
    # sites where B has lower abundance
    Bmin <- 1*(pmin(A,B) == B) * AneqB
    # the distance
    D.i <- sum(Amin) * sum(Bmin)
    D[combs[1,i], combs[2,i]] <- D.i
  }
  D <- as.dist(t(D))
  
  n <- ncol(m) 
  N <- (n*(n-1))/2
  
  D <- D/N
  
  do.fun <- is.function(fun)
  if(do.fun) D <- fun(D)  
  
  return(D)
}

# CA_Ulrich(t(comm$abund.table))
