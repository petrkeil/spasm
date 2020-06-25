################################################################################
#'One realization of the IT algorithm of Ulrich & Gotelli (2010)
#'
#'@description Uses the IT algorithm by Ulrich & Gotelli (2010). This is the slower
#'but reliable R implementation.
#'
#'@param m Community data matrix with species as rows and sites as columns. Can
#'contain either incidences (1/0) or abundances (natural numbers).


step_CA_IT <- function(m)
{
  # row and column marginals that will be used as probabilities in the sampling
  rsums <- rowSums(m)
  csums <- colSums(m)
  N <- sum(m)
  # empty matrix in thich the randomized number will be located
  res <- matrix(0, nrow=nrow(m), ncol=ncol(m))

  # use categorical distribution to sample row and column coordinates
  rs <- sample(x=1:nrow(m), size=N, replace=TRUE, prob=rsums)
  cs <- sample(x=1:ncol(m), size=N, replace=TRUE, prob=csums)

  for(i in 1:N)
  {
    # for each sampled coordinates add 1 to the resulting matrix
    res[rs[i], cs[i]] <- res[rs[i], cs[i]] + 1
  }

  # A CORRECTION FOR RANDOMIZED MATRICES WITH 0 ROWS (= species)
  # -> WE SIMPLY ADD 1
      # which species have 0 individuals?
      zero.spec <- which(rowSums(res) == 0)
      # where will random 1 be placed?
      one.loc <- sample(1:ncol(m), size = sum(zero.spec), replace=TRUE)
      res[zero.spec, one.loc] <- 1

  return(res)
}

# TESTING
# x <- data.Ulrich[[1]]
# x <- matrix(c(1,10, 100, 0, 5, 0,0,1,0,0,0,0,
#              6,12, 499, 0, 1000,0,0,0,10,0,0,0), byrow=TRUE, nrow=2)

# par(mfrow=c(1,3))
# plot(colSums(x), colSums(step_CA_IT(x))); abline(a = 0, b=1)
# plot(rowSums(x), rowSums(step_CA_IT(x))); abline(a = 0, b=1)
# plot(x, step_CA_IT(x)); abline(a = 0, b=1)

# sum(x); sum(step_CA_IT(x))





################################################################################
#'One realization of the an algorithm that resamples each row according to a
#'Poisson distribution
#'
#'@description First, for each species in the matrix, mean abunance is calculated.
#'Then N samples are taken for each species from a Poisson distribution, with lambda
#'equal to the mean abunance.
#'
#'@param m Community data matrix with species as rows and sites as columns. Can
#'contain either incidences (1/0) or abundances (natural numbers).

step_CA_Poisson <- function(m)
{
  lambdas <- rowMeans(m)
  N = ncol(m)
  res <- t(sapply(X = lambdas, FUN = rpois, n = N))

  # A CORRECTION FOR RANDOMIZED MATRICES WITH 0 ROWS
    # which species have 0 individuals
    zero.spec <- which(rowSums(res) == 0)
    # where will random 1 placed?
    one.loc <- sample(1:ncol(m), size = sum(zero.spec), replace=TRUE)
    res[zero.spec, one.loc] <- 1
    return(res)
}

#m <- data.Ulrich[[1]]; rowSums(step_CA_Poisson(m))


################################################################################
#'One realization of the an algorithm that resamples each row randomly
#'
#'@description For each row, take the total number of individuals, and
#'re-assignn each of them with equal probability to each column.
#' This should actually be equivalent to the Poisson algorithm (i.e. the
#' distribution of abundances of each species should be Poisson).
#'
#'@param m Community data matrix with species as rows and sites as columns. Can
#'contain either incidences (1/0) or abundances (natural numbers).

step_CA_rowrandom <- function(m)
{
  for(i in 1:nrow(m))
  {
    print(i)
    x <- m[i,] # extract one row of m (species)
    new.x <- rep(0, times = length(x))
    smp <- table(sample(x = 1:length(x), size = sum(x), replace = TRUE)  )
    smp <- as.data.frame(smp)
    new.x[smp[,1]] <- smp[,2]
    m[i,] <- new.x # replace the old row with the resampled row
  }
  return(m)
}

# m <- data.Ulrich[[3]]; rowSums(step_CA_rowrandom(m))



################################################################################
#'One realization of the sim2 algorithm from EcoSimR
#'
#' @description https://cran.r-project.org/web/packages/EcoSimR/vignettes/CoOccurrenceVignette.html
#'
#' @param m Community data matrix with species as rows and sites as columns. Can
#' contain either incidences (1/0) or abundances (natural numbers).


step_C_sim2 <- function(m)
{
  t(apply(X = m, MARGIN = 1, FUN = sample))
}

# TESTING
# null_step_C_sim2(m)


################################################################################
#'Standardized Z-score for a given null model algorithm and ISA metric (for pairwise metrics)

#' @param m Community data matrix with species as rows and sites as columns. Can
#' contain either incidences (1/0) or abundances (natural numbers).
#' @param alrgorithm Character string. Name of the algorithm to be executed.
#' @param metric Character string. Name of spasm's pairwise association metric function.
#' @param N.sim Number of algirthm runs.
#' @param ... Additional arguments to the `metric` function.
#' @export
#'
Z_score <- function(m, algorithm, N.sim, metric, ...)
{
  # get the extra parameters to be passed to the metric function
  extra.pars <- as.list(match.call(Z_score, expand.dots = FALSE))$...
  #print(extra.pars)

  res <- array(dim=c(nrow(m), nrow(m), N.sim))

  for(i in 1:N.sim)
  {
    # radnomize using a given algorithm
    null.m <- do.call(algorithm, list(m))
    # calculate the metric
    null.metric <- as.matrix(do.call(metric, c(list(null.m), extra.pars)))
    res[,,i] <- null.metric
  }

  # calculate the Z-score
  obs <- do.call(metric, c(list(m), extra.pars) )
  mean.null <- as.dist(apply(X = res, MARGIN = c(1,2), FUN = mean))
  sd.null <- as.dist(apply(X = res, MARGIN = c(1,2), FUN = sd))

  # the Z-score
  Z <- (obs - mean.null) / sd.null
  return(Z)
}

# TESTING
# m <- round( data.Ulrich[[3]])
# m <- m[rowSums(m) > 0,]; m <- m[,colSums(m) > 0]
# x <- Z_score(m,
#             algorithm="step_CA_IT",
#             metric="CA_cov_cor", N.sim = 100,
#             correlation = FALSE,
#             method = "pearson")
# x




