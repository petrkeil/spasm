################################################################################
#'One realization of the IT algorithm of Ulrich & Gotelli (2010)
#'
#'@description Uses the IT algorithm by Ulrich & Gotelli (2010). This is the slower
#'but reliable R implementation.
#'
#'@param m Community data matrix with species as rows and sites as columns. Can
#'contain either incidences (1/0) or abundances (natural numbers).


null_step_CA_IT <- function(m)
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

  # A CORRECTION FOR RANDOMIZED MATRICES WITH 0 ROWS
    # which species have 0 individuals
    zero.spec <- which(rowSums(res) == 0)
    # where will random 1 placed?
    one.loc <- sample(1:ncol(m), size = sum(zero.spec), replace=TRUE)
    res[zero.spec, one.loc] <- 1
    return(res)

  return(res)
}

# TESTING
#m <- data.Ulrich[[1]]; rowSums(null_step_CA_IT(m))

################################################################################
#'One realization of the an algorithm that resamples each row according to a Poisson distribution
#'
#'@description First, for each species in the matrix, mean abunance is calculated.
#'Then N samples are taken for each species from a Poisson distribution, with lambda
#'equal to the mean abunance.
#'
#'@param m Community data matrix with species as rows and sites as columns. Can
#'contain either incidences (1/0) or abundances (natural numbers).


null_step_CA_Poisson <- function(m)
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

#m <- data.Ulrich[[1]]; rowSums(null_step_CA_Poisson(m))


################################################################################
#'One realization of the sim2 algorithm from EcoSimR
#'
#' @description https://cran.r-project.org/web/packages/EcoSimR/vignettes/CoOccurrenceVignette.html
#'
#' @param m Community data matrix with species as rows and sites as columns. Can
#' contain either incidences (1/0) or abundances (natural numbers).


null_step_C_sim2 <- function(m)
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
#' @param metric Character string. Name of spasm's pairwise association metric.
#' @param N.sim Number of algirthm runs.
#' @export
#'
null_model_Z_score <- function(m, algorithm, metric, N.sim)
{
  res <- array(dim=c(nrow(m), nrow(m), N.sim))

  for(i in 1:N.sim)
  {
    # radnomize using a given algorithm
    null.m <- do.call(algorithm, list(m))
    # calculate the metric
    null.metric <- as.matrix(do.call(metric, list(null.m)))
    res[,,i] <- null.metric
  }

  # calculate the Z-score
  obs <- do.call(metric, list(m))
  mean.null <- as.dist(apply(X = res, MARGIN = c(1,2), FUN = mean))
  sd.null <- as.dist(apply(X = res, MARGIN = c(1,2), FUN = sum))

  Z <- (obs - mean.null) / sd.null
  return(Z)
}

# TESTING
#m <- data.Atmar[[1]]
#x <- null_model_Z_score(m,
#                        algorithm="null_step_C_sim2",
#                        metric="C_jacc", N.sim = 4)
#xobs <- CA_hell(m)
#plot(x, xobs)
