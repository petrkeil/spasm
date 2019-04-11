################################################################################
#'The IT algorithm of Ulrich & Gotelli (2010)
#'
#'@description Uses the IT algorithm by Ulrich & Gotelli (2010). This is the slower
#'but reliable R implementation.
#'
#'@param m Community data matrix with species as rows and sites as columns. Can
#'contain either incidences (1/0) or abundances (natural numbers).
#'@export

IT.r <- function(m)
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

  return(res)
}

################################################################################
#'Abundance null model
#'
#'@description Randomizes community matrix many times to create a null distribution.
#'This function uses the IT algorithm by Ulrich & Gotelli (2010)
#'
#'@param m Community data matrix with species as rows and sites as columns. Can
#'contain either incidences (1/0) or abundances (natural numbers).
#'@param metric An abundance-based association measure. Character string. Options are names
#'of CA_ functions of this package spasm.
#'@param nReps Integer. Number of samples from the null distribution.
#'@param summaryF Name of function that should be used to summarize pairwise similarity matrices.
#'For example mean, median, sum, ...
#'@param ... Further arguments to the CA_ metric function.
#'@export

null_model_abund_spasm <- function(m,
                                   metric,
                                   nReps = 1000,
                                   summaryF,
                                   ...)
{
  # convert the string to the function object
  metric <- get(metric)
  # observed metric
  obs <- summaryF(metric(m, ...))

  sims <- list()

  pb <- txtProgressBar(min = 1, max = nReps, char=">", style=3) # initialize progress bar
  for(i in 1:nReps)
  {
    # perm.m <- IT(m.rowsums, m.colsums, m.nrow, m.ncol, N) # the C++ implementation
    perm.m <- IT.r(m) # one instance of the randomization algorithm
    sims[i] <- summaryF(metric(perm.m, ...)) # summarize the results
    setTxtProgressBar(pb, i)
  }
  close(pb)

  sims <- unlist(sims)
  std_z_score <- (obs - sims)/sd(sims) # standardized deviation from the null model

  return(list(obs=obs, sims=sims, std_z_score=std_z_score))
}

# ----------------------------------------------
# Testing:
# require(spasm); m <- data.Ulrich[[1]]
#
# meanpos <- function(x) mean(x[x<0])
# x <- null_model_abund_spasm(m, "CA_cov_cor", summaryF = mean, nReps = 1000,
#                             method="pearson", transf = "hellinger")
# require(ggplot2)
# ggplot(data = data.frame(sims=x$std_z_score), aes(x = sims)) + geom_histogram()
# ggplot(data = data.frame(sims=x$sims), aes(x = sims)) + geom_histogram()  +
# geom_vline(xintercept=x$obs, colour="red")



################################################################################
#'Cooccurrence (binary presence-absence) null model - modified verstion from EcoSimR
#'
#'@description Create a Co-Occurrence null model. This is based on the function
#'"cooc_null_model" from EcoSimR, with has been modified to have
#'no restrictions on the metric that is applied to the community matrix.
#'
#'@param speciesData a dataframe in which rows are species, columns are sites,
#' and the entries indicate the absence (0) or presence (1) of a species in a
#' site. Empty rows and empty columns should not be included in the matrix.
#'@param algo the algorithm to use, must be "sim1", "sim2", "sim3", "sim4", "sim5", "sim6", "sim7", "sim8", "sim9", "sim10"; default is "sim9".
#'@param metric the metric used to calculate the null model: the choice is unrestricted
#'@param nReps the number of replicate null assemblages to create; default is 1000 replicates.
#'@param burn_in The number of burn_in iterations to use with the simFast algorithm; default is 500 burn-in replicates.
#'@param suppressProg TRUE or FALSE. If true, display of the progress bar in the console is suppressed; default is FALSE. This setting is useful for creating markdown documents with `knitr`.
#'@import EcoSimR
#'@export

null_model_cooc_spasm <- function(speciesData,
                                  algo = "sim9",
                                  metric,
                                  nReps = 1000,
                                  burn_in = 500,
                                  suppressProg = FALSE,
                                  summaryF,
                                  ...){

  # convert the matrix to binary, if it already does not come as such
  speciesData[speciesData>1] <- 1


  if(algo == "sim2")
  {
    sim <- sim2_spasm(m = speciesData,
                      metric = metric,
                      nReps = nReps,
                      summaryF = summaryF,
                      ...)
  }

  if(algo == "sim9")
  {
    params = list(speciesData = speciesData,
           metric = metric,
           nReps = nReps,
           burn_in = burn_in,
           suppressProg = suppressProg,
           summaryF = summaryF)
    sim <- do.call(sim9_spasm, params)
  }

  return(list(obs=sim$Obs,
              sims=sim$Sim,
              z_score = sim$Sim - sim$Obs,
              std_z_score = (sim$Sim - sim$Obs)/sd(sim$Sim)))
}


# ----------------------------------------------
# Testing:
# require(spasm); m <- data.Atmar[[1]]
# x <- null_model_cooc_spasm(m, metric="C_combo", summaryF = mean, nReps = 1000)
# require(ggplot2)
# ggplot(data = data.frame(sims=x$std_z_score), aes(x = sims)) + geom_histogram()
# ggplot(data = data.frame(sims=x$sims), aes(x = sims)) + geom_histogram()  +
# geom_vline(xintercept=x$obs, colour="red")



# ------------------------------------------------------------------------------
#' Sim2 Co-occurrence Randomization Algorithm
sim2_spasm <- function(m,
                       metric,
                       nReps,
                       summaryF,
                       ...)

{
  # convert the string to the function object
  metric <- get(metric)
  # observed metric
  obs <- summaryF(metric(m, ...))

  sims <- list()

  pb <- txtProgressBar(min = 1, max = nReps, char=">", style=3) # initialize progress bar
  for(i in 1:nReps)
  {
    perm.m <-  t(apply(m,1,sample)) # one instance of the randomization algorithm 2
    sims[i] <- summaryF(metric(perm.m, ...)) # summarize the results
    setTxtProgressBar(pb, i)
  }
  close(pb)

  sims <- unlist(sims)

  return(list(Obs=obs, Sim=sims))
}

# sim2_spasm(m, metric = "C_w", nReps = 10, summaryF = mean)

# ------------------------------------------------------------------------------
#' Modified sim9 algorithm to handle spasm's functions, taken from EcoSimR

sim9_spasm <- function (speciesData,
                        metric,
                        nReps = 1000,
                        burn_in = 0,
                        suppressProg = TRUE,
                        summaryF)
{
  ## Convert to matrix for type consistency
  if(!is.matrix(speciesData)){ speciesData <- as.matrix(speciesData)}

  metricF <- get(metric)
  summaryF <- summaryF # Petr Keil: define a summary function to be applied to spasm metric

  Obs <- summaryF(metricF(speciesData))

  #Trim the matrix to be just rowssums > 0
  msim <- speciesData[rowSums(speciesData) > 0, ]
  ifelse(burn_in == 0, burn_in <- max(1000,10*nrow(msim)),burn_in <- burn_in)
  burn.in.metric <- vector(mode="numeric",length = burn_in)
  simulated.metric <- vector(mode="numeric",length = nReps)

  # run sequential swap for burn in series
  if(suppressProg){
    bi_pb <- txtProgressBar(min = 0, max = nReps, style = 3, file = stderr())
  } else{
    cat("Burn-in Progress \n")
    bi_pb <- txtProgressBar(min = 0, max = nReps, style = 3)
  }
  for (i in 1:burn_in)
  {
    msim <-sim9_single(msim)
    burn.in.metric[i] <- summaryF(metricF(msim))
    setTxtProgressBar(bi_pb, i)
  }
  close(bi_pb)
  # run sequential swap for simulated series
  if(suppressProg){
    stat_pb <- txtProgressBar(min = 0, max = nReps, style = 3, file = stderr())
  } else{
    cat("Swap Progress \n")

    stat_pb <- txtProgressBar(min = 0, max = nReps, style = 3)
  }
  for (i in 1: nReps)
  {
    msim <-sim9_single(msim)
    simulated.metric[i] <- summaryF(metricF(msim))
    setTxtProgressBar(stat_pb, i)
  }
  close(stat_pb)

  sim9.fast.out <- list(Obs=Obs, Sim=simulated.metric)
  return(sim9.fast.out)
}


