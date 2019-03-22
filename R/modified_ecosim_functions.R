#'Co-Occurrence Null model - modified verstion from EcoSimR
#'
#'@description Create a Co-Occurrence null model. This is a modified version of the function
#'"cooc_null_model"
#'from EcoSimR, with no restrictions on the metric that is applied to the community matrix.
#'
#'@param speciesData a dataframe in which rows are species, columns are sites,
#' and the entries indicate the absence (0) or presence (1) of a species in a
#' site. Empty rows and empty columns should not be included in the matrix.
#'@param algo the algorithm to use, must be "sim1", "sim2", "sim3", "sim4", "sim5", "sim6", "sim7", "sim8", "sim9", "sim10"; default is "sim9".
#'@param metric the metric used to calculate the null model: the choice is unrestricted
#'@param nReps the number of replicate null assemblages to create; default is 1000 replicates.
#'@param saveSeed TRUE or FALSE. If TRUE the current seed is saved so the simulation can be repeated; default is FALSE.
#'@param burn_in The number of burn_in iterations to use with the simFast algorithm; default is 500 burn-in replicates.
#'@param algoOpts a list containing all the options for the specific algorithm you want to use.  Must match the algorithm given in the `algo` argument.
#'@param metricOpts a list containing all the options for the specific metric you want to use.  Must match the metric given in the `metric` argument.
#'@param suppressProg TRUE or FALSE. If true, display of the progress bar in the console is suppressed; default is FALSE. This setting is useful for creating markdown documents with `knitr`.
#'@import EcoSimR
#'@export

cooc_null_model_spasm <- function(speciesData,
                                  algo = "sim9",
                                  metric,
                                  nReps = 1000,
                                  saveSeed = FALSE,
                                  burn_in = 500,
                                  algoOpts = list(),
                                  metricOpts = list(),
                                  suppressProg = FALSE,
                                  summaryF){

  params <- list(speciesData = speciesData,
                 algo = algo,
                 metric = metric,
                 nReps = nReps,
                 saveSeed = saveSeed,
                 burn_in = burn_in,
                 suppressProg = suppressProg,
                 summaryF = summaryF)

  output <- do.call(sim9_spasm, params)
  class(output) <- "coocnullmod"
  return(output)
}


# ------------------------------------------------------------------------------
# modified sim9 algorithm to handle spasm's functions

sim9_spasm <- function (speciesData,
                  algo,
                  metric,
                  nReps = 1000,
                  saveSeed = FALSE,
                  burn_in = 0,
                  algoOpts = list(),
                  metricOpts = list(),
                  suppressProg = TRUE,
                  summaryF)
{

  if(saveSeed){
    randomSeed <- .Random.seed
  } else {
    randomSeed <- NULL
  }

  ## Convert to matrix for type consistency
  if(!is.matrix(speciesData)){ speciesData <- as.matrix(speciesData)}

  ### Check for row names hidden in the data frame and automagically strip them.

  if(suppressWarnings(is.na(as.numeric(speciesData[2,1])))){
    speciesData <- speciesData[,-1]
    class(speciesData) <- "numeric"
  }

  Start.Time <- Sys.time()
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

  Sim <- simulated.metric
  End.Time <- Sys.time()
  Elapsed.Time <- format(End.Time-Start.Time,digits=2)
  Time.Stamp <- date()

  sim9.fast.out <- list(Obs=Obs,Sim=Sim, Elapsed.Time=Elapsed.Time, Time.Stamp=Time.Stamp,Metric = metric, Algorithm = algo, N.Reps = nReps, SaveSeed = saveSeed, RandomSeed = randomSeed,Randomized.Data = msim , Data = speciesData,burn.in = burn_in,burn.in.metric= burn.in.metric)
  # plot to screen the trace function for the burn in

  class(sim9.fast.out) <- "nullmod"
  return(sim9.fast.out)
}
