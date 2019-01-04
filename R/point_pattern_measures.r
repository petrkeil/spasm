
################################################################################
#' Bivariate nearest neighbor distance function (D) and K-function (K)
#' 
#' @param X Object of class 'ppp' (from 'spatstat' packge). This is a marked 
#' point pattern in which each mark gives identity of a species. Thus, the object
#' must have at least 2 species.
#' @param spec.1 Character string indicating name of species 1 in the 'ppp' object.
#' @param spec.2 Character string indicating name of species 2 in the 'ppp' object.
#' @param n.iter Number of iterations for simulation from a homogeneous Poisson
#' point process.
#' @param radius Radius for which the D and K functions should be calculated 
#' 
#' @return A list.
#' @import spatstat
#' @references Wiegand T. et al. (2007) Species associations in a heterogeneous
#' Sri Lankan dipterocarp forest. American Naturalist 170: E77-E95.
#' @references Wang X. et al. (2010) Species associations in an old-growth temperate
#' forest in north-eastern China. Journal of Ecology, 95: 674-686.
#' @references Wiegand T. & Moloney K. A. (2014) Handbook of spatial point pattern
#' analysis in ecology. CRC press.
#' @export
#' 

KD <- function(X, spec.1, spec.2, n.iter, radius)
{
  # extract the 2 focal species from the marked point pattern
  X1 <- split(X)[[spec.1]]
  X2 <- split(X)[[spec.2]]
  
  # empty container for the simulation results
  D12.sims <- K12.sims <- list()
  
  # estimate lambda of the CSR model in X2
  lambda.X2 <- exp(ppm(X2)$coef)
  
  
  for(i in 1:n.iter)
  {
    # randomize the pattern of species 2
    repeat{
      X2.rnd <- rpoispp(lambda.X2, win = X2$window)
      if(X2.rnd$n > 0) break # only return X2.rnd if there are n > 0 points
    }
    
    X1_X2.rnd <-  suppressWarnings(superimpose(X1 = X1, X2.rnd = X2.rnd)) 
    
    # the bivariate K function
    K12.r <- Kcross(X1_X2.rnd, i = "X1", j = "X2.rnd", 
                    r = c(0, radius), correction = "border")
    
    # the bivariate D function
    D12.r <- Gcross(X1_X2.rnd, i = "X1", j = "X2.rnd", 
                    r = c(0, radius), correction = "rs")
    
    K12.sims[[i]] <- K12.r$border[2]
    D12.sims[[i]] <- D12.r$rs[2]
    cat("*")
  }
  cat("\n")
  
  # summarize the simulations by quantiles
  K12.sims <- quantile(do.call(c, K12.sims), probs = c(0.025, 0.5, 0.975))
  names(K12.sims) <- c("K12sim.2.5", "K12sim.50", "K12sim.97.5")
  D12.sims <- quantile(do.call(c, D12.sims), probs = c(0.025, 0.5, 0.975))
  names(D12.sims) <- c("D12sim.2.5", "D12sim.50", "D12sim.97.5")
  
  # the observed bivariate K-function
  K12.r.obs <- Kcross(X = X, i = spec.1, j = spec.2, r = c(0, radius), 
                      correction = "border")$border[2]
  
  # the observed bivariate D-function
  D12.r.obs <- Gcross(X = X, i = spec.1, j = spec.2, r = c(0, radius), 
                      correction = "rs")$rs[2]
  
  # calcuate the P and M values
  P <- D12.r.obs - D12.sims 
  names(P) <- c("P2.5", "P.50", "P97.5")
  M <- log(K12.r.obs) - log(K12.sims)
  names(M) <- c("M2.5", "M.50", "M97.5") 
  
  # coerce results to a list
  all.res <- list(spec.1 = spec.1, spec.2 = spec.2,
                  radius = radius,
                  P = P,
                  M = M,
                  K12 = c(K12obs = K12.r.obs, K12.sims), 
                  K12.signif = unname(K12.r.obs < K12.sims[1] | 
                                        K12.r.obs > K12.sims[3]),
                  D12 = c(D12obs = D12.r.obs, D12.sims),
                  D12.signif = unname(D12.r.obs < D12.sims[1] | 
                                        D12.r.obs > D12.sims[3]))
  return(all.res)
}


################################################################################
#' Multi-species P-M classfication of Wiegand et al. (2007)
#' 
#' @param X Object of class 'ppp' (from 'spatstat' packge). This is a marked 
#' point pattern in which each mark gives identity of a species. Thus, the object
#' must have at least 2 species.
#' @param n.iter Number of iterations for simulation from a homogeneous Poisson
#' point process.
#' @param r Radius for which the P and M values should be calculated.
#' 
#' @return A list.
#' @import spatstat
#' @references Wiegand T. et al. (2007) Species associations in a heterogeneous
#' Sri Lankan dipterocarp forest. American Naturalist 170: E77-E95.
#' @references Wang X. et al. (2010) Species associations in an old-growth temperate
#' forest in north-eastern China. Journal of Ecology, 95: 674-686.
#' @references Wiegand T. & Moloney K. A. (2014) Handbook of spatial point pattern
#' analysis in ecology. CRC press.
#' @export
#' 


PM.multispec <- function(X, n.iter, radius)
{
  spec.names <- as.character(unique(X$marks))
  spec.combos <- t(combn(spec.names, m = 2))
  
  res <- list()
  for(i in 1:nrow(spec.combos))
  {
    cat(paste(spec.combos[i,1], spec.combos[i,2], sep="-")); cat("\n")
    KD.i <- KD(X = X, 
               spec.1 = spec.combos[i,1], 
               spec.2 = spec.combos[i,2], 
               n.iter = n.iter, radius = radius)
    
    res.i <- data.frame(spec.1 = spec.combos[i,1],
                        spec.2 = spec.combos[i,2],
                        data.frame(t(KD.i$P)),
                        data.frame(t(KD.i$M)))
    res[[i]] <- res.i
  }
  res <- do.call(what = rbind, args = res)
  return(res)
}


################################################################################
#' Bivariate nearest neighbor distance function (D) and K-function (K)
#' 
#' @param PM A list object returned by the PM.multispec function.
#' 
#' @return A ggplot figure.
#' @import spatstat
#' @import ggplot2
#' @export
#' 

plot.PM <- function(PM, xlim = NULL, ylim = NULL, 
                    labs = FALSE, lab.shift = NULL,
                    main = NULL)
{
  require(ggplot2)
  
    PM <- data.frame(labels = paste(PM[,1], PM[,2], sep="-"),
                     PM)
    pp <- ggplot(data = PM, aes(x = P.50, y = M.50)) +
      geom_hline(yintercept = 0, colour="darkgrey") + 
      geom_vline(xintercept = 0, colour="darkgrey") +
      geom_pointrange(aes(x = P.50, y = M.50,
                          ymin = M2.5, ymax = M97.5)) +
      geom_errorbarh(aes(xmin = P2.5, xmax = P97.5)) +
      xlab("Classification axis P") +
      ylab("Classification axis M") +
      geom_point() +
      theme_bw()
    
    if(labs)
    {
      pp <- pp + geom_text(data = PM, 
                           aes(x = P.50, 
                               y = M.50, 
                               label = labels), 
                           nudge_x = lab.shift,
                           nudge_y = lab.shift)

    }
    
    if(is.null(main) == FALSE) { pp <- pp + ggtitle(main)  }
    
    if(length(xlim) == 2) { pp <- pp + xlim(xlim[1], xlim[2]) }
    if(length(ylim) == 2) { pp <- pp + ylim(ylim[1], ylim[2]) }
    return(pp)
}


################################################################################
#' Bivariate pair correlation function for a given distance r
#' 
#' @param X Object of class 'ppp' (from 'spatstat' packge). This is a marked 
#' point pattern in which each mark gives identity of a species. Thus, the object
#' must have at least 2 species.
#' @param spec.1 Character string indicating name of species 1 in the 'ppp' object.
#' @param spec.2 Character string indicating name of species 2 in the 'ppp' object.
#' @param n.iter Number of iterations for simulation from a homogeneous Poisson
#' point process.
#' @param radius.seq Numeric vector. Sequence of radii for which the pcf function 
#' should be calculated 
#' 
#' @return A list.
#' @import spatstat
#' @references Wiegand T. et al. (2007) Species associations in a heterogeneous
#' Sri Lankan dipterocarp forest. American Naturalist 170: E77-E95.
#' @references Wang X. et al. (2010) Species associations in an old-growth temperate
#' forest in north-eastern China. Journal of Ecology, 95: 674-686.
#' @references Wiegand T. & Moloney K. A. (2014) Handbook of spatial point pattern
#' analysis in ecology. CRC press.
#' @export
#' 

pcf.r <- function(X, spec.1, spec.2, n.iter, radius.seq)
{
  # extract the 2 focal species from the marked point pattern
  X1 <- split(X)[[spec.1]]
  X2 <- split(X)[[spec.2]]
  
  ### the observed bivariate pcf
  # ----------------------------------------------------------------------------    
  pcf12.obs <- pcfcross(X = X, 
                        i = spec.1, 
                        j = spec.2, 
                        r = radius.seq, 
                        correction = "translate")
  
  ### pcf for complete spatial independence under homogeneous Poisson pattern
  # ----------------------------------------------------------------------------    
  pcf12.hom <- envelope(Y = X, 
                        i = spec.1, 
                        j = spec.2, 
                        r = radius.seq,
                        fun = pcfcross, 
                        nsim = n.iter, 
                        correction="translate",
                        savefuns = TRUE)
  pcf12.hom <- data.frame(attr(pcf12.hom, "simfuns"))
  # remove the first r column
  pcf12.hom <- pcf12.hom[,-1]
  pcf12.hom.summary <- data.frame(t(apply(pcf12.hom, 
                                          MARGIN = 1, 
                                          FUN = quantile, 
                                          probs = c(0.025, 0.5, 0.975))))
  names(pcf12.hom.summary) <- c("pcf12hom.2.5", "pcf12hom.50", "pcf12hom.97.5")
  
  
  ### pcf for heterogeneous Poisson pattern
  # ----------------------------------------------------------------------------    
  # empty container for the simulation results
  pcf12.sims <- list()
  
  # estiamte optimal bandwidth for the kernel density 
  sigma = bw.CvL(X2)
  
  # estimate lambda of the CSR model in X2
  lambda.X2 <- density.ppp(X2, kernel = "gaussian", sigma = sigma)
  lambda.X2[lambda.X2 < 0] <- 0 # get rid of any negative intensities
  
  cat(paste("Generating", n.iter, "simulation of heterogenous poisson process ...\n" ))
  for(i in 1:n.iter)
  {
    # randomize the pattern of species 2
    repeat{
      X2.rnd <- rpoispp(lambda.X2) 
      if(X2.rnd$n > 0) break # only return X2.rnd if there are n > 0 points
    }
    
    X1_X2.rnd <-  suppressWarnings(superimpose(X1 = X1, X2.rnd = X2.rnd)) 
    
    # the bivariate Pair Correlation Function
    pcf12.r <- pcfcross(X1_X2.rnd, 
                        i = "X1", 
                        j = "X2.rnd", 
                        r = radius.seq, 
                        correction = "translate")
    
    pcf12.sims[[i]] <- pcf12.r$trans
    
    cat(paste(i, ", ", sep=""))
  }
  cat("\nDone\n")
  
  # convert the simulation list to data.frame
  sims.dat <- data.frame(pcf12.sims); names(sims.dat) <- 1:n.iter
  
  # summarize the simulations by quantiles
  pcf12.het.summary <- data.frame(t(apply(sims.dat, 
                                          MARGIN = 1, 
                                          FUN = quantile, 
                                          probs = c(0.025, 0.5, 0.975))))
  names(pcf12.het.summary) <- c("pcf12het.2.5", "pcf12het.50", "pcf12het.97.5")
  
  
  ### coerce correlation results to one data frame
  res <- data.frame(radius = radius.seq, 
                    pcf12obs = pcf12.obs$trans, 
                    pcf12.hom.summary,
                    pcf12.het.summary)
  
  # remove the first line with r=0 and Inf pcf
  res <- res[-1,]
  
  # coerce results to a list
  all.res <- list(spec.1 = spec.1, spec.2 = spec.2,
                  sigma = unname(sigma),
                  simulations = res)
  return(all.res)
}




################################################################################


plot.pp.pcf <- function(x)
{
  col <- c("#e41a1c","#377eb8")
  alpha = 0.4
  
  p <- ggplot(data = x$simulations, aes(x = radius, y = pcf12het.50)) +
          # homogeneous
          geom_ribbon(aes(ymin = pcf12hom.2.5, ymax = pcf12hom.97.5), fill=col[1], 
                      alpha = alpha) +
          geom_line(aes(x = radius, y = pcf12hom.50), colour = col[1]) +
          geom_point(aes(x = radius, y = pcf12hom.50), colour = col[1]) +
          # heterogeneous
          geom_ribbon(aes(ymin = pcf12het.2.5, ymax = pcf12het.97.5), fill=col[2], 
                      alpha = alpha) +
          geom_line(colour = col[2]) + 
          geom_point(colour = col[2]) +
          # obsered
          geom_line(aes(x = radius, y = pcf12obs), colour = "black") + 
          geom_point(aes(x = radius, y = pcf12obs), colour = "black") +
          ylab("Bivariate PCF") +
          ggtitle(paste(x$spec.1, " - ", x$spec.2)) +
          theme_bw() 
  p
}


