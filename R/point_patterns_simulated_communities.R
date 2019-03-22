#' Bivariate probablity density function in a 2D square

#' @param x x Coordinates of points whose PD should be evaluated.
#' @param y y Coordinates of points whose PD should be evaluated.
#' @param var Variance of the PDF, which is the same in x and y direction. Note: covariance is 0.
#' @param x.centr x coordinate of the mean in the 2D space
#' @param y.centr y coordinate of the mean in the 2D space
#' @import spatstat
#' @import mvtnorm

dpoint.MVN <- function(x, y, var, x.centr, y.centr)
{
  sigma = matrix(c(var, 0, 0, var), nrow=2, ncol=2)
  dmvnorm(x=cbind(x,y), mean = c(x.centr, y.centr), sigma)
}

# Test:
# dpoint.MVN (0.1, 0.1, var=0.1, x.centr=0.1, y.centr=0.1)


# ------------------------------------------------------------------------------
#' Funciton that generates the MVN density as a pixel image
#'
#' @inheritParams dpoint.MVN

dpoint.MVN.image <- function(var, x.centr, y.centr)
{
  W <- square()
  Z <- as.im(dpoint.MVN, var = var, x.centr = x.centr, y.centr = y.centr, W)
  return(Z)
}

# Test:
# plot(dpoint.MVN.image(var=0.1, 0.5, 0.5))


# ------------------------------------------------------------------------------
#' Function that generates n points in a 2D space with a bivariate normal PDF,
#' given the variance, and xy coordinates of the centre of the PDF.

#' @inheritParams dpoint.MVN
#' @param n number of points

rpoint.MVN <- function(n, var, x.centr, y.centr)
{
  spec <- rpoint(n = n, f = dpoint.MVN ,
                 var = var, x.centr = x.centr, y.centr=y.centr)
  return(spec)
}

# Test:
# a <- rpoint.MVN(n = 100, var = 1, x.centr=0.5, y.centr=0.5); plot(a)

# ------------------------------------------------------------------------------
#' Function that simulates a "centered" community

#' The simulation proceeds in 2 steps:
#' 1. Genreate set of 'mother points' from the bivariate normal, with a given var.intersp.
#' 2. Treat each 'mother point' as a centre of a cluster of points of a given species,
#'    with var.consp as the "spread" of the cluster.

#' @param abund.vect vector of (integer) abundances of species in a community
#' @param var.consp conspecific variance, i.e. spread of points in a single species
#' @param var.intersp interspecific variance, i.e. spread of the mother points
#' @import spatstat
#' @import mvtnorm
#' @export

sim.comm.centered <- function(abund.vect, var.consp, var.intersp)
{
  sp.centres <- rpoint.MVN(n=length(abund.vect),
                           var=var.intersp,
                           x.centr=0.5, y.centr=0.5)

  sp.centres.coords <- coords(sp.centres)

  # generate density for each species, with its unique center
  comm <- list()
  for(i in 1:length(abund.vect))
  {
    sp <- rpoint.MVN(n = abund.vect[i],
                     var = var.consp,
                     x.centr = sp.centres.coords[i,1],
                     y.centr = sp.centres.coords[i,2])
    comm[[i]] <- data.frame(coords(sp), species = i)
  }
  comm <- do.call("rbind", comm)
  comm <- ppp(x = comm$x,
              y=comm$y,
              marks = as.factor(paste("sp", comm$species, sep="")),
              window = square())
  return(comm)
}

# Test:
# a <- sim.comm.centered(abund.vect  = c(2,4, 8, 16, 32, 64, 128),
#                       var.consp   = 0.01,
#                       var.intersp = 0.1); plot(split(a))


# ------------------------------------------------------------------------------
#' Function that simulates a "jittered" community
#'
#' The simulation proceeds as:
#' 1. Generate point pattern for the most abundant (master) species,
#'    with a given MV variance 'var.consp'
#' 2. Generate smoothed density of the master points, with a
#'    given gaussian kernel width 'var.intersp'
#' 3. For each other species in the vector, draw points from the smoothed
#'    density surface and a given abundance
#' @param abund.vect vector of (integer) abundances of species in a community
#' @param var.consp conspecific variance, i.e. spread of points in the master species
#' @param var.intersp interspecific variance, i.e. width of the density kernel
#' @import spatstat
#' @export

sim.comm.jittered <- function(abund.vect, var.consp, var.intersp)
{
  # spread the points of the most abundant "master" species
  master.species <- rpoint.MVN(n=max(abund.vect),
                               var=var.consp,
                               x.centr=0.5, y.centr=0.5)
  # generate desnity surface of the master points, with a kernel width
  # that is given by var.intersp
  master.dens <- density(master.species,
                         sigma = var.intersp,
                         kernel="gaussian",
                         positive = TRUE)

  comm <- list()
  comm[[1]] <- data.frame(coords(master.species), species = 1)

  abund.vect <- sort(abund.vect, decreasing = TRUE)

  # for each species except for the most abundant one
  for(i in 2:(length(abund.vect)))
  {
    # generate random points with a given abundance and with master density
    sp <- rpoint(n = abund.vect[i],
                 f = master.dens)
    comm[[i]] <- data.frame(coords(sp), species = i)
  }
  comm <- do.call("rbind", comm)
  comm <- ppp(x = comm$x,
              y=comm$y,
              marks = as.factor(paste("sp", comm$species, sep="")),
              window = square())
  return(comm)
}


# Test:
# a <- sim.com.jittered(abund.vect  = c(2,4, 8, 16, 32, 64, 128),
#                       var.consp   = 0.1,
#                       var.intersp = 0.1); plot(split(a))


# ------------------------------------------------------------------------------


fractal.norm.raster <- function(grains = c(64, 32, 16, 8, 4, 2), seed = FALSE)
{

  # block the random number generator
  if(seed){ set.seed(seed)}

  sigma <- 1

  # initialize the algorithm at the coarsest resolution
  coarsest.grain <- min(grains)
  N <- coarsest.grain^2
  M = raster(matrix(rnorm(n=N, mean=0, sd = sigma), coarsest.grain, coarsest.grain))

  res <- list(M)

  for(i in 2:length(grains))
  {
    M <- res[[i-1]]
    M <- disaggregate(M, fact = 2)
    N <- nrow(M)*nrow(M)
    X <- M
    X[] <- rnorm(N, mean = M[], sd = sigma)
    res[[i]] <- X
  }

  X <- res[[length(grains)]]
  X <- as.im(as.matrix(X), W = square())
  # make everything on a positive scale and between 0 and 1
  X <- (X - min(X)) / max(X - min(X))
  return(X)
}

# Example:
# x <- fractal.norm.raster(); par(mfrow=c(1,2)); hist(x); plot(x)



sim.comm.fractal <- function(abund.vect, var.consp, var.intersp)
{
  require(truncnorm)

  # generate the fractal landscape
  x <- fractal.norm.raster(seed = 222)

  # sample niche optima for all species from a truncated normal distribution
  optima <- rtruncnorm(n=length(abund.vect),
                       a=0, b=1,
                       mean = 0.7, sd = var.intersp) # NOTE THE FIXED MEAN!!!!!

  comm <- list()

  # for each species
  for(i in 1:(length(abund.vect)))
  {
    habitat.suitability <- x
    habitat.suitability[] <- dnorm(x[], mean = optima[i], sd = var.consp)

    # generate random points with a given abundance and density given
    # by the habitat suitability
    sp <- rpoint(n = abund.vect[i],
                 f = habitat.suitability)
    comm[[i]] <- data.frame(coords(sp), species = i)
  }

  comm <- do.call("rbind", comm)
  comm <- ppp(x = comm$x,
              y=comm$y,
              marks = as.factor(paste("sp", comm$species, sep="")),
              window = square())
  return(comm)
}

# a <- sim.com.fractal(abund.vect  = c(2,4, 8, 16, 32, 64, 128),
#                      var.consp   = 0.001,
#                      var.intersp = 0.001)
# plot(split(a))

# ------------------------------------------------------------------------------
# Function that converts the marked ppp object generated by multi.spec.MVN (above)
# to a siteXspec matrix. Grain is the number of grid cells along side of the
# big square.

ppp.to.siteXspec <- function(comm.ppp, grain)
{
  comm.ppp <- split(comm.ppp)
  comm.im  <- lapply(comm.ppp, FUN = as.im, dimyx = c(grain, grain))
  comm.num <- lapply(comm.im,  FUN = as.numeric)
  siteXspec<- do.call("cbind", comm.num)
  return(siteXspec)
}

