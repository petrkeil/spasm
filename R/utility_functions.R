#' Conversion of point pattern to community matrix Y
#'
#' @param comm.ppp A marked point pattern of spatstat's 'ppp', where marks identify each
#' point (individual) to a species.
#' @param dim.yx Number of cells along the Y and X dimension of the raster.
#' @param plotstack Logical. Should the raster for each species be plotted?
#' @param clean Logical. Should only sites with positive abundance of at least one
#' species be returned?
#'
#' @return A list containing the raster stack with rasters for each species,
#' the abundance matrix Y, spatial coordinates of each site (grid cell), and
#' number of cells.
#' @import spatstat
#' @import raster
#' @export
#'

ppp.to.comm <- function(comm.ppp, dim.yx, plotstack = FALSE, clean = FALSE)
{
  require(raster)
  pixellated <- lapply(X = split(comm.ppp), FUN = pixellate, dimyx = dim.yx)
  rasterized <- lapply(pixellated, FUN = raster)
  final.stack <- stack(rasterized)
  if(plotstack) plot(final.stack, box=FALSE, frame = FALSE, axes=FALSE)

  abundance <- as.data.frame(final.stack)
  coords <- data.frame(coordinates(final.stack))

  N.nonempty <- ncell(final.stack) - sum(is.na(final.stack[[1]][]))

  if(clean)
  {
    good <- (rowSums(abundance) == 0 | is.na(rowSums(abundance))) == FALSE
    abundance <- abundance[good,]
    coords <- coords[good,]
  }

  res <- list(rasters = final.stack,
              abundance = t(abundance),
              coords = coords,
              N.cells = N.nonempty)
  return(res)
}

# ------------------------------------------------------------------------------
#' Function that removes, from point pattern, species with N lower than a specificed threshold
#'
#' @param comm.ppp A marked point pattern of spatstat's 'ppp', where marks identify each
#' point (individual) to a species.
#' @param min.abu The threshold. All species with N >= min.abu will be removed.
#' @return ppp object
#' @import spatstat
#' @import raster
#' @export
#'

ppp.trim <- function(comm.ppp, min.abu)
{
  smr <- summary(comm.ppp)$marks
  good.sp <- rownames(smr)[smr$frequency >= min.abu]
  subset.ppp(comm.ppp, marks %in% good.sp, drop = TRUE)
}

# ------------------------------------------------------------------------------
ppp.plot.2.spec <- function(comm.ppp, spec.1, spec.2)
{
  X1 <- split(comm.ppp)[[spec.1]]
  X2 <- split(comm.ppp)[[spec.2]]
  plot(X1, main = paste(spec.1, spec.2, sep = " - "))
  plot(X2, add=T, col="red")
}

#plot.2.spec.ppp (comm.ppp = Bsubp,
#                 spec.1 = "Ambrosia deltoidea",
#                 spec.2 = "Haplopappus tenuisectus")


# ------------------------------------------------------------------------------
#' Conversion of distance matrix to list

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
