library(ggplot2)
library(dplyr)
library(vegan)
library(tidyr)
library(gridExtra)
library(raster)
library(spasm)
library(spatstat)

# ------------------------------------------------------------------------------
# ILLUSTRATIONS

alpha.vec <- c(20, 8, 0, -8, -20)
params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1),
                      alpha = alpha.vec)

res.pts <- list()

for(i in 1:nrow(params))
{
  set.seed(123)
  a <- sim.pair(abund.vect = c(100, 100),
                var.consp = params$var.consp[i],
                alpha = params$alpha[i],
                plot.comm = FALSE)

  res.pts[[i]] <- data.frame(x = a$x,
                             y = a$y,
                             Species= as.character(a$marks),
                             ISA = params$alpha[i],
                             CSA = params$var.consp[i])
}
res.pts <- do.call("rbind", res.pts)


pp.pairs <-  ggplot(data = res.pts, aes(x = x, y = y)) +
              geom_point(aes(colour = Species), shape = 19, alpha = 0.5) +
              facet_grid(CSA ~ ISA, labeller = "label_both") +
              theme_bw() + theme(panel.grid = element_blank()) +
              theme(#legend.position = "bottom",
                    legend.position = c(0.95, 0.78),
                    legend.background = element_rect(fill = "lightgray")) +
             scale_x_continuous(breaks = c(0, 1)) +
             scale_y_continuous(breaks = c(0, 1)) +
            labs(title = "(a)")

print(pp.pairs)



res <- list()
for(i in 1:length(alpha.vec))
{
  d <- seq(0, 1, by = 0.01)
  a <- PDFtexp(d,
               alpha = alpha.vec[[i]],
               dlim =c(0, 1))
  res[[i]] <- data.frame(Distance = d, pdf = a, ISA = alpha.vec[i])
}
res <- do.call("rbind", res)

pdfs <- ggplot(data = res, aes(x = Distance, y = pdf)) +
  geom_line() +
  facet_grid(.~ISA, labeller = "label_both") + theme_bw() +
  theme(panel.grid = element_blank(),
        plot.margin=unit(c(5.5,25,5.5,4), "points")) +
  scale_x_continuous(breaks = c(0, 1)) +
  ylab(expression(f[sp2](D))) + xlab("D") +
  labs(title = "(b)")
pdfs

png("figures/simulations_ISA_vs_CSA_examples.png",
    width=2100, height = 2000, res=250 )

grid.arrange(pp.pairs, pdfs, ncol=1, nrow=2, heights = c(1, 0.4))

grid::grid.text("repulsion, seggregation", x = 0.23, y = 0.98)
grid::grid.text("independence", x = 0.5, y = 0.98)
grid::grid.text("attraction", x = 0.78, y = 0.98)

grid::grid.text("repulsion, seggregation", x = 0.23, y = 0.267)
grid::grid.text("independence", x = 0.5, y = 0.267)
grid::grid.text("attraction", x = 0.78, y = 0.267)

dev.off()

# ------------------------------------------------------------------------------
# Plotting simulation procedure

png("figures/sumulation_steps.png", width = 2000, height = 1200, res= 300 )

par(mfrow=c(2,3), mai = c(0.3, 0.3, 0.3, 0.3), adj=0)

  set.seed(123)
  x <- runif(1, 0, 1)
  y <- runif(1, 0, 1)
  sigma <- 0.01
  alpha = -5

  mu <- ppp(x,y)
  plot(mu, main = expression(paste("(a) ", mu)))

  f.sp1 <- dpoint.MVN.image(var=sigma, x.centr=x, y.centr=y)
  f.sp1 <- f.sp1/sum(f.sp1)
  plot(f.sp1, main = expression(paste("(b) ", f[sp1](mu, Sigma) )))
  contour(f.sp1, add=T, col="white")

  sp1 <- rpoint.MVN(n=100,
                    var=sigma,
                    x.centr=x, y.centr= y)
  plot(sp1, main = expression(paste("(c) sp1")))

  sp1.dist <- sp1.prob <- distmap(sp1)
  sp1.prob[] <- PDFtexp(sp1.dist[], alpha)
  sp1.prob <- sp1.prob/sum(sp1.prob)

  plot(sp1.dist, main = expression("(d) D"))
  contour(sp1.dist, add=T, col="white")

  plot(sp1.prob, main = expression(paste("(e) ", f[sp2], "(D)")))
  contour(sp1.prob, add=T, col="white")

  sp2 <- rpoint(n = 100, f = sp1.prob)
  plot(sp2, main = expression(paste("(c) sp2")))

  #plot(sp1); plot(sp2, add=TRUE, col= "red")

dev.off()

# ------------------------------------------------------------------------------
# SIMULATIONS - VARIOGRAMS

library(vario)
library(devtools)
#remove.packages("vario", lib="~/R/x86_64-pc-linux-gnu-library/3.5")
#install('/home/pk33loci/Dropbox/Interspecific_aggregation/vario')
 library.dynam(chname="vario", package = "vario", lib.loc=.libPaths())

params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1),
                      alphas = alpha.vec)

res.vario <- list()
for(i in 1:nrow(params))
{
  set.seed(12345)
  comm <- sim.pair(abund.vect = c(200, 200),
                   var.consp = params$var.consp[i],
                   alpha = params$alphas[i],
                   plot.comm = FALSE)

  comm <- ppp.to.comm(comm, dim.yx = c(20, 20))
  abu <- t(comm$abundance)

  # remove double zeroes
  good.cells <- rowSums(abu) != 0
  abu <- abu[good.cells,]
  xy <- comm$coords
  xy <- xy[good.cells,]

  v.obs <- vario(x = as.matrix(abu),
                 coord = xy,
                 grain = 0.05,
                 pos.neg = TRUE,
                 hmax = 1,
                 binary = FALSE)

  pos <- data.frame(Dist = v.obs$vario$Dist,
                    Covariance = v.obs$vario$pos,
                    Type = "Positive")
  neg <- data.frame(Dist = v.obs$vario$Dist,
                    Covariance = abs(v.obs$vario$neg),
                    Type = "| Negative |")
  posneg <- rbind(neg, pos)

  res.vario[[i]] <- data.frame(posneg,
                               ISA = params$alphas[i],
                               CSA = params$var.consp[i])
}

res.vario <- do.call(rbind, res.vario)

png("figures/simulations_vario.png",
    width=2400, height = 1300, res=250 )
ggplot(data=res.vario, aes(x=Dist, Covariance)) +
  geom_line(aes(linetype = Type)) +
  #geom_point(aes(colour = Sign, shape= Sign)) +
  scale_colour_brewer(palette= "Set1") +
  facet_grid(CSA~ISA, scales = "free_y", labeller = "label_both") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        #legend.position = "none",
        legend.background = element_rect(fill = "lightgray"),
        plot.margin=unit(c(20,5.5,5.5,5.5), "points")) +
  scale_x_continuous(breaks = c(0, 1), limits=c(0,1)) +
  xlab("Distance")

grid::grid.text("repulsion, seggregation", x = 0.2, y = 0.98)
grid::grid.text("independence", x = 0.44, y = 0.98)
grid::grid.text("attraction", x = 0.68, y = 0.98)
dev.off()


################################################################################
# POINT PATTERN BIVARIATE PCFs

params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1),
                      alphas = alpha.vec)

distance <- seq(0,1, by=0.01)

res.pcf <- list()
for(i in 1:nrow(params))
{
  set.seed(12345)
  comm <- sim.pair(abund.vect = c(200, 200),
                   var.consp = params$var.consp[i],
                   alpha = params$alphas[i],
                   plot.comm = FALSE)

  pcf <- pcfcross(X = comm,
                  i = "sp1",
                  j = "sp2",
                  correction="translate",
                  r = distance)
  res.pcf[[i]] <- data.frame(distance, g = pcf$trans,
                             ISA = params$alphas[i],
                             CSA = params$var.consp[i])
}

res.pcf <- do.call(rbind, res.pcf)
res.pcf <- res.pcf[res.pcf$g != Inf,]
res.pcf <- res.pcf[is.na(res.pcf[,1]) == FALSE,]

png("figures/simulations_PCF.png",
    width=2400, height = 1300, res=250 )
ggplot(data=res.pcf, aes(x=distance, g)) +
  #geom_hline(yintercept=1, colour = "grey") +
  geom_line() +
  facet_grid(CSA~ISA, scales = "free_y", labeller = "label_both") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        #legend.position = "none",
        plot.margin=unit(c(20,90,5.5,5.4), "points"),
        legend.background = element_rect(fill = "lightgray")) +
  scale_x_continuous(breaks = c(0, 1), limits=c(0,1)) +
  ylab("Bivariate PCF") + xlab("Distance")

grid::grid.text("repulsion, seggregation", x = 0.2, y = 0.98)
grid::grid.text("independence", x = 0.44, y = 0.98)
grid::grid.text("attraction", x = 0.68, y = 0.98)
dev.off()

# ------------------------------------------------------------------------------







