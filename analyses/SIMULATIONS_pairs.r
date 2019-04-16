library(ggplot2)
library(dplyr)
library(vegan)
library(tidyr)
library(gridExtra)
library(raster)
library(spasm)

# ------------------------------------------------------------------------------
# SIMULATIONS - ILLUSTRATIONS

alpha.vec <- c(20, 8, 0, -8, -20)
params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1),
                      alpha = alpha.vec)

res <- list()

for(i in 1:nrow(params))
{
  set.seed(123)
  a <- sim.pair(abund.vect = c(100, 100),
                var.consp = params$var.consp[i],
                alpha = params$alpha[i],
                plot.comm = FALSE)

  res[[i]] <- data.frame(x = a$x,
                         y = a$y,
                         Species= as.character(a$marks),
                         alpha = params$alpha[i],
                         var.consp = params$var.consp[i])
}
res <- do.call("rbind", res)


pp.pairs <-  ggplot(data = res, aes(x = x, y = y)) +
              geom_point(aes(colour = Species), shape = 19, alpha = 0.5) +
              facet_grid(var.consp ~ alpha, labeller = "label_both") +
              theme_bw() + theme(panel.grid = element_blank()) +
              theme(legend.position = "bottom",
                    legend.background = element_rect(fill = "lightgray")) +
             scale_x_continuous(breaks = c(0, 1)) +
             scale_y_continuous(breaks = c(0, 1))
pp.pairs

# ------------------------------------


res <- list()
for(i in 1:length(alpha.vec))
{
  d <- seq(0, 1, by = 0.01)
  a <- PDFtexp(d,
               alpha = alpha.vec[[i]],
               dlim =c(0, 1))
  res[[i]] <- data.frame(Distance = d, pdf = a, alpha = alpha.vec[i])
}
res <- do.call("rbind", res)

pdfs <- ggplot(data = res, aes(x = Distance, y = pdf)) +
        geom_line() +
        facet_grid(.~alpha, labeller = "label_both") + theme_bw() +
        theme(panel.grid = element_blank()) +
        scale_x_continuous(breaks = c(0, 1)) +
        ylab("Probability density")
pdfs

png("../Figures/ISA_vs_CSA_pairs_simulations.png",
    width=2100, height = 2000, res=250 )
grid.arrange(pdfs, pp.pairs, ncol=1, nrow=2, heights = c(0.4, 1.1))
dev.off()



# ------------------------------------------------------------------------------
# SIMULATIONS - METHOD COMPARISON

params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1),
                      alpha = seq(-20, 20, by=2.5),
                      #alpha = c(-100, -10, -1, 0, 1, 10, 100)
                      grain = c(64, 32, 16, 8, 4),
                      N1 = c(10, 100, 1000),
                      N2 = c(10, 100, 1000))

res <- list()
for(i in 1:nrow(params))
{
  message(paste("Simulation", i, "out of", nrow(params)))

  a <- sim.pair(abund.vect = c(params$N1[i], params$N2[i]),
                var.consp = params$var.consp[i],
                alpha = params$alpha[i],
                plot.comm = FALSE)

  n.cells.side <- sqrt(params[i, 'grain'])

  m <- spasm::ppp.to.comm(a, dim.yx = c(n.cells.side, n.cells.side))$abundance

  # convert to binary matrix
  m.bin <-  m
  m.bin[m >= 1] <- 1


  res[[i]] <- data.frame(C_segSc = mean( spasm::C_seg(m.bin, scale=TRUE)),
                         C_togSc = mean( spasm::C_tog(m.bin, scale=TRUE)),
                         C_jacc = mean( spasm::C_jacc(m.bin)),
                         C_sor = mean( spasm::C_sor(m.bin)),
                         C_forbes = mean( spasm::C_forbes(m.bin)),
                         C_pears = mean(spasm::C_pears(m.bin)),
                         CA_cov = mean(spasm::CA_cov_cor(m, correlation=FALSE)),
                         CA_cov_hell = mean(spasm::CA_cov_cor(m, correlation=FALSE, transf="hellinger")),
                         CA_cor = mean(spasm::CA_cov_cor(m, correlation=TRUE)),
                         CA_cor_hell = mean(spasm::CA_cov_cor(m, correlation=TRUE, transf = "hellinger")),
                         CA_hell = mean(spasm::CA_hell(m)),
                         CA_tau = mean(spasm::CA_cov_cor(m, correlation=TRUE, transf=NULL, method="kendall")),
                         CA_bray = mean(spasm::CA_bray(m)),
                         CA_ruz = mean(spasm::CA_ruz(m)),
                         CA_chi = mean(spasm::CA_chi(m)),
                              params[i,])
}

res.all <- do.call("rbind", res) %>%
           mutate(Area = 1/(grain^2))

#write.csv(res.all, file="../Data/sim_results.csv", row.names = FALSE)

#res.all <- read.csv("../Data/sim_results.csv")

# ------------------------------------------------------------------------------



C.measures <- names(res.all)[grep(x=names(res.all), pattern="C_")]
CA.measures <- names(res.all)[grep(x=names(res.all), pattern="CA_")]

library(tidyr)

res <- reshape2::melt(data = res.all, measure.vars = c(C.measures,
                                                       CA.measures))


res <- res[is.na(res$value) == FALSE,]
write.csv(res, file = "sim_results.csv", row.names=FALSE)
# ------------------------------------------------------------------------------

res <- read.csv("sim_results.csv")
cor.res <- plyr::ddply(.data=res,
                       .variables=c("variable"),
                       .fun=summarize,
                       Kendall = cor(alpha, value))
type <- as.character(cor.res$variable)
type[grep(type, pattern="C_")] <- "Incidence-based"
type[grep(type, pattern="CA_")] <- "Abundance-based"
cor.res <- data.frame(cor.res, type)


png("../Figures/simulations_scatterplots.png", width = 2000, height = 1600, res = 150)
ggplot(data=res, aes(x=alpha, y = value)) +
  geom_vline(xintercept= 0) +
  geom_point(shape = 1, alpha = 0.3) +
  #geom_point(shape = 1, aes(colour = var.consp)) +
  facet_wrap(.~variable, scales="free") +
  geom_smooth(se = FALSE, span = 0.4, colour = "red") +
  theme_bw()
dev.off()

png("../Figures/simulations_bars.png", width = 800, height = 500, res = 150)
ggplot(data=cor.res, aes(x= reorder(variable, abs(Kendall)), y = abs(Kendall))) +
  geom_col(aes(fill = type)) + ylim(c(0, 1)) + theme_bw() + coord_flip() +
  ylab("abs (Kendall correlation with truth)") +
  xlab("")
dev.off()

# ------------------------------------------------------------------------------

# SIMULATIONS - SPATIALLY EXPLICIT METHODS


alphas <- c(-20, -10, 0, 10, 20)
distance <- seq(0,1, by=0.001)

res.f <- list()
res.pp <- list()

for(i in 1:length(alphas))
{
  alpha = alphas[i]

  pdf <- PDFtexp(distance,
                 alpha = alpha,
                 dlim =c(0, 1))
  #plot(distance, pdf, type = "l")

  set.seed(12345)
  a <- sim.pair(abund.vect = c(200, 200),
                var.consp = 0.01,
                alpha = alpha,
                plot.comm = FALSE)

  pcf <- pcfcross(X = a,
                  i = "sp1",
                  j = "sp2",
                  correction="translate",
                  r = distance)


  res.f[[i]] <- data.frame(distance, pdf, g = pcf$trans, alpha)
  res.pp[[i]]<- data.frame(x = a$x,
                           y = a$y,
                           Species= as.character(a$marks),
                           alpha)
  #plot(pcf, xlim=c(0,1), ylim = c(0,max(pdf)))
  #lines(distance, pdf)
}
res.f <- do.call(rbind, res.f)
res.f <- res.f[res.f$g != Inf,]
res.pp <- do.call(rbind, res.pp)

pdf <- ggplot(data=res.f, aes(x=distance, y=pdf)) +
  geom_line() + facet_grid(.~alpha)
pdf

g <- ggplot(data=res.f, aes(x=distance, y=g)) +
  facet_grid(.~alpha, labeller = "label_both") +
  geom_line(aes(x = distance, y = pdf), colour = "black") +
  geom_hline(yintercept=1, linetype = 2, colour="grey") +
  theme_bw() + theme(panel.grid = element_blank()) +
  geom_line(colour = "red")
g

pp.pairs <-  ggplot(data = res.pp, aes(x = x, y = y)) +
  geom_point(aes(colour = Species), shape = 19, alpha = 0.5) +
  facet_grid(. ~ alpha, labeller = "label_both") +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "lightgray")) +
  scale_x_continuous(breaks = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 1))
pp.pairs

png("../Figures/ISA_vs_CSA_PCF.png",
    width=2100, height = 1300, res=250 )
grid.arrange(g, pp.pairs, ncol=1, nrow=2, heights = c(1, 1.2))
dev.off()

# ------------------------------------------------------------------------------
# SIMULATIONS - VARIOGRAMS

library(vario)
library(devtools)
#remove.packages("vario", lib="~/R/x86_64-pc-linux-gnu-library/3.5")
#install('/home/pk33loci/Dropbox/Interspecific_aggregation/vario')
library.dynam(chname="vario", package = "vario", lib.loc=.libPaths())


comm <- sim.pair(abund.vect = c(200, 200),
              var.consp = 0.1,
              alpha = -10,
              plot.comm = TRUE)

comm <- ppp.to.comm(comm, dim.yx = c(10,10))
abu <- t(comm$abundance)#[,1:2]
xy <- comm$coords

v.obs <- vario(x = as.matrix(abu),
               coord = xy,
               grain = 0.1,
               pos.neg = TRUE)
ggplot(data=v.obs$vario, aes(x=Dist, y = pos + neg)) +
   geom_line()




