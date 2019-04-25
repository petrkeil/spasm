library(ggplot2)
library(dplyr)
library(vegan)
library(tidyr)
library(gridExtra)
library(raster)
library(spasm)
library(spatstat)

# ------------------------------------------------------------------------------
# SIMULATIONS - ILLUSTRATIONS

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




# ------------------------------------


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
        ylab("pdf") + xlab("Distance from sp1") +
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
# SIMULATIONS - METHOD COMPARISON

params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1),
                      alpha = seq(-20, 20, by=2.5),
                      #alpha = c(-100, -10, -1, 0, 1, 10, 100)
                      grain = c(128, 64, 32, 16, 8, 4),
                      N1 = c(10, 100, 1000, 10000),
                      N2 = c(10, 100, 1000, 10000))

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
# ------------------------------------------------------------------------------

#res.all <- read.csv("../Data/sim_results.csv")
# ------------------------------------------------------------------------------

C.measures <- names(res.all)[grep(x=names(res.all), pattern="C_")]
CA.measures <- names(res.all)[grep(x=names(res.all), pattern="CA_")]

library(tidyr)

res <- reshape2::melt(data = res.all, measure.vars = c(C.measures,
                                                       CA.measures))

res$variable <- as.character(res$variable)
res <- res[is.na(res$value) == FALSE,]
write.csv(res, file = "sim_results.csv", row.names=FALSE)
# ------------------------------------------------------------------------------

res <- read.csv("sim_results.csv")
chr <- as.character(res$variable)
#res$variable <- factor(x= res$variable, levels = unique(chr))
res$variable <- ordered(x= res$variable, levels = unique(chr))


cor.res <- plyr::ddply(.data=res,
                       .variables=c("variable"),
                       .fun=summarize,
                       Kendall = cor(alpha, value, method = "kendall"))
type <- as.character(cor.res$variable)
type[grep(type, pattern="C_")] <- "Incidence-based"
type[grep(type, pattern="CA_")] <- "Abundance-based"
cor.res <- data.frame(cor.res, type)

png("figures/simulations_scatterplots.png",
    width = 2000, height = 1600, res = 200)
ggplot(data=res, aes(x=alpha, y = value)) +
  geom_vline(xintercept= 0) +
  geom_point(shape = 1, alpha = 0.3, colour = "darkgrey") +
  geom_smooth(se = FALSE, colour = "#FF7570") +
  theme_bw() + xlab("ISA") +
  facet_wrap(~variable, scales="free")
dev.off()

png("figures/simulations_bars.png", width = 800, height = 500, res = 150)
ggplot(data=cor.res, aes(x= reorder(variable, abs(Kendall)), y = abs(Kendall))) +
  geom_col(aes(fill = type)) + ylim(c(0, 1)) + theme_bw() + coord_flip() +
  ylab("| Kendall correlation with ISA |") +
  xlab("")
dev.off()

# -------------------------summarize by scale
cor.res <- plyr::ddply(.data=res,
                       .variables=c("variable", "grain"),
                       .fun=summarize,
                       Kendall = abs(cor(alpha, value, method = "kendall")))
cor.res$variable <- as.factor(cor.res$variable)
#cor.res$grain <- (1/cor.res$grain)^2

gr <- ggplot(data = cor.res, aes(x = grain, y = Kendall)) +
             geom_point(aes(colour = variable)) +
             geom_line(aes(colour = variable)) +
             scale_x_continuous(trans = "log2") +
             theme_bw() + ggtitle("(a)") +
             theme(legend.position="none") +
             xlab("Grain [# pixels along one side]") + ylab("| Kendall correlation with ISA |")

cor.res <- plyr::ddply(.data=res,
                       .variables=c("variable", "var.consp"),
                       .fun=summarize,
                       Kendall = abs(cor(alpha, value, method = "kendall")))
cor.res$variable <- as.factor(cor.res$variable)

csa <-  ggplot(data = cor.res, aes(x = var.consp, y = Kendall)) +
               geom_point(aes(colour = variable)) +
               geom_line(aes(colour = variable)) +
               scale_x_continuous(trans = "log10") +
               theme_bw() + ggtitle("(b)") +
               xlab("CSA") + ylab("| Kendall correlation with ISA |")

png("figures/simulations_kendall_by_params.png",
    width = 1400, height= 600, res=130)
grid.arrange(gr, csa, nrow=1, ncol=2, widths = c(1, 1.3))
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
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_line(colour = "red") +
  scale_x_continuous(breaks = c(0, 1)) +
  labs(title = "(b)", x = "Distance from sp1", y = "pdf (black) or pcf (red)")
g

pp.pairs <-  ggplot(data = res.pp, aes(x = x, y = y)) +
  geom_point(aes(colour = Species), shape = 19, alpha = 0.5) +
  facet_grid(. ~ alpha, labeller = "label_both") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.9, 0.1),
        #legend.position = "none",
        legend.background = element_rect(fill = "lightgray"),
        plot.margin=unit(c(5.5,5.5,5.5,10), "points")) +
  scale_x_continuous(breaks = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 1)) +
  labs(title = "(a)")
pp.pairs

png("figures/simulations_pcf.png",
    width=2200, height = 1300, res=250 )
grid.arrange(pp.pairs, g, ncol=1, nrow=2, heights = c(1, 1))
grid::grid.text("repulsion, seggregation", x = 0.23, y = 0.97)
grid::grid.text("independence", x = 0.52, y = 0.97)
grid::grid.text("attraction", x = 0.8, y = 0.97)

dev.off()

# ------------------------------------------------------------------------------
# SIMULATIONS - VARIOGRAMS

library(vario)
library(devtools)
#remove.packages("vario", lib="~/R/x86_64-pc-linux-gnu-library/3.5")
#install('/home/pk33loci/Dropbox/Interspecific_aggregation/vario')
library.dynam(chname="vario", package = "vario", lib.loc=.libPaths())


comm <- sim.pair(abund.vect = c(200, 200),
              var.consp = 0.01,
              alpha = 20,
              plot.comm = FALSE)

comm <- ppp.to.comm(comm, dim.yx = c(20,20))
abu <- t(comm$abundance)#[,1:2]
xy <- comm$coords

v.obs <- vario(x = as.matrix(abu),
               coord = xy,
               grain = 0.05,
               pos.neg = TRUE,
               hmax = 1,
               binary = FALSE)

neg <- data.frame(value = v.obs$vario$neg, sign = "Negative covariance")
pos <- data.frame(value = v.obs$vario$pos, sign = "Positive covariance")
posneg <- rbind(neg, pos)


ggplot(data=v.obs$vario, aes(x=Dist, y = -1*(neg))) +
   geom_line(colour = "red") +
   geom_line(aes(x=Dist, y = pos), colour = "green") +
  xlim(c(0,1))

# ------------------------------------------------------------------------------
params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1),
                      alphas = c(-20, -10, 0, 10, 20))

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
                 binary = TRUE)

  pos <- data.frame(Dist = v.obs$vario$Dist,
                    Covariance = v.obs$vario$pos,
                    Sign = "Positive")
  neg <- data.frame(Dist = v.obs$vario$Dist,
                    Covariance = abs(v.obs$vario$neg),
                    Sign = "| Negative |")
  posneg <- rbind(neg, pos)

  res.vario[[i]] <- data.frame(posneg,
                               ISA = params$alphas[i],
                               CSA = params$var.consp[i])
}

res.vario <- do.call(rbind, res.vario)

ggplot(data=res.vario, aes(x=Dist, Covariance)) +
  geom_line(aes(colour = Sign)) +
  geom_point(aes(colour = Sign)) +
  facet_grid(CSA~ISA, scales = "free_y", labeller = "label_both") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        #legend.position = "none",
        legend.background = element_rect(fill = "lightgray"),
        plot.margin=unit(c(5.5,5.5,5.5,10), "points")) +
  scale_x_continuous(breaks = c(0, 1), limits=c(0,1)) +
  geom_hline(yintercept=0)



