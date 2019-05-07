library(ggplot2)
library(dplyr)
library(vegan)
library(tidyr)
library(gridExtra)
library(raster)
library(spasm)
library(spatstat)

# load Aniko Toth's functions from GitHub
source("https://raw.githubusercontent.com/anikobtoth/FCW/master/Pair_Functions.R")

################################################################################

params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1),
                      alpha = seq(-20, 20, by=2.5),
                      #alpha = c(-100, -10, -1, 0, 1, 10, 100)
                      grain = c(128, 64, 32, 16, 8, 4),
                      N1 = c(10, 100, 1000, 10000),
                      N2 = c(10, 100, 1000, 10000))

params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1),
                      alpha = seq(-20, 20, by=5),
                      #alpha = c(-100, -10, -1, 0, 1, 10, 100)
                      grain = c(64, 32, 16, 8, 4),
                      N1 = c(200, 1000),
                      N2 = c(200, 1000))


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
                         C_sorZ = spasm::C_sorZ(m.bin),
                         C_forbes = mean( spasm::C_forbes(m.bin)),
                         C_FETmP = mean(FETmP_Pairwise(m.bin)),
                         C_FET = mean(FET_Pairwise(m.bin)),
                         C_pears = mean(spasm::C_pears(m.bin)),
                         CA_cov = mean(spasm::CA_cov_cor(m, correlation=FALSE)),
                         CA_cov_hell = mean(spasm::CA_cov_cor(m, correlation=FALSE, transf="hellinger")),
                         CA_cor = mean(spasm::CA_cov_cor(m, correlation=TRUE)),
                         CA_cor_hell = mean(spasm::CA_cov_cor(m, correlation=TRUE, transf = "hellinger")),
                         CA_hell = mean(spasm::CA_hell(m)),
                         CA_tau = mean(spasm::CA_cov_cor(m, correlation=TRUE, transf=NULL, method="kendall")),
                         CA_tauZ = spasm::CA_tauZ(m),
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

# ------------------------------------------------------------------------------
#  write.csv(res, file = "sim_results.csv", row.names=FALSE)
# ------------------------------------------------------------------------------


chr <- as.character(res$variable)
#res$variable <- factor(x= res$variable, levels = unique(chr))
res$variable <- ordered(x= res$variable, levels = unique(chr))
res$sign <- rep("independence", times=nrow(res))
res$sign[res$alpha > 0] <- "attraction"
res$sign[res$alpha < 0] <- "repulsion"


png("figures/simulations_scatterplots.png",
    width = 2000, height = 1600, res = 200)
ggplot(data=res, aes(x=alpha, y = value)) +
  geom_vline(xintercept= 0) +
  geom_point(shape = 1, alpha = 0.3, colour = "darkgrey") +
  geom_smooth(se = FALSE, colour = "#FF7570") +
  theme_bw() + xlab("ISA") +
  facet_wrap(~variable, scales="free")
dev.off()

# ------------------------------------------------------------------------------
# overall correlation

cor.res <- plyr::ddply(.data=res,
                       .variables=c("variable"),
                       .fun=summarize,
                       Kendall = cor(alpha, value, method = "kendall"))
type <- as.character(cor.res$variable)
type[grep(type, pattern="C_")] <- "Incidence-based"
type[grep(type, pattern="CA_")] <- "Abundance-based"
cor.res <- data.frame(cor.res, type)

# ------------------------------------------------------------------------------
# by negative or positive ISA

cor.posneg <- plyr::ddply(.data=res,
                          .variables=c("variable", "sign"),
                          .fun=summarize,
                          Kendall = cor(alpha, value, method = "kendall"))
type <- as.character(cor.posneg$variable)
type[grep(type, pattern="C_")] <- "Incidence-based"
type[grep(type, pattern="CA_")] <- "Abundance-based"
cor.posneg <- data.frame(cor.posneg, type)

cor.pos <- cor.posneg[cor.posneg$sign == "attraction",]
cor.neg <- cor.posneg[cor.posneg$sign == "repulsion",]



# -----------------------------------------------------------------------------
# barplots
bar.all <- ggplot(data=cor.res,
                  aes(x= reorder(variable, abs(cor.res$Kendall)), y = abs(Kendall))) +
            geom_col(aes(fill = type)) +
            ylim(c(0, 1)) + theme_bw() + coord_flip() +
            ylab("| Kendall correlation with overall ISA |") +
            theme(legend.position="none") +
            xlab("") + ggtitle("(a)")

bar.neg <- ggplot(data=cor.neg, aes(x= reorder(variable, abs(cor.res$Kendall)), y = abs(Kendall))) +
  geom_col(aes(fill = type)) + ylim(c(0, 1)) +
  theme_bw() + coord_flip() +
  ylab("| Kendall correlation with repulsion |") +
  theme(legend.position="none") +
  xlab("")+ ggtitle("(b)")

bar.pos <- ggplot(data=cor.pos, aes(x= reorder(variable, abs(cor.res$Kendall)), y = abs(Kendall))) +
  geom_col(aes(fill = type)) + ylim(c(0, 1)) +
  theme_bw() + coord_flip() +
  ylab("| Kendall correlation with attraction |") +
  theme(legend.position=c(0.7, 0.3)) +
  xlab("") + ggtitle("(c)")


png("figures/simulations_bars.png", width = 1600, height = 600, res = 150)
  grid.arrange(bar.all, bar.neg, bar.pos, nrow=1, ncol=3)
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
  xlab("CSA") + ylab("| Kendall correlation with overall ISA |")

png("figures/simulations_kendall_by_params.png",
    width = 1400, height= 600, res=130)
grid.arrange(gr, csa, nrow=1, ncol=2, widths = c(1, 1.3))
dev.off()




