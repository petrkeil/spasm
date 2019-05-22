library(ggplot2)
library(dplyr)
library(vegan)
library(tidyr)
library(gridExtra)
library(raster)
library(spasm)
library(spatstat)
library(ggrepel)

# load Aniko Toth's functions from GitHub
source("https://raw.githubusercontent.com/anikobtoth/FCW/master/Pair_Functions.R")

################################################################################

params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1),
                      alpha = seq(-20, 20, by=2.5),
                      #alpha = c(-100, -10, -1, 0, 1, 10, 100)
                      grain = c(64, 32, 16, 8, 4),
                      N1 = c(10, 100, 1000, 10000),
                      N2 = c(10, 100, 1000, 10000))

#params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1),
#                      alpha = seq(-20, 20, by=5),
#                      grain = c(32, 16, 8, 4),
#                      N1 = c(200, 1000),
#                      N2 = c(200, 1000))


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


  res[[i]] <- data.frame(
                         C_segSc = mean( spasm::C_seg(m.bin)),
                         C_segSc_Z = mean(spasm::Z_score(m.bin, "step_C_sim2", "C_seg", N.sim=100)),
                         C_togSc = mean( spasm::C_tog(m.bin)),
                         C_togSc_Z = mean(spasm::Z_score(m.bin, "step_C_sim2", "C_tog", N.sim=100)),
                         C_jacc = mean( spasm::C_jacc(m.bin)),
                         C_jacc_Z = mean(spasm::Z_score(m.bin, "step_C_sim2", "C_jacc", N.sim=100)),
                         C_sor = mean( spasm::C_sor(m.bin)),
                         C_sor_Z = mean(spasm::Z_score(m.bin, "step_C_sim2", "C_sor", N.sim=100)),
                         C_forbes = mean( spasm::C_forbes(m.bin)),
                         C_FETmP = mean(FETmP_Pairwise(m.bin)),
                         #C_FET = mean(FET_Pairwise(m.bin)),
                         C_pears = mean(spasm::C_pears(m.bin)),
                         C_match = mean(spasm::C_match(m.bin)),
                         C_match_Z = mean(spasm::Z_score(m.bin, "step_C_sim2", "C_match", N.sim=100)),
                         CA_cov = mean(spasm::CA_cov_cor(m, correlation=FALSE)),
                         CA_cov_hell = mean(spasm::CA_cov_cor(m, correlation=FALSE, transf="hellinger")),
                         CA_cor = mean(spasm::CA_cov_cor(m, correlation=TRUE, method="pearson")),
                         CA_cor_hell = mean(spasm::CA_cov_cor(m, correlation=TRUE, transf = "hellinger", method="pearson")),
                         CA_hell = mean(spasm::CA_hell(m)),
                         CA_hell_Z = mean(spasm::Z_score(m, "step_CA_rowrandom", "CA_hell", N.sim=100)),
                         CA_tau = mean(spasm::CA_cov_cor(m, correlation=TRUE, method="kendall")),
                         CA_tau_Z = mean(spasm::Z_score(m, "step_CA_rowrandom", "CA_cov_cor", N.sim=100)),
                         CA_bray = mean(spasm::CA_bray(m)),
                         CA_bray_Z = mean(spasm::Z_score(m, "step_CA_rowrandom", "CA_bray", N.sim=100)),
                         CA_ruz = mean(spasm::CA_ruz(m)),
                         CA_ruz_Z = mean(spasm::Z_score(m, "step_CA_rowrandom", "CA_ruz", N.sim=100)),
                         CA_chi = mean(spasm::CA_chi(m)),
                         CA_chi_Z = mean(spasm::Z_score(m, "step_CA_rowrandom", "CA_chi", N.sim=100)),
                         params[i,]
                        )
}

res.all <- do.call("rbind", res) %>%
  mutate(Area = 1/(grain^2))

write.csv(res.all, file="sim_results.csv", row.names = FALSE)
# ------------------------------------------------------------------------------

res.all <- read.csv("sim_results.csv")
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


png("figures/simulations_bars.png", width = 1600, height = 900, res = 150)
  grid.arrange(bar.all, bar.neg, bar.pos, nrow=1, ncol=3)
dev.off()



################################################################################

# -------------------------summarize by scale
scale <- plyr::ddply(.data=res,
                       .variables=c("variable", "grain"),
                       .fun=summarize,
                       Kendall = abs(cor(alpha, value, method = "kendall")))
scale$variable <- as.factor(scale$variable)

# abundance or incidence?
type <- as.character(scale$variable)
type[grep(type, pattern="C_")] <- "Incidence-based"
type[grep(type, pattern="CA_")] <- "Abundance-based"
scale <- data.frame(scale, type)

# raw or Z-score?
Z <- rep("Raw", times = nrow(scale))
Z[grep(scale$variable, pattern="_Z")] <- "Z-score"
scale <- data.frame(scale, Z)

# colour palette
cols <- data.frame(variable = cor.res$variable,
                   rnd.col = as.factor(sample(1:nrow(cor.res))) )
scale <- left_join(x=scale, y = cols, by="variable")

# label positions
scale.labs <- scale[scale$grain == 32,]

gr <- ggplot(data = scale, aes(x = grain, y = Kendall)) +
  geom_point(aes(colour = rnd.col)) +
  geom_line(aes(colour = rnd.col)) +
  scale_x_continuous(trans = "log2") +
  theme_bw() +
  facet_grid(Z~type) +
  geom_label_repel(data = scale.labs,
                   aes(x = grain, y = Kendall, label=variable,
                       colour = rnd.col)) +
  theme(legend.position="none") +
  xlab("Grain [# pixels along one side]") +
  ylab("| Kendall correlation with ISA |")
gr

png("figures/simulations_kendall_by_grain.png",
    width = 1400, height= 1400, res=200)
 gr
dev.off()


# -------------------------summarize by CSA

csa <- plyr::ddply(.data=res,
                       .variables=c("variable", "var.consp"),
                       .fun=summarize,
                       Kendall = abs(cor(alpha, value, method = "kendall")))
csa$variable <- as.factor(csa$variable)
names(csa)[2] <- "CSA"

# abundance or incidence?
type <- as.character(csa$variable)
type[grep(type, pattern="C_")] <- "Incidence-based"
type[grep(type, pattern="CA_")] <- "Abundance-based"
csa <- data.frame(csa, type)

# raw or Z-score?
Z <- rep("Raw", times = nrow(csa))
Z[grep(csa$variable, pattern="_Z")] <- "Z-score"
csa <- data.frame(csa, Z)

# colour palette
csa <- left_join(x=csa, y = cols, by="variable")

# label positions
scale.labs <- csa[csa$CSA == 0.01,]

cs <- ggplot(data = csa, aes(x = CSA, y = Kendall)) +
  geom_point(aes(colour = rnd.col)) +
  geom_line(aes(colour = rnd.col)) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  facet_grid(Z~type) +
  geom_label_repel(data = scale.labs,
                   aes(x = CSA, y = Kendall, label=variable,
                       colour = rnd.col)) +
  theme(legend.position="none") +
  xlab("Con-specific aggregation (CSA)") +
  ylab("| Kendall correlation with ISA |")
cs

png("figures/simulations_kendall_by_CSA.png",
    width = 1400, height= 1400, res=200)
cs
dev.off()
