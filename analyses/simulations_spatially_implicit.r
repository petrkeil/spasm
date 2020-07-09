################################################################################
#
# Code for the simulations behind spatially implicit simulations and the comparison
# of the spatially implicit metrics of ISA
#
#
# Petr Keil
#
################################################################################

library(ggplot2)
library(dplyr)
library(vegan)
library(tidyr)
library(gridExtra)
library(raster)
library(spasm)
library(spatstat)
library(ggrepel)

# load Aniko Toth's functions from GitHub (for C_FETmP)
source("https://raw.githubusercontent.com/anikobtoth/FCW/master/Pair_Functions.R")

################################################################################


# Set simulation parameters
params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1),  # CSA
                      alpha = seq(-20, 20, by=2.5),       # ISA
                      grain = c(32, 16, 8, 4),        # number of grid cells along one side
                      N1 = c(10, 100, 1000, 10000),       # number of individuals of species 1
                      N2 = c(10, 100, 1000, 10000))       # number of individuals of species 2


# ------------------------------------------------------------------------------
# The grand loop that does all the simulations
# ------------------------------------------------------------------------------

set.seed(12345)

res <- list() # empy container for results

for(i in 1:nrow(params)) # for each combination of parameters
{
  message(paste("Simulation", i, "out of", nrow(params)))

  # simulate point patterns of the two species
  a <- sim.pair(abund.vect = c(params$N1[i], params$N2[i]),
                var.consp = params$var.consp[i],
                alpha = params$alpha[i],
                plot.comm = FALSE)

  n.cells.side <- params[i, 'grain']

  m <- spasm::ppp.to.comm(a, dim.yx = c(n.cells.side, n.cells.side))$abundance

  # create a binary matrix as well
  m.bin <-  m
  m.bin[m >= 1] <- 1

  # which Z-score randomization algorithm to use for abundance data?
  algo.abu <- "step_CA_IT"

  # how many null model simulations?
  N = 200

  # all of the metrics
  res[[i]] <- data.frame(
                         C_segSc = mean( spasm::C_seg(m.bin)),
                         C_segSc_Z = mean(spasm::Z_score(m.bin, "step_C_sim2", "C_seg", N.sim=N)),
                         C_togSc = mean( spasm::C_tog(m.bin)),
                         C_togSc_Z = mean(spasm::Z_score(m.bin, "step_C_sim2", "C_tog", N.sim=N)),
                         C_jacc = mean( spasm::C_jacc(m.bin)),
                         C_jacc_Z = mean(spasm::Z_score(m.bin, "step_C_sim2", "C_jacc", N.sim=N)),
                         C_sor = mean( spasm::C_sor(m.bin)),
                         C_sor_Z = mean(spasm::Z_score(m.bin, "step_C_sim2", "C_sor", N.sim=N)),
                         C_forbes = mean( spasm::C_forbes(m.bin)),
                         C_alroy = mean(spasm::C_alroy(m.bin)),
                         C_FETmP = mean(FETmP_Pairwise(m.bin)),
                         C_pears = mean(spasm::C_pears(m.bin)),
                         C_match = mean(spasm::C_match(m.bin)),
                         C_match_Z = mean(spasm::Z_score(m.bin, "step_C_sim2", "C_match", N.sim=N)),
                         CA_cov = mean(spasm::CA_cov_cor(m, correlation=FALSE)),
                         CA_cov_hell = mean(spasm::CA_cov_cor(m, correlation=FALSE, transf="hellinger")),
                         CA_cor = mean(spasm::CA_cov_cor(m, correlation=TRUE, method="pearson")),
                         CA_cor_hell = mean(spasm::CA_cov_cor(m,
                                                              correlation=TRUE,
                                                              transf = "hellinger",
                                                              method="pearson")),
                         CA_cor_hell_Z = mean(spasm::Z_score(m, "step_CA_IT", "CA_cov_cor",
                                                             correlation = TRUE,
                                                             method = "pearson",
                                                             transf = "hellinger",
                                                             N.sim=N)),
                         CA_hell = mean(spasm::CA_hell(m)),
                         CA_hell_Z = mean(spasm::Z_score(m, "step_CA_IT", "CA_hell", N.sim=N)),
                         CA_rho = mean(spasm::CA_cov_cor(m, correlation=TRUE, method="spearman")),
                         CA_rho_Z = mean(spasm::Z_score(m, "step_CA_IT", "CA_cov_cor",
                                                        method="spearman", N.sim=N)),
                         CA_bray = mean(spasm::CA_bray(m)),
                         CA_bray_Z = mean(spasm::Z_score(m, "step_CA_IT", "CA_bray", N.sim=N)),
                         CA_ruz = mean(spasm::CA_ruz(m)),
                         CA_ruz_Z = mean(spasm::Z_score(m, "step_CA_IT", "CA_ruz", N.sim=N)),
                         CA_chi = mean(spasm::CA_chi(m)),
                         CA_chi_Z = mean(spasm::Z_score(m, "step_CA_IT", "CA_chi", N.sim=N)),
                         params[i,] # saving the simulation parameters as well
                        )
}

# convert the list to a data frame
res.all <- do.call("rbind", res) %>%
  mutate(Area = 1/(grain^2))

# save the results to a file (so that I don't have to re-run the loop all the time)
write.csv(res.all, file="sim_results.csv", row.names = FALSE)
# ------------------------------------------------------------------------------

# read the results from the file
res.all <- read.csv("sim_results.csv")
# ------------------------------------------------------------------------------

# classify abundance- and incidence-based measures
C.measures <- names(res.all)[grep(x=names(res.all), pattern="C_")]
CA.measures <- names(res.all)[grep(x=names(res.all), pattern="CA_")]
res <- reshape2::melt(data = res.all, measure.vars = c(C.measures,
                                                       CA.measures))

# remove outliers, NAs, infinity values and similar mess
res$variable <- as.character(res$variable)
res <- res[is.na(res$value) == FALSE,] # remove NAs
res <- res[res$value != Inf,]    # remove infinity
res <- res[res$value != (-Inf),] # remove negative infinity
#res <- res[res$value < 1000,]    # remove extreme values - not used for the paper
#res <- res[res$value > -1000,]   # remove extreme values - not used for the paper

# some other transfomrations and additions to the data
chr <- as.character(res$variable)
res$variable <- ordered(x= res$variable, levels = unique(chr))
res$sign <- rep("independence", times=nrow(res))
res$sign[res$alpha > 0] <- "attraction"
res$sign[res$alpha < 0] <- "repulsion"

################################################################################
# FIGURES - SCATTERPLOTS
################################################################################

# figure with all the details on the relationships
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
# overall correlation between the metric and the simulated ISA (represented by alpha)
cor.res <- plyr::ddply(.data=res,
                       .variables=c("variable"),
                       .fun=summarize,
                       Spearman = cor(alpha, value, method = "spearman"))

# some further data transformations and classifications
data <- as.character(cor.res$variable)
data[grep(data, pattern="C_")] <- "Binary"
data[grep(data, pattern="CA_")] <- "Abundance"

type <- as.character(cor.res$variable)
type[] <- "Raw"
type[grep(as.character(cor.res$variable), pattern="_Z")] <- "Z-score"

cor.res <- data.frame(cor.res, data, type)

# ------------------------------------------------------------------------------
# sorting the results by negative or positive ISA

cor.posneg <- plyr::ddply(.data=res,
                          .variables=c("variable", "sign"),
                          .fun=summarize,
                          Spearman = cor(alpha, value, method = "spearman"))
data <- as.character(cor.posneg$variable)
data[grep(data, pattern="C_")] <- "Binary"
data[grep(data, pattern="CA_")] <- "Abundance"

type <- as.character(cor.posneg$variable)
type[] <- "Raw"
type[grep(as.character(cor.posneg$variable), pattern="_Z")] <- "Z-score"

cor.posneg <- data.frame(cor.posneg, data, type)

cor.pos <- cor.posneg[cor.posneg$sign == "attraction",]
cor.neg <- cor.posneg[cor.posneg$sign == "repulsion",]


cor.res <- cor.res[cor.res$type == "Raw", ]
cor.pos <- cor.pos[cor.pos$type == "Raw", ]
cor.neg <- cor.neg[cor.neg$type == "Raw", ]

################################################################################
# FIGURES - BARPLOT WITH ONLY RAW METRICS
################################################################################

# barplot - overall correlation
bar.all <- ggplot(data=cor.res,
                  aes(x= reorder(variable, abs(cor.res$Spearman)), y = abs(Spearman))) + # note the abs
  geom_col(aes(fill = data)) +
  ylim(c(0, 1)) + theme_bw() + coord_flip() +
  ylab("| Correlation with overall ISA |") +
  theme(legend.position="none") +
  # facet_grid(.~type) +
  xlab("") + ggtitle("(a)")

# barplot - negative correlation
bar.neg <- ggplot(data=cor.neg, aes(x= reorder(variable, abs(cor.res$Spearman)), y = abs(Spearman))) +
  geom_col(aes(fill = data)) + ylim(c(0, 1)) +
  theme_bw() + coord_flip() +
  ylab("| Correlation with repulsion |") +
  theme(legend.position="none") +
  xlab("")+ ggtitle("(b)")

# barplot - positive correlation
bar.pos <- ggplot(data=cor.pos, aes(x= reorder(variable, abs(cor.res$Spearman)), y = abs(Spearman))) +
  geom_col(aes(fill = data)) + ylim(c(0, 1)) +
  theme_bw() + coord_flip() +
  ylab("| Correlation with attraction |") +
  theme(legend.position=c(0.7, 0.3)) +
  xlab("") + ggtitle("(c)")

# export the figure
png("figures/simulations_bars_raw.png", width = 1600, height = 900, res = 150)
  grid.arrange(bar.all, bar.neg, bar.pos, nrow=1, ncol=3)
dev.off()

pdf("figures/simulations_bars_raw.pdf", width = 10, height = 5)
grid.arrange(bar.all, bar.neg, bar.pos, nrow=1, ncol=3)
dev.off()



################################################################################
# FIGURES - BARPLOT WITH RAW METRICS AND Z-SCORES TOGETHER
################################################################################

bar.all <- ggplot(data=cor.res,
                  aes(x= reorder(variable, abs(cor.res$Spearman)), y = abs(Spearman))) +
            geom_col(aes(fill = data)) +
            ylim(c(0, 1)) + theme_bw() + coord_flip() +
            ylab("| Correlation with overall ISA |") +
            theme(legend.position="none") +
           # facet_grid(.~type) +
            xlab("") + ggtitle("(a)")

bar.neg <- ggplot(data=cor.neg, aes(x= reorder(variable, abs(cor.res$Spearman)), y = abs(Spearman))) +
  geom_col(aes(fill = data)) + ylim(c(0, 1)) +
  theme_bw() + coord_flip() +
  ylab("| Correlation with repulsion |") +
  theme(legend.position="none") +
  xlab("")+ ggtitle("(b)")

bar.pos <- ggplot(data=cor.pos, aes(x= reorder(variable, abs(cor.res$Spearman)), y = abs(Spearman))) +
  geom_col(aes(fill = data)) + ylim(c(0, 1)) +
  theme_bw() + coord_flip() +
  ylab("| Correlation with attraction |") +
  theme(legend.position=c(0.7, 0.3)) +
  xlab("") + ggtitle("(c)")


pdf("figures/simulations_bars.pdf", width = 10, height = 5)
  grid.arrange(bar.all, bar.neg, bar.pos, nrow=1, ncol=3)
dev.off()





################################################################################
# FIGURES - MEAN ISA AS A FUNCTION OF GRAIN
################################################################################

# -------------------------summarize by scale
scale <- plyr::ddply(.data=res,
                     .variables=c("variable", "grain"),
                     .fun=summarize,
                     Spearman = abs(cor(alpha, value, method = "spearman")))
scale$variable <- as.character(scale$variable)

# abundance or incidence?
type <- as.character(scale$variable)
type[grep(type, pattern="C_")] <- "Binary data"
type[grep(type, pattern="CA_")] <- "Abundance data"
scale <- data.frame(scale, type)

# raw or Z-score?
Z <- rep("Raw", times = nrow(scale))
Z[grep(scale$variable, pattern="_Z")] <- "Z-score"
scale <- data.frame(scale, Z)

# colour palette
var.names <- as.character(unique(scale$variable))
cols <- data.frame(variable = as.character(var.names),
                   rnd.col = as.factor(sample(1:length(var.names))))
scale <- left_join(x=scale, y = cols, by="variable")

# label positions
scale.labs <- scale[scale$grain == 32,]

gr <- ggplot(data = scale, aes(x = grain, y = Spearman)) +
  geom_point(aes(colour = rnd.col)) +
  geom_line(aes(colour = rnd.col)) +
  scale_x_continuous(trans = "log2") +
  theme_bw() +
  facet_grid(Z~type) +
  geom_label_repel(data = scale.labs,
                   aes(x = grain, y = Spearman, label=variable,
                       colour = rnd.col)) +
  theme(legend.position="none") +
  xlab("Grain [# pixels along one side]") +
  ylab("Mean metric")
gr

png("figures/simulations_mean_by_GRAIN.png", width = 1400, height= 1400, res=200)
print(gr)
dev.off()

################################################################################
# FIGURES - CORRELATIONS SUMMARIZED BY CSA
################################################################################

csa <- plyr::ddply(.data=res,
                   .variables=c("variable", "var.consp"),
                   .fun=summarize,
                   Spearman = abs(cor(alpha, value, method = "spearman")))
csa$variable <- as.character(csa$variable)
names(csa)[2] <- "CSA"

# abundance or incidence?
type <- as.character(csa$variable)
type[grep(type, pattern="C_")] <- "Binary data"
type[grep(type, pattern="CA_")] <- "Abundance data"
csa <- data.frame(csa, type)

# colour palette
var.names <- as.character(unique(csa$variable))
cols <- data.frame(variable = as.character(var.names),
                   rnd.col = as.factor(sample(1:length(var.names))))
csa <- left_join(x=csa, y = cols, by="variable")

# raw or Z-score?
Z <- rep("Raw", times = nrow(csa))
Z[grep(csa$variable, pattern="_Z")] <- "Z-score"
csa <- data.frame(csa, Z)

# label positions
scale.labs <- csa[csa$CSA == 0.01,]

cs <- ggplot(data = csa, aes(x = CSA, y = Spearman)) +
  geom_point(aes(colour = rnd.col)) +
  geom_line(aes(colour = rnd.col)) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  facet_grid(Z~type) +
  geom_label_repel(data = scale.labs,
                   aes(x = CSA, y = Spearman, label=variable,
                       colour = rnd.col)) +
  theme(legend.position="none") +
  xlab("Con-specific aggregation (CSA)") +
  ylab("| Correlation with ISA |")
cs

png("figures/simulations_spearman_by_CSA.png", width = 1400, height= 1400, res=200)
print(cs)
dev.off()
