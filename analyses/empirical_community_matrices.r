################################################################################
#
# Analysis of empirical community matrices behind Figure 5
#
# Petr Keil
#
################################################################################

library(psych)
library(tidyverse)
library(gridExtra)
library(spasm)
library(qgraph)
library(factoextra)

# load Aniko Toth's functions from GitHub
source("https://raw.githubusercontent.com/anikobtoth/FCW/master/Pair_Functions.R")

# In case spasm isn't installed:
# library(devtools)
# install_github(repo="petrkeil/spasm")

# ------------------------------------------------------------------------------

# function removing -Inf, Inf, NaN and NA from a distance matrix
cleaner <- function(D)
{
  D <- D[]
  D <- na.omit(D)
  D <- D[D != Inf]
  D <- D[D != -Inf]
  return(D)
}

################################################################################
# INCIDENCE-BASED MEASURES
################################################################################

# prepare the data - they come with the package "spasm"
dat <- c(data.Ulrich, data.Atmar)

N = 200 # how many null model simulations?
set.seed(12345) # set seed for reproducibility

res <- list()
for(i in 1:length(dat))
{
  message(i)
  m <- dat[[i]]
  m.bin <- m
  m.bin[m.bin > 0] <- 1 # convert to binary matrix

  # calculate the metrics
  res.i <- c(N = sum(colSums(m.bin) > 0),
             Gamma = sum(rowSums(m.bin) > 0),
             Tot.incid = sum(m.bin),
             C_segSc = mean( spasm::C_seg(m.bin)),
             C_togSc = mean(spasm::C_tog(m.bin)),
             C_jacc = mean( spasm::C_jacc(m.bin)),
             C_sor = mean( spasm::C_sor(m.bin)),
             C_forbes = mean( spasm::C_forbes(m.bin)),
             C_alroy = mean(spasm::C_alroy(m.bin)),
             C_FETmP = mean(FETmP_Pairwise(m.bin)),
             C_pears = mean(na.omit(spasm::C_pears(m.bin))),
             C_match = mean(spasm::C_match(m.bin)),
             C_w = C_w(m.bin),
             C_ratio = C_ratio(m.bin),
             C_combo = C_combo(m.bin),
             C_checker= C_checker(m.bin),
             C_conn = C_conn(m.bin),
             C_segSc_Z = mean(cleaner(spasm::Z_score(m.bin, "step_C_sim2", "C_seg", N.sim=N))),
             C_togSc_Z = mean(cleaner(spasm::Z_score(m.bin, "step_C_sim2", "C_tog", N.sim=N))),
             C_jacc_Z = mean(cleaner(spasm::Z_score(m.bin, "step_C_sim2", "C_jacc", N.sim=N))),
             C_sor_Z = mean(cleaner(spasm::Z_score(m.bin, "step_C_sim2", "C_sor", N.sim=N))),
             C_match_Z = mean(cleaner(spasm::Z_score(m.bin, "step_C_sim2", "C_match", N.sim=N)))
             )

  res[[names(dat[i])]] <- res.i
}

res <- do.call("rbind", res)
res <- data.frame(res)

# save the results
write.csv(res, file = "empirical_results_binary.csv", row.names=FALSE)

# ------------------------------------------------------------------------------
# read the results
res <- read.csv("empirical_results_binary.csv")

# check distributions of the measures
par(mfrow=c(5,5))
for(i in 1:ncol(res))
{
  hist(res[,i], xlab = names(res)[i], breaks=20, col="grey", main = NULL)
}
apply(X=res, MARGIN=2, FUN=range)

# transformation of some of the measures to make the relationships more symmetrical
res2 <- mutate(res,
               C_forbes = log(C_forbes),
               C_segSc = log(C_segSc + 1),
               C_togSc = log(C_togSc + 1),
               C_w = log(C_w),
               C_checker = log(C_checker + 1),
               C_combo = log(C_combo),
               C_ratio = log(C_ratio),
               N = log(N),
               Gamma = log(Gamma),
               Tot.incid = log(Tot.incid))

res <- rename(res, n = N, gamma = Gamma)

# plot histogram of each variable
par(mfrow=c(5,5))
for(i in 1:ncol(res2))
{
  hist(res2[,i], xlab = names(res2)[i], breaks=20, col="grey", main = NULL)
}

# remove -Inf and NA values
res2[res2 == -Inf] <- NA
res2[res2 == Inf] <- NA
res2 <- na.omit(res2)

# remove Z-scores
res2 <- res2[,(1:ncol(res2) %in% grep(names(res2), pattern="_Z")) == FALSE]

################################################################################
# PAIRPLOT OF THE INCIDENCE-BASED MEASURES
################################################################################

png(file = "figures/binary_pairwise_pairplot2.png", width = 1400, height = 1400, res=150)

    par(xaxt = "n", yaxt = "n")

    panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
    {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      txt <- round(cor(x, y),2)
      plotrix::draw.circle(x = 0.5, y = 0.5, radius = abs(txt)/2,
                           col = grey(1 - abs(txt)),
                           border = NA)
      text(0.5, 0.5, txt, cex = abs(txt) + 0.2, col = "white")
    }

    panel.xy <- function(x, y){points(x, y, cex = 0.1)}

    pairs(res2, upper.panel = panel.cor, lower.panel = panel.xy, gap = 0)

dev.off()

################################################################################
# PCA ORDINATION PLOT OF THE INCIDENCE-BASED MEASURES
################################################################################

res.pca <- prcomp(res2, scale = TRUE, center = TRUE)
cols <- c("S or N", "S or N", "S or N",  rep("", times=ncol(res2)-3))

cols2 <- c("#F8766D", "#F8766D", "#F8766D",
           rep("#FFFFFF", times=ncol(res2)-3))

inc.all <-  factoextra::fviz_pca_var(res.pca, repel=TRUE, label="var",
                                     col.ind = "grey",
                                     col.circle= "darkgrey", fill.var = "white",
                                     col.var = cols2) +
  scale_colour_manual(values=c("#F8766D", "black")) +
  theme_minimal() +
  theme(legend.position='none', plot.title=element_blank()) #+


png(file = "figures/binary_pairwise_PCA.png", width = 1000, height = 1000, res=200)
  inc.all
dev.off()

################################################################################
# NETWORK GRAPH OF THE INCIDENCE-BASED MEASURES
################################################################################

png(file = "figures/binary_pairwise_graph.png", width = 1000, height = 1000, res=300)
    qgraph(cor(res2), layout = "spring",
           labels = colnames(cor(res2)),
           edge.color = c("black"),
           color = cols2,
           label.cex = 1.4,
           shape = "heart",
           node.label.offset = c(0.5, 0.1))
dev.off()


################################################################################
# ABUNDANCE-BASED MEASURES
################################################################################

# load the data from the spasm package
dat <- c(data.Ulrich)

# removing 8 studies with extremely low values (bad for rounding)
good.studies <- lapply(X = data.Ulrich, FUN = max) > 5
dat <- dat[good.studies]

# round the abundnaces to integers (in order for the null models to work)
dat <- lapply(X = dat, FUN = round)

res <- list() # empty container for the results

pos.fun <- function(D) mean(D[D>0]) # mean of the positive values in D
neg.fun <- function(D) mean(D[D<0]) # mean of the negative values in D

N = 200 # how many null model simulations?

for(i in 1:length(dat))
{
  message(i)
  m <- dat[[i]]

  # remove zero columns and rows
  m <- m[rowSums(m) > 0, ]
  m <- m[, colSums(m) > 0]

  # the indices
  res.i <- c(N = sum(colSums(m) >= 1),
             Gamma = sum(rowSums(m) >= 1),
             Tot.abu = sum(m),
             CA_tau = mean(spasm::CA_cov_cor(m, correlation=TRUE, method="kendall")),
             CA_cov = mean(spasm::CA_cov_cor(m, correlation=FALSE)),
             CA_cov_hell = mean(spasm::CA_cov_cor(m, correlation=FALSE, transf="hellinger")),
             CA_cor = mean(spasm::CA_cov_cor(m, correlation=TRUE, method="pearson")),
             CA_cor_hell = mean(spasm::CA_cov_cor(m, correlation=TRUE, transf = "hellinger", method="pearson")),
             CA_bray = mean(CA_bray(m)),
             CA_hell = mean(CA_hell(m)),
             CA_ruz = mean(CA_ruz(m)),
             CA_chi = mean(CA_chi(m)),
             CA_ratio = C_ratio(m),
             CA_cor_hell_Z = mean(cleaner(spasm::Z_score(m, "step_CA_IT", "CA_cov_cor",
                                                         N.sim=N, transf="hellinger", method="pearson"))),
             CA_hell_Z = mean(cleaner(spasm::Z_score(m, "step_CA_IT", "CA_hell", N.sim=N))),
             CA_tau_Z = mean(cleaner(spasm::Z_score(m, "step_CA_IT", "CA_cov_cor", N.sim=N))),
             CA_bray_Z = mean(cleaner(spasm::Z_score(m, "step_CA_IT", "CA_bray", N.sim=N))),
             CA_ruz_Z = mean(cleaner(spasm::Z_score(m, "step_CA_IT", "CA_ruz", N.sim=N))),
             CA_chi_Z = mean(cleaner(spasm::Z_score(m, "step_CA_IT", "CA_chi", N.sim=N)))
             )

  print(res.i)
  res[[names(dat[i])]] <- res.i
}

res <- do.call("rbind", res)
res <- data.frame(res)

# save the results
write.csv(res, row.names=FALSE, file = "empirical_results_abundance.csv")

################################################################################

# load the results
res <- read.csv(file = "empirical_results_abundance.csv")

# transform some of the measures to make the realtionhips more symmetrical
res <- mutate(res,
              N = log(N),
              Gamma = log(Gamma),
              Tot.abu = log(Tot.abu),
              CA_ratio = log(CA_ratio),
              CA_cov = sqrt(CA_cov),
              CA_cov_hell = sqrt(CA_cov_hell))

res <- rename(res, n = N, gamma = Gamma)

# remove Z-scores
res <- res[,(1:ncol(res) %in% grep(names(res), pattern="_Z")) == FALSE]

# remove -Inf and NA values
res[res == -Inf] <- NA
res <- na.omit(res)

cols <- c("S or N", "S or N", "S or N", rep("", times=ncol(res)-3))

# removing extreme values
res <- res[res$CA_cov < 200,]

# define color scheme
cols2 <- c("#F8766D", "#F8766D",  "#F8766D",
           rep("#FFFFFF", times=ncol(res)-3))

# check distributions of the measures
par(mfrow=c(5,4))
for(i in 1:ncol(res))
{
  hist(res[,i], xlab = names(res)[i], breaks=20, col="grey", main = NULL)
}

################################################################################
# PAIRPLOT OF ABUNDANCE-BASED MEASURES
################################################################################


png(file = "figures/abundance_pairwise_pairplot.png", width = 1400, height = 1400, res=150)
    psych::pairs.panels(res, hist.col = "grey",
                 smooth=FALSE, ellipses = FALSE, density = FALSE, scale=TRUE)
dev.off()


png(file = "figures/abundance_pairwise_pairplot2.png", width = 1000, height = 1000, res=150)

    par(xaxt = "n", yaxt = "n")

    panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
    {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      txt <- round(cor(x, y),2)
      plotrix::draw.circle(x = 0.5, y = 0.5, radius = abs(txt)/2,
                           col = grey(1 - abs(txt)),
                           border = NA)
      text(0.5, 0.5, txt, cex = abs(txt) + 0.2, col = "white")
    }

    panel.xy <- function(x, y){points(x, y, cex = 0.1)}
    pairs(res, upper.panel = panel.cor, lower.panel = panel.xy, gap = 0)

dev.off()

################################################################################
# NETWORK GRAPH OF ABUNDANCE-BASED MEASURES
################################################################################


png(file = "figures/abundance_pairwise_graph.png", width = 1200, height = 1200, res=100)

    qgraph(cor(res), layout = "spring",
           labels = colnames(cor(res)),
           edge.color = c("black"),
           color = cols2,
           label.cex = 1.4,
           shape = "heart",
           node.label.offset = c(0.5, 0.1))

dev.off()


################################################################################
# PCA ORDINATION PLOT OF ABUNDANCE-BASED MEASURES
################################################################################

res.pca <- prcomp(res, scale = TRUE, center = TRUE)

png(file = "figures/abundance_pairwise_PCA.png", width = 1000, height = 1000, res=200)

    factoextra::fviz_pca_var(res.pca, repel=TRUE, label="var",
                             col.ind = "grey", #col.var = "black",
                             col.circle= "darkgrey", fill.var = "white",
                             col.var = cols2) +
      scale_colour_manual(values=c("#F8766D", "black")) +
      theme_minimal() + ggtitle("") +
      theme(legend.position="none", plot.title=element_blank())

dev.off()
