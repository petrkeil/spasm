library(psych)
library(tidyverse)
library(gridExtra)
library(devtools)
#install_github(repo="petrkeil/spasm")
library(spasm)
library(viridis)
library(RColorBrewer)
library(qgraph)

# ANALYZING INCIDENCE-BASED MEASURES
dat <- c(data.Ulrich, data.Atmar)



res <- list()
for(i in 1:length(dat))
{
  message(i)
  m <- dat[[i]]
  m.bin <- m
  m.bin[m.bin > 0] <- 1 # convert to binary matrix

  res.i <- c(N = sum(colSums(m.bin) > 0),
             S = sum(rowSums(m.bin) > 0),
             C_segSc = mean( spasm::C_seg(m.bin)),
             C_togSc = mean(spasm::C_tog(m.bin)),
             C_jacc = mean( spasm::C_jacc(m.bin)),
             C_sor = mean( spasm::C_sor(m.bin)),
             C_forbes = mean( spasm::C_forbes(m.bin)),
             C_FETmP = mean(FETmP_Pairwise(m.bin)),
             C_pears = mean(na.omit(spasm::C_pears(m.bin))),
             C_match = mean(spasm::C_match(m.bin)),
             C_w = C_w(m.bin),
             C_ratio = C_ratio(m.bin),
             C_combo = C_combo(m.bin),
             C_checker= C_checker(m.bin),
             C_conn = C_conn(m.bin),
             C_segSc_Z = mean(na.omit(spasm::Z_score(m.bin, "step_C_sim2", "C_seg", N.sim=100))),
             C_togSc_Z = mean(na.omit(spasm::Z_score(m.bin, "step_C_sim2", "C_tog", N.sim=100))),
             C_jacc_Z = mean(na.omit(spasm::Z_score(m.bin, "step_C_sim2", "C_jacc", N.sim=100))),
             C_sor_Z = mean(na.omit(spasm::Z_score(m.bin, "step_C_sim2", "C_sor", N.sim=100))),
             C_match_Z = mean(na.omit(spasm::Z_score(m.bin, "step_C_sim2", "C_match", N.sim=100)))
             )

 # print(res.i)
  res[[names(dat[i])]] <- res.i
}

res <- do.call("rbind", res)
res <- data.frame(res)
write.csv(res, file = "empirical_results_binary.csv", row.names=FALSE)

# ------------------------------------------------------------------------------
res <- read.csv("empirical_results_binary.csv")

# check distributions of the measures
par(mfrow=c(5,4))
for(i in 1:ncol(res))
{
  hist(res[,i], xlab = names(res)[i], breaks=20, col="grey", main = NULL)
}
apply(X=res, MARGIN=2, FUN=range)

res2 <- mutate(res,
               C_forbes = log(C_forbes),
               C_segSc = log(C_segSc + 1),
               C_togSc = log(C_togSc + 1),
               C_w = log(C_w),
               C_checker = log(C_checker + 1),
               C_combo = log(C_combo),
               C_ratio = log(C_ratio),
               N = log(N),
               S = log(S))

# plot histogram of each variable
par(mfrow=c(5,4))
for(i in 1:ncol(res2))
{
  hist(res2[,i], xlab = names(res2)[i], breaks=20, col="grey", main = NULL)
}

# remove -Inf and NA values
res2[res2 == -Inf] <- NA
res2[res2 == Inf] <- NA
res2 <- na.omit(res2)


# ------------------------------------------------------------------------------


png(file = "figures/binary_pairwise_pairplot2.png",
    width = 1400, height = 1400, res=150)
par(xaxt = "n", yaxt = "n")
#pairs.panels(res2, hist.col = "grey",
#             smooth=FALSE, ellipses = FALSE,
##             density = FALSE, scale=TRUE,
#            pch = 1, gap=0.1, cex.cor = 1.4, cex = 0.01)

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

#title("a", adj = 0.01,  cex.main = 2, line = 3)
dev.off()

# ------------------------------------------------------------------------------


res.pca <- prcomp(res2, scale = TRUE, center = TRUE)
cols <- c("S or N", "S or N", rep("", times=ncol(res2)-2))



inc.all <-  factoextra::fviz_pca_var(res.pca, repel=TRUE, label="var",
                                     col.ind = "grey", #col.var = "black",
                                     col.circle= "darkgrey", fill.var = "white",
                                     col.var = cols) +
  scale_colour_manual(values=c("black", "#F8766D")) +
  theme_minimal() +
  theme(legend.position='none', plot.title=element_blank()) #+
 # ggtitle("b")


png(file = "figures/binary_pairwise_PCA.png",
    width = 1000, height = 1000, res=200)
  inc.all
dev.off()

png(file = "figures/binary_pairwise_graph.png",
    width = 1000, height = 1000, res=300)
qgraph(cor(res2), layout = "spring",
       labels = colnames(cor(res2)),
       edge.color = c("black"))#,
       #title="c", title.cex = 2)
dev.off()


################################################################################
# Abundance-based measures

# function removing -Inf, Inf, NaN and NA from a distance matrix
cleaner <- function(D)
{
  D <- D[]
  D <- na.omit(D)
  D <- D[D != Inf]
  D <- D[D != -Inf]
  return(D)
}


# ANALYZING MEAN VALUE
dat <- c(data.Ulrich)

# removing 8 studies with extremely low values (bad for rounding)
good.studies <- lapply(X = data.Ulrich, FUN = max) > 5
dat <- dat[good.studies]
# round the abundnaces to integers (in order for the null models to work)
dat <- lapply(X = dat, FUN = round)


res <- list()

pos.fun <- function(D) mean(D[D>0])
neg.fun <- function(D) mean(D[D<0])

for(i in 1:length(dat))
{
  message(i)
  m <- dat[[i]]
  # remove zero columns and rows
  m <- m[rowSums(m) > 0, ]
  m <- m[, colSums(m) > 0]

  # other statistics
  res.i <- c(N = sum(colSums(m) >= 1),
             S = sum(rowSums(m) >= 1),
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
             CA_cor_hell_Z = mean(cleaner(spasm::Z_score(m, "step_CA_rowrandom", "CA_cov_cor", N.sim=100, transf="hellinger", method="pearson"))),
             CA_hell_Z = mean(cleaner(spasm::Z_score(m, "step_CA_rowrandom", "CA_hell", N.sim=100))),
             CA_tau_Z = mean(cleaner(spasm::Z_score(m, "step_CA_rowrandom", "CA_cov_cor", N.sim=100))),
             CA_bray_Z = mean(cleaner(spasm::Z_score(m, "step_CA_rowrandom", "CA_bray", N.sim=100))),
             CA_ruz_Z = mean(cleaner(spasm::Z_score(m, "step_CA_rowrandom", "CA_ruz", N.sim=100))),
             CA_chi_Z = mean(cleaner(spasm::Z_score(m, "step_CA_rowrandom", "CA_chi", N.sim=100)))
             )

  print(res.i)
  res[[names(dat[i])]] <- res.i
}

res <- do.call("rbind", res)
res <- data.frame(res)


res <- mutate(res,
              N = log(N),
              S = log(S),
              Tot.abu = log(Tot.abu),
              CA_ratio = log(CA_ratio))





# remove -Inf and NA values
res[res == -Inf] <- NA
res <- na.omit(res)

cols <- c("S or N", "S or N", "S or N", rep("", times=ncol(res)-3))



# check distributions of the measures
par(mfrow=c(5,4))
for(i in 1:ncol(res))
{
  hist(res[,i], xlab = names(res)[i], breaks=20, col="grey", main = NULL)
}

require(psych)
png(file = "figures/abundance_pairwise_pairplot.png",
    width = 1400, height = 1400, res=150)
pairs.panels(res, hist.col = "grey",
             smooth=FALSE, ellipses = FALSE, density = FALSE, scale=TRUE)
dev.off()


png(file = "figures/abundance_pairwise_pairplot2.png",
    width = 1000, height = 1000, res=150)
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



png(file = "figures/abundance_pairwise_graph.png",
    width = 1000, height = 1000, res=300)
qgraph(cor(res), layout = "spring",
       labels = colnames(cor(res)),
       edge.color = c("black"))#,
#title="c", title.cex = 2)
dev.off()



res.pca <- prcomp(res, scale = TRUE, center = TRUE)

png(file = "figures/abundance_pairwise_PCA.png",
    width = 1000, height = 1000, res=200)
factoextra::fviz_pca_var(res.pca, repel=TRUE, label="var",
                         col.ind = "grey", #col.var = "black",
                         col.circle= "darkgrey", fill.var = "white",
                         col.var = cols) +
  scale_colour_manual(values=c("black", "#F8766D")) +
  theme_minimal() + ggtitle("") +
  theme(legend.position="none", plot.title=element_blank())

dev.off()
