library(psych)
library(tidyverse)
library(gridExtra)
library(devtools)
#install_github(repo="petrkeil/spasm")
library(spasm)
library(viridis)
library(RColorBrewer)
library(qgraph)

my.functions <- c("C_forbes", "C_pears", "C_jacc", "C_sor", "C_sim",
                  "C_w", "C_checker", "C_combo", "C_conn", "C_ratio")

# ANALYZING INCIDENCE-BASED MEASURES
dat <- c(data.Ulrich, data.Atmar)



res <- list()
for(i in 1:length(dat))
{
  message(i)
  m <- dat[[i]]

  # the distance matrix-based measures
  res.i <- list()
  for(j in my.functions)
  {
    res.i[[j]] <- mean(do.call(j, list(m = m)))
  }
  res.i <- unlist(res.i)

  # other statistics
  res.i <- c(N = sum(colSums(m) >= 1),
             S = sum(rowSums(m) >= 1),
             C_togSc = mean(C_tog(m, scale= TRUE)),
             C_segSc = mean(C_seg(m, scale = TRUE)),
             res.i)

 # print(res.i)
  res[[names(dat[i])]] <- res.i
}

res <- do.call("rbind", res)
res <- data.frame(res)

# check distributions of the measures
par(mfrow=c(4,4))
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
par(mfrow=c(4,4))
for(i in 1:ncol(res2))
{
  hist(res2[,i], xlab = names(res2)[i], breaks=20, col="grey", main = NULL)
}

# remove -Inf and NA values
res2[res2 == -Inf] <- NA
res2 <- na.omit(res2)


# ------------------------------------------------------------------------------


png(file = "figures/binary_pairwise_pairplot2.png",
    width = 1000, height = 1000, res=150)
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


png(file = "figures/binary_corrplot.png", width = 1000, height = 1000, res=200)
corrplot::corrplot.mixed(corr = cor(res2), lower="ellipse",
                         upper = "number", order="FPC",
                         tl.col = "black",
                         tl.cex = 0.7,
                         tl.pos = "lt",
                         tl.srt = 40,
                         number.cex = 0.8,
                         lower.col = viridis(200),
                         upper.col = viridis(200))
title("a", adj = 0,  cex.main = 2, line = 2.5)
dev.off()


res.pca <- prcomp(res2, scale = TRUE, center = TRUE)
cols <- c("S or N", "S or N", rep("", times=ncol(res2)-2))



inc.all <-  factoextra::fviz_pca_var(res.pca, repel=TRUE, label="var",
                                     col.ind = "grey", #col.var = "black",
                                     col.circle= "darkgrey", fill.var = "white",
                                     col.var = cols) +
  scale_colour_manual(values=c("black", "red")) +
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
       theme = "colorblind")#,
       #title="c", title.cex = 2)
dev.off()


################################################################################
# Abundance-based measures


# ANALYZING MEAN VALUE
dat <- c(data.Ulrich)
res <- list()

pos.fun <- function(D) mean(D[D>0])
neg.fun <- function(D) mean(D[D<0])


for(i in 1:length(dat))
{
  message(i)
  m <- dat[[i]]

  # other statistics
  res.i <- c(N = sum(colSums(m) >= 1),
             S = sum(rowSums(m) >= 1),
             Tot.abu = sum(m),
            #CA_gower = CA_gower(m, fun = mean),
             CA_bray = CA_bray(m, fun = mean),
             CA_hell = CA_hell(m, fun = mean),
             CA_ruz = CA_ruz(m, fun=mean),
             CA_chi = CA_chi(m, fun=mean),
             CA_cor = CA_cov_cor(m, fun = mean,
                                    transf = "hellinger",
                                    correlation = TRUE,
                                    method = "pearson"),
             CA_tau = CA_cov_cor(m, fun=mean, method = "kendall", correlation=TRUE),
             #CA_pos = CA_cov_cor(m, fun = pos.fun,
            #                              transf = "hellinger",
             #                             correlation = TRUE,
            #                              method = "pearson"),
            # CA_neg = CA_cov_cor(m, fun = neg.fun,
            #                              transf = "hellinger",
            #                              correlation = TRUE,
            #                              method = "pearson",
            CA_ratio = C_ratio(m))

  # print(res.i)
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
par(mfrow=c(3,4))
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
       labels = colnames(cor(res)), theme = "colorblind")#,
#title="c", title.cex = 2)
dev.off()



res.pca <- prcomp(res, scale = TRUE, center = TRUE)

png(file = "figures/abundance_pairwise_PCA.png",
    width = 1000, height = 1000, res=200)
factoextra::fviz_pca_var(res.pca, repel=TRUE, label="var",
                         col.ind = "grey", #col.var = "black",
                         col.circle= "darkgrey", fill.var = "white",
                         col.var = cols) +
  scale_colour_manual(values=c("black", "red")) +
  theme_minimal() + ggtitle("") +
  theme(legend.position="none", plot.title=element_blank())

dev.off()
