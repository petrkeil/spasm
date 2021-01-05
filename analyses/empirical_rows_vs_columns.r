################################################################################
#
# Analysis comparing the relationship between mean Jaccard index calculated over
# rows of the community matrix, with mean pairwise Jaccard calculated over columns.
# Compared is (a) the raw index, (b) the Z-score calculated using the sim2 algorithm from
# EcoSimR, and (c) the IT algorithm from Ulrich & Gotelli ()
#
# Author: Petr Keil
#
################################################################################

library(tidyverse)
library(spasm)
library(gridExtra)

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
#dat <- c(data.Atmar)


N = 400 # how many null model simulations?
set.seed(12345) # set seed for reproducibility

res <- list()
for(i in 1:length(dat))
{
  message(i)
  m <- dat[[i]]
  m.bin <- m
  m.bin[m.bin > 0] <- 1 # convert to binary matrix

  # calculate the metrics
  res.i <- c(n = sum(colSums(m.bin) > 0),
             S = sum(rowSums(m.bin) > 0),
             Tot.incid = sum(m.bin),
             C_jacc = mean( spasm::C_jacc(m.bin)),
             C_jacc_Z_sim2 = mean(cleaner(spasm::Z_score(m.bin, "step_C_sim2", "C_jacc", N.sim=N))),
             C_jacc_Z_CA = mean(cleaner(spasm::Z_score(m.bin, "step_CA_IT", "C_jacc", N.sim=N))),
             B_jacc =  mean( spasm::C_jacc(t(m.bin))),
             B_jacc_Z_sim2 = mean(cleaner(spasm::Z_score(t(m.bin), "step_C_sim2", "C_jacc", N.sim=N))),
             B_jacc_Z_CA = mean(cleaner(spasm::Z_score(t(m.bin), "step_CA_IT", "C_jacc", N.sim=N))))

  res[[names(dat[i])]] <- res.i
}

res <- do.call("rbind", res)
res <- data.frame(res)

#write.csv(res, file = "empirical_results_row_vs_col_jaccard.csv", row.names=FALSE)

res <- read.csv(file = "empirical_results_row_vs_col_jaccard.csv")


# --------------------------------------------------------------
# PLOT THE RESULTS

# raw Jaccard indices
raw.plot <- ggplot(data = res, aes(x = C_jacc, y = B_jacc)) +
  geom_point(shape = 1) +
  geom_abline(intercept=0, slope=1) +
  theme_bw() +#gtitle(label=paste ("(a)")) +
  labs(title = "(a)",
       subtitle = paste("r = ", round(cor(res$C_jacc, res$B_jacc), 2))) +
  xlab("C_Jaccard") + ylab("Beta_Jaccard")

# Z-scores using the sim2 algorithm
Z.plot.sim2 <- ggplot(data = res, aes(x = C_jacc_Z_sim2, y = B_jacc_Z_sim2)) +
  geom_vline(xintercept=0, colour = "grey") + geom_hline(yintercept=0, colour = "grey") +
  geom_point(shape = 1) +
  geom_abline(intercept=0, slope=1) +
  theme_bw()+
  labs(title = "(b)",
       subtitle = paste("r = ", round(cor(res$C_jacc_Z_sim2, res$B_jacc_Z_sim2, use = "pairwise.complete.obs"), 2))) +
  xlab("Z-score of C_Jaccard (sim2)") + ylab("Z-score of Beta_Jaccard (sim2)")

# Z-scores using the IT algorithm
Z.plot.CA <- ggplot(data = res, aes(x = C_jacc_Z_CA, y = B_jacc_Z_CA)) +
  geom_vline(xintercept=0, colour = "grey") + geom_hline(yintercept=0, colour = "grey") +
  geom_point(shape = 1) +
  geom_abline(intercept=0, slope=1) +
  theme_bw()+
  labs(title = "(c)",
       subtitle = paste("r = ", round(cor(res$C_jacc_Z_CA, res$B_jacc_Z_CA), 2))) +
  xlab("Z-score of C_Jaccard (IT)") + ylab("Z-score of Beta_Jaccard (IT)")

# export the figure
pdf("figures/rows_vs_columns.pdf", width=11, height= 4)
grid.arrange(raw.plot, Z.plot.sim2, Z.plot.CA, ncol = 3)
dev.off()






