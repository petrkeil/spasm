################################################################################
#
# Code for Figure S1, the analysis of Web of Science trends
#
# Petr Keil
#
################################################################################

library(ggplot2)
library(ggrepel) # for non-overlapping plot labels
library(plyr)
library(gridExtra)

# ------------------------------------------------------------------------------
# read the data
data.WOS <- read.csv("../data/wos_trends.csv")

# export the data for packaging purpose
save(data.WOS, file = "../data/Web_of_Science_search.RData")

dat <- data.WOS # this is what I am going to work with

dat$Type <- as.character(dat$Type)

# ------------------------------------------------------------------------------

# create labels for search terms
labs <- plyr::ddply(.data = dat,
           .fun=summarize,
           .variables="Term",
           Count = max(Count),
           Year = max(Year),
           Type = Type[1])

ISA <- dat[dat$Type == "ISA",]
other <- dat[dat$Type == "other",]
interaction <- dat[dat$Type == "interaction",]

# ------------------------------------------------------------------------------

# plot temporal trends
P.temp <- ggplot(data = dat, aes(x = Year, y = Count)) +
  geom_point(aes(colour = Type)) +
  geom_line(aes(colour = Type, linetype = Term), size = 1, alpha = 0.5) +
  scale_linetype_manual(values=rep(1, times = length(unique(dat$Term)))) +
  scale_y_log10(breaks=c(1,10,100, 1000)) +
  geom_label_repel(data = labs, aes(x = Year, y = Count, label = Term, colour = Type),
                   nudge_x = 0.1)+
  theme_bw()+
  theme(legend.position = "none") +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks=c(1995, 2000, 2005, 2010, 2015)) +
  facet_grid(.~Type) +
  ylab("Number of papers")

# export the plot
pdf("../analyses/figures/WOS_trends.pdf", width=12, height = 7)
 P.temp
dev.off()


