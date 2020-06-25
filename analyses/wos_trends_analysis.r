################################################################################
#
# Code for Appendix S2
#
# Petr Keil
#
################################################################################

library(ggplot2)
library(ggrepel) # for non-overlapping plot labels
library(plyr)
library(gridExtra)
library(spasm)
library(plyr)
library(dplyr)
library(tidyverse)
library(reshape2)



############################
# 1. Web of Science search #
############################

# read the raw data
data.WOS <- read.csv("../data/wos_trends.csv")
# export the data as an R object for packaging purpose
save(data.WOS, file = "../data/Web_of_Science_search.RData")


# some minor tweaks
dat <- data.WOS
dat$Type <- as.character(dat$Type)


# create rows for year with 0 papers, since these are not in the original data
years <- 1996:2019
terms <- unique(dat$Term)
x <- expand.grid(years, terms)
names(x) <- c("Year", "Term")
x <- left_join(unique(data.frame(Term = dat$Term, Type = dat$Type)), x, by=  "Term")
dat <- select(dat, - Type)
dat <- join(x, dat, by = c("Term", "Year"))
dat$Count[is.na(dat$Count)] <- 0


# create labels for search terms to be plotted next to the time series
labs <- plyr::ddply(.data = dat,
           .fun=summarize,
           .variables="Term",
           Count = max(Count),
           Year = max(Year),
           Type = Type[1])


# plot temporal trends using ggplot2
P.temp <- ggplot(data = dat, aes(x = Year, y = Count)) +
  geom_point(aes(colour = Type)) +
  geom_line(aes(colour = Type, linetype = Term), size = 1, alpha = 0.5) +
  scale_linetype_manual(values=rep(1, times = length(unique(dat$Term)))) +
  scale_y_log10(breaks=c(1,10,100, 1000)) +
  geom_label_repel(data = labs,
                   aes(x = Year, y = Count, label = Term, colour = Type),
                   nudge_x = 0.1)+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks=c(1995, 2000, 2005, 2010, 2015)) +
  facet_grid(.~Type) +
  ylab("Number of papers")

P.temp


# export the plot as a .pdf
pdf("../analyses/figures/WOS_trends.pdf", width=14, height = 7)
 P.temp
dev.off()





########################################################################
# 2. Manual search through Ecology, Ecography, and American Naturalist #
########################################################################

# convert journal_search.csv file into an .RData object for packaging purposes:
data.journals <- read.csv("../data/journal_search.csv")
save(data.journals, file = "../data/Journal_Search.RData")


# all three journals pooled together
journals <- select(.data = data.journals, -Journal, -First_author, -Title, - Volume, -Pages, -DOI,
                   -EF, -BEF, -ITA, -CTA, -S_comp_change,  -Dist_decay) %>%
  group_by(Year) %>%
  summarise_all(funs(sum)) %>%
  melt(id.vars = c("Year")) %>% data.frame()

journals <- data.frame(journals, Type = "Other")
journals$Type <- as.character(journals$Type)
journals$Type[journals$variable == "ISA"] <- "ISA"
journals$Type[journals$variable == "ISA"] <- "ISA"
reordered.variable <- reorder(journals$variable, journals$value, FUN = median)
journals$variable <- reordered.variable

ISA.only <- journals[journals$variable == "ISA",] %>% select(-variable, -Type)

J.temp.all <- ggplot(data = journals[journals$variable != "ISA",],
                     aes(x = Year, y = value)) +
  geom_point(colour = "black") +
  geom_line(colour = "black") +
  geom_point(data = ISA.only, aes(x = Year, y = value), colour = "grey") +
  geom_line(data = ISA.only, aes(x = Year, y = value), colour = "grey") +
  scale_y_log10(breaks=c(1,10,100, 1000)) +
  scale_x_continuous(breaks=c(1995, 2003, 2011, 2019)) +
  facet_grid(.~variable) + theme_bw() +
  labs(y ="# of papers", title = "(a) All 3 journals pooled") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
J.temp.all



journals <- select(.data = data.journals, -First_author, -Title, - Volume, -Pages, -DOI,
                           -EF, -BEF, -ITA, -CTA, -S_comp_change, -Dist_decay) %>%
   group_by( Journal, Year) %>%
   summarise_all(funs(sum)) %>%
   melt(id.vars = c("Journal", "Year")) %>% data.frame()

journals <- data.frame(journals, Type = "Other")
journals$Type <- as.character(journals$Type)
journals$Type[journals$variable == "ISA"] <- "ISA"
journals$Type[journals$variable == "ISA"] <- "ISA"
journals$variable <- factor(journals$variable, levels = levels(reordered.variable))

ISA.only <- journals[journals$variable == "ISA",] %>% select(-variable, -Type)

J.temp <- ggplot(data = journals[journals$variable != "ISA",],
                 aes(x = Year, y = value)) +
  geom_point(colour = "black") +
  geom_line(colour = "black") +
  geom_point(data = ISA.only, aes(x = Year, y = value), colour = "grey") +
  geom_line(data = ISA.only, aes(x = Year, y = value), colour = "grey") +
  scale_y_log10(breaks=c(1,10,100, 1000)) +
  scale_x_continuous(breaks=c(1995, 2003, 2011, 2019)) +
  facet_grid(Journal~variable) + theme_bw() +
   labs(y ="# of papers", title = "(b) Each journal separaterly") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
J.temp

# export the figures
pdf("../analyses/figures/journal_trends.pdf", width=15, height = 10)
grid.arrange(J.temp.all, J.temp, nrow = 2, ncol=1, heights = c(1, 2.2))
dev.off()



