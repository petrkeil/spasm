library(ggplot2)

x = 3
gamma <- c(1,10,100,1000)
tot <- (gamma*(gamma - 1))/2
sig <- gamma*x

dat <- dara.frame(gamma, tot, sig)


