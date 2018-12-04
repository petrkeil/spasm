hf <- read.csv("hf253-03-trees-2014.csv")

head(hf)
hf <- hf[hf$df.status == "alive", ]
hf <- hf[hf$dbh > 5, ]

plot(hf$gx, hf$gy, cex = 0.2)