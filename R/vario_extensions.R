
# ------------------------------------------------------------------------------
vario.summarize <- function(v.obs, v.rand)
{
  S <- v.rand$parms$S

  vr.pos <- v.rand$vario[,2,]
  vr.pos.q <- t(apply(vr.pos, MARGIN = 1, FUN = quantile, c(0.025, 0.975)))
  vr.pos.q <- data.frame(vr.pos.q)
  names(vr.pos.q) <- c("Pos.min", "Pos.max")

  vr.neg <- v.rand$vario[,3,]
  vr.neg.q <- t(apply(vr.neg, MARGIN = 1, FUN = quantile, c(0.025, 0.975)))
  vr.neg.q <- data.frame(vr.neg.q)
  names(vr.neg.q) <- c("Neg.min", "Neg.max")

  vr.con <- v.rand$vario[,1,]
  vr.con.q <- t(apply(vr.con, MARGIN = 1, FUN = quantile, c(0.025, 0.975)))
  vr.con.q <- data.frame(vr.con.q)
  names(vr.con.q) <- c("Con.min", "Con.max")

  vr.all <- data.frame(H = v.obs$vario$H,
                       Dist = v.obs$vario$Dist,
                       n = v.obs$vario$n,
                       S = v.obs$parms$S,
                       N = v.obs$parms$N,
                       Pos.covar = v.obs$vario$pos,
                       vr.pos.q,
                       Neg.covar = v.obs$vario$neg,
                       vr.neg.q,
                       Con.var = v.obs$vario$exp.var,
                       vr.con.q)
  return(vr.all)
}
# vr.all <- summarize.vario(v.obs, v.rand)

# ------------------------------------------------------------------------------
vario.plot <- function(vr.all)
{

  pos.color <- "#4daf4a"
  neg.color <- "#e41a1c"
  con.color <- "#377eb8"

  inter <- ggplot(data = vr.all, aes(x = Dist, y = Pos.covar)) +
    geom_hline(yintercept = 0, colour = "darkgrey") +
    geom_ribbon(aes(x = Dist, ymin = Pos.min, ymax = Pos.max),
                fill = pos.color, alpha = 0.5) +
    geom_ribbon(aes(x = Dist, ymin = Neg.min, ymax = Neg.max),
                fill = neg.color, alpha = 0.5) +
    geom_line(colour = pos.color, size= 1) +
    geom_line(aes(x = Dist, y = Neg.covar), colour = neg.color, size = 1) +
    ylab("Covariance") + xlab("Distance") +
    #ylim(c(-50, 50)) +
    ggtitle("Interspecific aggregation and seggregation") +
    theme_bw()

  con <- ggplot(data = vr.all, aes(x = Dist, y = Con.var)) +
    geom_ribbon(aes(x = Dist, ymin = Con.min, ymax = Con.max),
                fill = con.color, alpha = 0.5) +
    geom_line(colour = con.color, size= 1) +
    ylab("Variance") + xlab("Distance") +
    #ylim(c(0, 120)) +
    ggtitle("Conspecific aggregation") +
    theme_bw()

  gridExtra::grid.arrange(con, inter, ncol=2)

}

#plot.vr(vr.all)


