library(mcmcOutput)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

setwd("C:/Users/Quresh.Latif/files/projects/BLM/ARIM")
load("Data_compiled.RData")

#__________ Script inputs _____________#
mod.nam <- "mod_path"
scripts.loc <- "ARIM-analysis/"
mod <- loadObject(mod.nam)
nsims <- dim(mod$mcmcOutput)[1]

  # Data processing
source(str_c(scripts.loc, "Data_processing.R"))
#______________________________________#


## Stratum covariate values ##
source(str_c(scripts.loc, "Calculate_mech_path_covariate_values.R"))

## Guild membership ##
guilds <- c("All", sort(unique(spp.out$Guild)))
guild.mem <- matrix(NA, nrow = length(spp.list), ncol = length(guilds),
                    dimnames = list(spp.list, guilds))
guild.mem[,"All"] <- TRUE
for(g in 2:length(guilds)) guild.mem[,g] <- spp.out$Guild == guilds[g]
# Classify sagebrush species also as shrubland species.
guild.mem[which(guild.mem[,"Sagebrush"]), "Shrubland"] <- TRUE

#__________ Table compilation functions __________#
plot.table.grid.fn <- function(w, B0, B1, X.B.hi, X.B.lo, X.B.bg, dev.B,
                          D0, D1, X.D.hi, X.D.lo, X.D.bg, X.trend) {
  # Setup #
  out.plot <- expand.grid(Development = c("High", "Low", "Background"),
                          Year = 2010:2019, stringsAsFactors = F) %>%
    mutate(SR.yr = 0,
           SR.yr.lo = 0,
           SR.yr.hi = 0,
           SR.pred = 0,
           SR.pred.lo = 0,
           SR.pred.hi = 0)
  nsp <- dim(B0)[[2]]
  nyr <- dim(dev.B)[[3]]

  # High development #
  BETA <- B0
  DELTA <- D0
  for(sp in 1:nsp) {
    BETA[,sp] <- BETA[,sp] + apply(B1[,sp,] * X.B.hi, 1, sum)
    DELTA[,sp] <- DELTA[,sp] + apply(D1[,sp,] * X.D.hi, 1, sum)
  }
  PSI.pred <- PSI.yr <- array(NA, dim = c(nsims, nsp, nyr))
  for(t in 1:nyr) {
    PSI.pred[,,t] <- QSLpersonal::expit(BETA + DELTA * X.trend[t]) * w
    PSI.yr[,,t] <- QSLpersonal::expit(BETA + dev.B[,,t] +
                                       DELTA * X.trend[t]) * w
  }
  out.plot$SR.yr[which(out.plot$Development == "High")] <-
    PSI.yr %>% apply(c(1, 3), sum) %>% apply(2, median)
  out.plot$SR.yr.lo[which(out.plot$Development == "High")] <-
    PSI.yr %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.1, type = 8)
  out.plot$SR.yr.hi[which(out.plot$Development == "High")] <-
    PSI.yr %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.9, type = 8)
  out.plot$SR.pred[which(out.plot$Development == "High")] <-
    PSI.pred %>% apply(c(1, 3), sum) %>% apply(2, median)
  out.plot$SR.pred.lo[which(out.plot$Development == "High")] <-
    PSI.pred %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.1, type = 8)
  out.plot$SR.pred.hi[which(out.plot$Development == "High")] <-
    PSI.pred %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.9, type = 8)

  # Low development #
  BETA <- B0
  DELTA <- D0
  for(sp in 1:nsp) {
    BETA[,sp] <- BETA[,sp] + apply(B1[,sp,] * X.B.lo, 1, sum)
    DELTA[,sp] <- DELTA[,sp] + apply(D1[,sp,] * X.D.lo, 1, sum)
  }
  PSI.pred <- PSI.yr <- array(NA, dim = c(nsims, nsp, nyr))
  for(t in 1:nyr) {
    PSI.pred[,,t] <- QSLpersonal::expit(BETA + DELTA * X.trend[t]) * w
    PSI.yr[,,t] <- QSLpersonal::expit(BETA + dev.B[,,t] +
                                        DELTA * X.trend[t]) * w
  }
  out.plot$SR.yr[which(out.plot$Development == "Low")] <-
    PSI.yr %>% apply(c(1, 3), sum) %>% apply(2, median)
  out.plot$SR.yr.lo[which(out.plot$Development == "Low")] <-
    PSI.yr %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.1, type = 8)
  out.plot$SR.yr.hi[which(out.plot$Development == "Low")] <-
    PSI.yr %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.9, type = 8)
  out.plot$SR.pred[which(out.plot$Development == "Low")] <-
    PSI.pred %>% apply(c(1, 3), sum) %>% apply(2, median)
  out.plot$SR.pred.lo[which(out.plot$Development == "Low")] <-
    PSI.pred %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.1, type = 8)
  out.plot$SR.pred.hi[which(out.plot$Development == "Low")] <-
    PSI.pred %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.9, type = 8)
  
  # Background #
  BETA <- B0
  DELTA <- D0
  for(sp in 1:nsp) {
    BETA[,sp] <- BETA[,sp] + apply(B1[,sp,] * X.B.bg, 1, sum)
    DELTA[,sp] <- DELTA[,sp] + apply(D1[,sp,] * X.D.bg, 1, sum)
  }
  PSI.pred <- PSI.yr <- array(NA, dim = c(nsims, nsp, nyr))
  for(t in 1:nyr) {
    PSI.pred[,,t] <- QSLpersonal::expit(BETA + DELTA * X.trend[t]) * w
    PSI.yr[,,t] <- QSLpersonal::expit(BETA + dev.B[,,t] +
                                        DELTA * X.trend[t]) * w
  }
  out.plot$SR.yr[which(out.plot$Development == "Background")] <-
    PSI.yr %>% apply(c(1, 3), sum) %>% apply(2, median)
  out.plot$SR.yr.lo[which(out.plot$Development == "Background")] <-
    PSI.yr %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.1, type = 8)
  out.plot$SR.yr.hi[which(out.plot$Development == "Background")] <-
    PSI.yr %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.9, type = 8)
  out.plot$SR.pred[which(out.plot$Development == "Background")] <-
    PSI.pred %>% apply(c(1, 3), sum) %>% apply(2, median)
  out.plot$SR.pred.lo[which(out.plot$Development == "Background")] <-
    PSI.pred %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.1, type = 8)
  out.plot$SR.pred.hi[which(out.plot$Development == "Background")] <-
    PSI.pred %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.9, type = 8)
  
  return(out.plot)
}

plot.table.point.fn <- function(w, B0, B1, X.B.hi, X.B.lo, X.B.bg, dev.B,
                               D0, D1, X.D.hi, X.D.lo, X.D.bg,
                               b0, b1, X.b.hi, X.b.lo, X.b.bg, dev.b,
                               d0, d1, X.d.hi, X.d.lo, X.d.bg, X.trend) {
  # Setup #
  out.plot <- expand.grid(Development = c("High", "Low", "Background"),
                          Year = 2010:2019, stringsAsFactors = F) %>%
    mutate(SR.yr = 0,
           SR.yr.lo = 0,
           SR.yr.hi = 0,
           SR.pred = 0,
           SR.pred.lo = 0,
           SR.pred.hi = 0)
  nsp <- dim(B0)[[2]]
  nyr <- dim(dev.B)[[3]]
  
  # High development #
  BETA <- B0
  DELTA <- D0
  for(sp in 1:nsp) {
    BETA[,sp] <- BETA[,sp] + apply(B1[,sp,] * X.B.hi, 1, sum)
    DELTA[,sp] <- DELTA[,sp] + apply(D1[,sp,] * X.D.hi, 1, sum)
  }
  PSI.pred <- PSI.yr <- array(NA, dim = c(nsims, nsp, nyr))
  for(t in 1:nyr) {
    PSI.pred[,,t] <- QSLpersonal::expit(BETA + DELTA * X.trend[t]) * w
    PSI.yr[,,t] <- QSLpersonal::expit(BETA + dev.B[,,t] +
                                        DELTA * X.trend[t]) * w
  }
  beta <- b0
  delta <- d0
  for(sp in 1:nsp) {
    beta[,sp] <- beta[,sp] + apply(b1[,sp,] * X.b.hi, 1, sum)
    delta[,sp] <- delta[,sp] + apply(d1[,sp,] * X.d.hi, 1, sum)
  }
  psi.pred <- psi.yr <- array(NA, dim = c(nsims, nsp, nyr))
  for(t in 1:nyr) {
    psi.pred[,,t] <- QSLpersonal::expit(beta + delta * X.trend[t])
    psi.yr[,,t] <- QSLpersonal::expit(beta + dev.b[,,t] +
                                        delta * X.trend[t])
  }
  out.plot$SR.yr[which(out.plot$Development == "High")] <-
    (PSI.yr * psi.yr) %>% apply(c(1, 3), sum) %>% apply(2, median)
  out.plot$SR.yr.lo[which(out.plot$Development == "High")] <-
    (PSI.yr * psi.yr) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.1, type = 8)
  out.plot$SR.yr.hi[which(out.plot$Development == "High")] <-
    (PSI.yr * psi.yr) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.9, type = 8)
  out.plot$SR.pred[which(out.plot$Development == "High")] <-
    (PSI.pred * psi.pred) %>% apply(c(1, 3), sum) %>% apply(2, median)
  out.plot$SR.pred.lo[which(out.plot$Development == "High")] <-
    (PSI.pred * psi.pred) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.1, type = 8)
  out.plot$SR.pred.hi[which(out.plot$Development == "High")] <-
    (PSI.pred * psi.pred) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.9, type = 8)
  
  # Low development #
  BETA <- B0
  DELTA <- D0
  for(sp in 1:nsp) {
    BETA[,sp] <- BETA[,sp] + apply(B1[,sp,] * X.B.lo, 1, sum)
    DELTA[,sp] <- DELTA[,sp] + apply(D1[,sp,] * X.D.lo, 1, sum)
  }
  PSI.pred <- PSI.yr <- array(NA, dim = c(nsims, nsp, nyr))
  for(t in 1:nyr) {
    PSI.pred[,,t] <- QSLpersonal::expit(BETA + DELTA * X.trend[t]) * w
    PSI.yr[,,t] <- QSLpersonal::expit(BETA + dev.B[,,t] +
                                        DELTA * X.trend[t]) * w
  }
  beta <- b0
  delta <- d0
  for(sp in 1:nsp) {
    beta[,sp] <- beta[,sp] + apply(b1[,sp,] * X.b.lo, 1, sum)
    delta[,sp] <- delta[,sp] + apply(d1[,sp,] * X.d.lo, 1, sum)
  }
  psi.pred <- psi.yr <- array(NA, dim = c(nsims, nsp, nyr))
  for(t in 1:nyr) {
    psi.pred[,,t] <- QSLpersonal::expit(beta + delta * X.trend[t])
    psi.yr[,,t] <- QSLpersonal::expit(beta + dev.b[,,t] +
                                        delta * X.trend[t])
  }
  out.plot$SR.yr[which(out.plot$Development == "Low")] <-
    (PSI.yr * psi.yr) %>% apply(c(1, 3), sum) %>% apply(2, median)
  out.plot$SR.yr.lo[which(out.plot$Development == "Low")] <-
    (PSI.yr * psi.yr) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.1, type = 8)
  out.plot$SR.yr.hi[which(out.plot$Development == "Low")] <-
    (PSI.yr * psi.yr) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.9, type = 8)
  out.plot$SR.pred[which(out.plot$Development == "Low")] <-
    (PSI.pred * psi.pred) %>% apply(c(1, 3), sum) %>% apply(2, median)
  out.plot$SR.pred.lo[which(out.plot$Development == "Low")] <-
    (PSI.pred * psi.pred) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.1, type = 8)
  out.plot$SR.pred.hi[which(out.plot$Development == "Low")] <-
    (PSI.pred * psi.pred) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.9, type = 8)
  
  # Background #
  BETA <- B0
  DELTA <- D0
  for(sp in 1:nsp) {
    BETA[,sp] <- BETA[,sp] + apply(B1[,sp,] * X.B.bg, 1, sum)
    DELTA[,sp] <- DELTA[,sp] + apply(D1[,sp,] * X.D.bg, 1, sum)
  }
  PSI.pred <- PSI.yr <- array(NA, dim = c(nsims, nsp, nyr))
  for(t in 1:nyr) {
    PSI.pred[,,t] <- QSLpersonal::expit(BETA + DELTA * X.trend[t]) * w
    PSI.yr[,,t] <- QSLpersonal::expit(BETA + dev.B[,,t] +
                                        DELTA * X.trend[t]) * w
  }
  beta <- b0
  delta <- d0
  for(sp in 1:nsp) {
    beta[,sp] <- beta[,sp] + apply(b1[,sp,] * X.b.bg, 1, sum)
    delta[,sp] <- delta[,sp] + apply(d1[,sp,] * X.d.bg, 1, sum)
  }
  psi.pred <- psi.yr <- array(NA, dim = c(nsims, nsp, nyr))
  for(t in 1:nyr) {
    psi.pred[,,t] <- QSLpersonal::expit(beta + delta * X.trend[t])
    psi.yr[,,t] <- QSLpersonal::expit(beta + dev.b[,,t] +
                                        delta * X.trend[t])
  }
  out.plot$SR.yr[which(out.plot$Development == "Background")] <-
    (PSI.yr * psi.yr) %>% apply(c(1, 3), sum) %>% apply(2, median)
  out.plot$SR.yr.lo[which(out.plot$Development == "Background")] <-
    (PSI.yr * psi.yr) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.1, type = 8)
  out.plot$SR.yr.hi[which(out.plot$Development == "Background")] <-
    (PSI.yr * psi.yr) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.9, type = 8)
  out.plot$SR.pred[which(out.plot$Development == "Background")] <-
    (PSI.pred * psi.pred) %>% apply(c(1, 3), sum) %>% apply(2, median)
  out.plot$SR.pred.lo[which(out.plot$Development == "Background")] <-
    (PSI.pred * psi.pred) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.1, type = 8)
  out.plot$SR.pred.hi[which(out.plot$Development == "Background")] <-
    (PSI.pred * psi.pred) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.9, type = 8)
  
  return(out.plot)
}
#_____________________________________________#

# Generate plots for each guild #
for(g in 1:length(guilds)[-which(guilds == "Shrubland")]) {
  gld <- guilds[g]
  gld.nam <- str_c(gld, " species")
  ind.spp <- which(guild.mem[,gld])
  
  # Grid cell richness #
  dat.plot <- plot.table.grid.fn(w = mod$mcmcOutput$w[,ind.spp],
                                 B0 = mod$mcmcOutput$BETA0[,ind.spp],
                                 B1 = mod$mcmcOutput$BETA1[,ind.spp,],
                                 X.B.hi = X.PSI.pred.hi,
                                 X.B.lo = X.PSI.pred.lo,
                                 X.B.bg = X.PSI.pred.bg,
                                 dev.B = mod$mcmcOutput$dev.BETA[,ind.spp,],
                                 D0 = mod$mcmcOutput$DELTA0[,ind.spp],
                                 D1 = mod$mcmcOutput$DELTA1[,ind.spp,],
                                 X.D.hi = X.LAMBDA.pred.hi,
                                 X.D.lo = X.LAMBDA.pred.lo,
                                 X.D.bg = X.LAMBDA.pred.bg,
                                 X.trend = X.trend) %>%
    mutate(Development = factor(Development, levels = c("High", "Low", "Background")),
           Year.jitter = ifelse(Development == "High", Year + 0.1,
                                ifelse(Development == "Background", Year - 0.1,
                                       Year)))
  
  p.grid <- ggplot(dat.plot, aes(x = Year, y = SR.pred)) +
    geom_ribbon(aes(ymin = SR.pred.lo, ymax = SR.pred.hi, fill = Development),
                size = 0, alpha = 0.3) +
    geom_line(aes(color = Development), size = 1) +
    geom_errorbar(aes(x = Year.jitter, ymin = SR.yr.lo, ymax = SR.yr.hi,
                      color = Development), width = 0.1) +
    geom_point(aes(x = Year.jitter, y = SR.yr,
                   color = Development, shape = Development)) +
    scale_color_manual(values = c("#D55E00", "#009E73", "#000000")) +
    scale_fill_manual(values = c("#D55E00", "#009E73", "#000000")) +
    scale_shape_manual(values = c(15, 17, 16)) +
    scale_x_continuous(breaks = seq(2010, 2019, by = 2))
  if(gld == "All") {
    p.grid <- p.grid +
      ylim(5,26) +
      theme(legend.position = c(1,0), legend.justification = c(1,0)) +
      xlab(NULL) + ylab("Grid-cell richness")  
  } else {
    p.grid <- p.grid +
      guides(color = F, fill = F, shape = F) +
      xlab(NULL) + ylab("Grid-cell richness")  
  }

  # Point occupancy #
  dat.plot <- plot.table.point.fn(w = mod$mcmcOutput$w[,ind.spp],
                            B0 = mod$mcmcOutput$BETA0[,ind.spp],
                            B1 = mod$mcmcOutput$BETA1[,ind.spp,],
                            X.B.hi = X.PSI.pred.hi,
                            X.B.lo = X.PSI.pred.lo,
                            X.B.bg = X.PSI.pred.bg,
                            dev.B = mod$mcmcOutput$dev.BETA[,ind.spp,],
                            D0 = mod$mcmcOutput$DELTA0[,ind.spp],
                            D1 = mod$mcmcOutput$DELTA1[,ind.spp,],
                            X.D.hi = X.LAMBDA.pred.hi,
                            X.D.lo = X.LAMBDA.pred.lo,
                            X.D.bg = X.LAMBDA.pred.bg,
                            b0 = mod$mcmcOutput$beta0[,ind.spp],
                            b1 = mod$mcmcOutput$beta1[,ind.spp,],
                            X.b.hi = X.psi.pred.hi,
                            X.b.lo = X.psi.pred.lo,
                            X.b.bg = X.psi.pred.bg,
                            dev.b = mod$mcmcOutput$dev.beta[,ind.spp,],
                            d0 = mod$mcmcOutput$delta0[,ind.spp],
                            d1 = mod$mcmcOutput$delta1[,ind.spp,],
                            X.d.hi = X.lambda.pred.hi,
                            X.d.lo = X.lambda.pred.lo,
                            X.d.bg = X.lambda.pred.bg,
                            X.trend = X.trend) %>%
    mutate(Development = factor(Development, levels = c("High", "Low", "Background")),
           Year.jitter = ifelse(Development == "High", Year + 0.1,
                                ifelse(Development == "Background", Year - 0.1,
                                       Year)))
  
  p.point <- ggplot(dat.plot, aes(x = Year, y = SR.pred)) +
    geom_ribbon(aes(ymin = SR.pred.lo, ymax = SR.pred.hi, fill = Development),
                size = 0, alpha = 0.3) +
    geom_line(aes(color = Development), size = 1) +
    geom_errorbar(aes(x = Year.jitter, ymin = SR.yr.lo, ymax = SR.yr.hi,
                      color = Development), width = 0.1) +
    geom_point(aes(x = Year.jitter, y = SR.yr,
                   color = Development, shape = Development)) +
    scale_color_manual(values = c("#D55E00", "#009E73", "#000000")) +
    scale_fill_manual(values = c("#D55E00", "#009E73", "#000000")) +
    scale_shape_manual(values = c(15, 17, 16)) +
    scale_x_continuous(breaks = seq(2010, 2019, by = 2)) +
    #theme(legend.position = c(1,0), legend.justification = c(1,0)) +
    guides(color = F, fill = F, shape = F) +
    xlab(NULL) + ylab("Point richness")
  
  # Put everything together
  p.gld <- ggdraw() +
    draw_plot(p.grid,  x = 0,   y = 0, width = 0.5, height = 0.95) +
    draw_plot(p.point, x = 0.5, y = 0, width = 0.5, height = 0.95) +
    draw_plot_label(gld.nam, x = 0.45, y = 1, angle = 0, hjust = 0, size = 15)
  assign(str_c("p.", gld), p.gld)
}

# Put it all together #
p <- ggdraw() +
  draw_plot(p.All,        x = 0,   y = 0.6833333, width = 0.5, height = 0.3166667) +
  draw_plot(p.Sagebrush,  x = 0.5, y = 0.6833333, width = 0.5, height = 0.3166667) +
  draw_plot(p.Grassland,  x = 0,   y = 0.3666667, width = 0.5, height = 0.3166667) +
  draw_plot(p.Generalist, x = 0.5, y = 0.3666667, width = 0.5, height = 0.3166667) +
  draw_plot(p.Montane,    x = 0,   y = 0.05,      width = 0.5, height = 0.3166667) +
  draw_plot(p.Riparian,   x = 0.5, y = 0.05,      width = 0.5, height = 0.3166667) +
  draw_plot_label("Year", x = 0.5, y = 0.05, angle = 0, hjust = 0)

save_plot("Figure_guild_trends_supported.jpg", p, ncol = 2.5, nrow = 4.5, dpi = 600)

# p <- ggdraw() +
#   draw_plot(p.All,        x = 0, y = 0.8944444, width = 1, height = 0.1055556) +
#   draw_plot(p.Sagebrush,  x = 0, y = 0.7888889, width = 1, height = 0.1055556) +
#   draw_plot(p.Shrubland,  x = 0, y = 0.6833333, width = 1, height = 0.1055556) +
#   draw_plot(p.Grassland,  x = 0, y = 0.5777778, width = 1, height = 0.1055556) +
#   draw_plot(p.Generalist, x = 0, y = 0.4722222, width = 1, height = 0.1055556) +
#   draw_plot(p.Montane,    x = 0, y = 0.3666667, width = 1, height = 0.1055556) +
#   draw_plot(p.Riparian,   x = 0, y = 0.2611111, width = 1, height = 0.1055556) +
#   draw_plot(p.Wetland,    x = 0, y = 0.1555556, width = 1, height = 0.1055556) +
#   draw_plot(p.Woodland,   x = 0, y = 0.05,      width = 1, height = 0.1055556) +
#   draw_plot_label("Year", x = 0.5, y = 0.05, angle = 0, hjust = 0)
# 
# save_plot("Figure_guild_trends_all.jpg", p, ncol = 2, nrow = 8, dpi = 600)
