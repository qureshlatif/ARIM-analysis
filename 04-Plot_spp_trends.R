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

plot.table.fn <- function(B0, B1, X.B.hi, X.B.lo, X.B.bg, dev.B,
                          D0, D1, X.D.hi, X.D.lo, X.D.bg, X.trend) {
  # Setup #
  out.plot <- expand.grid(Development = c("High", "Low", "Background"), Year = 2010:2019, stringsAsFactors = F) %>%
    mutate(PSI.yr = 0,
           PSI.yr.lo = 0,
           PSI.yr.hi = 0,
           PSI.pred = 0,
           PSI.pred.lo = 0,
           PSI.pred.hi = 0)
  nyr <- ncol(dev.B)

  # High development #
  BETA <- B0 + apply(B1 * X.B.hi, 1, sum)
  DELTA <- D0 + apply(D1 * X.D.hi, 1, sum)
  PSI.pred <- PSI.yr <- matrix(NA, nrow = nsims, ncol = nyr)
  for(t in 1:nyr) {
    PSI.pred[,t] <- QSLpersonal::expit(BETA + DELTA * X.trend[t])
    PSI.yr[,t] <- QSLpersonal::expit(BETA + dev.B[,t] +
                                       DELTA * X.trend[t])
  }
  out.plot$PSI.yr[which(out.plot$Development == "High")] <-
    apply(PSI.yr, 2, median)
  out.plot$PSI.yr.lo[which(out.plot$Development == "High")] <-
    apply(PSI.yr, 2, quantile, prob = 0.05, type = 8)
  out.plot$PSI.yr.hi[which(out.plot$Development == "High")] <-
    apply(PSI.yr, 2, quantile, prob = 0.95, type = 8)
  out.plot$PSI.pred[which(out.plot$Development == "High")] <-
    apply(PSI.pred, 2, median)
  out.plot$PSI.pred.lo[which(out.plot$Development == "High")] <-
    apply(PSI.pred, 2, quantile, prob = 0.05, type = 8)
  out.plot$PSI.pred.hi[which(out.plot$Development == "High")] <-
    apply(PSI.pred, 2, quantile, prob = 0.95, type = 8)

  # Low development #
  BETA <- B0 + apply(B1 * X.B.lo, 1, sum)
  DELTA <- D0 + apply(D1 * X.D.lo, 1, sum)
  PSI.pred <- PSI.yr <- matrix(NA, nrow = nsims, ncol = nyr)
  for(t in 1:nyr) {
    PSI.pred[,t] <- QSLpersonal::expit(BETA + DELTA * X.trend[t])
    PSI.yr[,t] <- QSLpersonal::expit(BETA + dev.B[,t] +
                                       DELTA * X.trend[t])
  }
  out.plot$PSI.yr[which(out.plot$Development == "Low")] <-
    apply(PSI.yr, 2, median)
  out.plot$PSI.yr.lo[which(out.plot$Development == "Low")] <-
    apply(PSI.yr, 2, quantile, prob = 0.05, type = 8)
  out.plot$PSI.yr.hi[which(out.plot$Development == "Low")] <-
    apply(PSI.yr, 2, quantile, prob = 0.95, type = 8)
  out.plot$PSI.pred[which(out.plot$Development == "Low")] <-
    apply(PSI.pred, 2, median)
  out.plot$PSI.pred.lo[which(out.plot$Development == "Low")] <-
    apply(PSI.pred, 2, quantile, prob = 0.05, type = 8)
  out.plot$PSI.pred.hi[which(out.plot$Development == "Low")] <-
    apply(PSI.pred, 2, quantile, prob = 0.95, type = 8)
  
  # Background #
  BETA <- B0 + apply(B1 * X.B.bg, 1, sum)
  DELTA <- D0 + apply(D1 * X.D.bg, 1, sum)
  PSI.pred <- PSI.yr <- matrix(NA, nrow = nsims, ncol = nyr)
  for(t in 1:nyr) {
    PSI.pred[,t] <- QSLpersonal::expit(BETA + DELTA * X.trend[t])
    PSI.yr[,t] <- QSLpersonal::expit(BETA + dev.B[,t] +
                                       DELTA * X.trend[t])
  }
  out.plot$PSI.yr[which(out.plot$Development == "Background")] <-
    apply(PSI.yr, 2, median)
  out.plot$PSI.yr.lo[which(out.plot$Development == "Background")] <-
    apply(PSI.yr, 2, quantile, prob = 0.05, type = 8)
  out.plot$PSI.yr.hi[which(out.plot$Development == "Background")] <-
    apply(PSI.yr, 2, quantile, prob = 0.95, type = 8)
  out.plot$PSI.pred[which(out.plot$Development == "Background")] <-
    apply(PSI.pred, 2, median)
  out.plot$PSI.pred.lo[which(out.plot$Development == "Background")] <-
    apply(PSI.pred, 2, quantile, prob = 0.05, type = 8)
  out.plot$PSI.pred.hi[which(out.plot$Development == "Background")] <-
    apply(PSI.pred, 2, quantile, prob = 0.95, type = 8)
  
  return(out.plot)
}

spp.plot <- R.utils::loadObject("Spp_plot_trends")
spp.plot <- c(spp.plot, "GTTO") # Add Green-tailed Towhee
spp.names <- read.csv("C:/Users/Quresh.Latif/files/data/Alpha_codes_tax_20220324.csv", stringsAsFactors = F) %>%
  filter(SPEC %in% spp.plot) %>%
  pull(COMMONNAME)
names(spp.names) <- read.csv("C:/Users/Quresh.Latif/files/data/Alpha_codes_tax_20220324.csv", stringsAsFactors = F) %>%
  filter(SPEC %in% spp.plot) %>%
  pull(SPEC)
spp.names <- spp.names[spp.plot]

for(sp in 1:length(spp.plot)) {
  spp.name <- spp.names[sp]
  spp <- spp.plot[sp]
  ind.spp <- which(spp.list == spp)
  
  # Grid cell occupancy #
  dat.plot <- plot.table.fn(B0 = mod$mcmcOutput$BETA0[,ind.spp],
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
  
  p.grid <- ggplot(dat.plot, aes(x = Year, y = PSI.pred)) +
    geom_ribbon(aes(ymin = PSI.pred.lo, ymax = PSI.pred.hi, fill = Development),
                size = 0, alpha = 0.3) +
    geom_line(aes(color = Development), size = 1) +
    geom_errorbar(aes(x = Year.jitter, ymin = PSI.yr.lo, ymax = PSI.yr.hi,
                      color = Development), width = 0.1) +
    geom_point(aes(x = Year.jitter, y = PSI.yr,
                   color = Development, shape = Development)) +
    scale_color_manual(values = c("#D55E00", "#009E73", "#000000")) +
    scale_fill_manual(values = c("#D55E00", "#009E73", "#000000")) +
    scale_shape_manual(values = c(15, 17, 16)) +
    scale_x_continuous(breaks = seq(2010, 2019, by = 2)) +
    ylim(0, 1)
  if(spp %in% c("HOLA", "BRSP")) {
    p.grid <- p.grid +
      theme(legend.position = c(1,0), legend.justification = c(1,0)) +
      xlab(NULL) + ylab("Grid-cell occupancy")  
  } else {
    p.grid <- p.grid +
      guides(color = F, fill = F, shape = F) +
      xlab(NULL) + ylab("Grid-cell occupancy")  
  }

  # Point occupancy #
  dat.plot <- plot.table.fn(B0 = mod$mcmcOutput$beta0[,ind.spp],
                            B1 = mod$mcmcOutput$beta1[,ind.spp,],
                            X.B.hi = X.psi.pred.hi,
                            X.B.lo = X.psi.pred.lo,
                            X.B.bg = X.psi.pred.bg,
                            dev.B = mod$mcmcOutput$dev.beta[,ind.spp,],
                            D0 = mod$mcmcOutput$delta0[,ind.spp],
                            D1 = mod$mcmcOutput$delta1[,ind.spp,],
                            X.D.hi = X.lambda.pred.hi,
                            X.D.lo = X.lambda.pred.lo,
                            X.D.bg = X.lambda.pred.bg,
                            X.trend = X.trend) %>%
    mutate(Development = factor(Development, levels = c("High", "Low", "Background")),
           Year.jitter = ifelse(Development == "High", Year + 0.1,
                                ifelse(Development == "Background", Year - 0.1,
                                       Year)))
  
  p.point <- ggplot(dat.plot, aes(x = Year, y = PSI.pred)) +
    geom_ribbon(aes(ymin = PSI.pred.lo, ymax = PSI.pred.hi, fill = Development),
                size = 0, alpha = 0.3) +
    geom_line(aes(color = Development), size = 1) +
    geom_errorbar(aes(x = Year.jitter, ymin = PSI.yr.lo, ymax = PSI.yr.hi,
                      color = Development), width = 0.1) +
    geom_point(aes(x = Year.jitter, y = PSI.yr,
                   color = Development, shape = Development)) +
    scale_color_manual(values = c("#D55E00", "#009E73", "#000000")) +
    scale_fill_manual(values = c("#D55E00", "#009E73", "#000000")) +
    scale_shape_manual(values = c(15, 17, 16)) +
    scale_x_continuous(breaks = seq(2010, 2019, by = 2)) +
    ylim(0, 1) +
    #theme(legend.position = c(1,0), legend.justification = c(1,0)) +
    guides(color = F, fill = F, shape = F) +
    xlab(NULL) + ylab("Point occupancy")
  
  # Put everything together
  p.spp <- ggdraw() +
    draw_plot(p.grid,  x = 0,   y = 0, width = 0.5, height = 0.95) +
    draw_plot(p.point, x = 0.5, y = 0, width = 0.5, height = 0.95) +
    draw_plot_label(spp.name, x = 0.45, y = 1, angle = 0, hjust = 0, size = 15)
  assign(str_c("p.", spp), p.spp)
}

p <- ggdraw() +
  draw_plot(p.BRSP, x = 0, y = 0.7625, width = 1, height = 0.2375) +
  draw_plot(p.SABS, x = 0, y = 0.5250, width = 1, height = 0.2375) +
  draw_plot(p.GTTO, x = 0, y = 0.2875, width = 1, height = 0.2375) +
  draw_plot(p.SATH, x = 0, y = 0.05,   width = 1, height = 0.2375) +
  draw_plot_label("Year", x = 0.5, y = 0.05, angle = 0, hjust = 0)

save_plot("Figure_focal_spp_trends.jpg", p, ncol = 1.5, nrow = 3.5, dpi = 600)

p <- ggdraw() +
  draw_plot(p.HOLA, x = 0,   y = 0.8416667, width = 0.5, height = 0.1583333) +
  draw_plot(p.WEME, x = 0,   y = 0.6833333, width = 0.5, height = 0.1583333) +
  draw_plot(p.CONI, x = 0,   y = 0.5250000, width = 0.5, height = 0.1583333) +
  draw_plot(p.KILL, x = 0,   y = 0.3666667, width = 0.5, height = 0.1583333) +
  draw_plot(p.BRBL, x = 0,   y = 0.2083333, width = 0.5, height = 0.1583333) +
  draw_plot(p.BHCO, x = 0,   y = 0.05,      width = 0.5, height = 0.1583333) +
  draw_plot(p.HOWR, x = 0.5, y = 0.8416667, width = 0.5, height = 0.1583333) +
  draw_plot(p.COGR, x = 0.5, y = 0.6833333, width = 0.5, height = 0.1583333) +
  draw_plot(p.CORA, x = 0.5, y = 0.5250000, width = 0.5, height = 0.1583333) +
  draw_plot(p.AMRO, x = 0.5, y = 0.3666667, width = 0.5, height = 0.1583333) +
  draw_plot(p.VGSW, x = 0.5, y = 0.2083333, width = 0.5, height = 0.1583333) +
  draw_plot(p.ROWR, x = 0.5, y = 0.05,      width = 0.5, height = 0.1583333) +
  draw_plot_label("Year", x = 0.5, y = 0.05, angle = 0, hjust = 0)

save_plot("Figure_other_spp_neg_trends.jpg", p, ncol = 3, nrow = 5, dpi = 600)
