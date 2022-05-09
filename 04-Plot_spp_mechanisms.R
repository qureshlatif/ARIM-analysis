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

plot.table.fn <- function(B0, B1, X.B, X.B.ind,
                          D0, D1, X.D, X.D.ind, X.trend, years,
                          x.max, x.min,
                          xmx.label, xmn.label) {
  # Setup #
  out.plot <- expand.grid(Level = c(xmx.label, xmn.label), Year = years, stringsAsFactors = F) %>%
    mutate(PSI.pred = 0,
           PSI.pred.lo = 0,
           PSI.pred.hi = 0)
  nyr <- length(years)
  
  X.B.xmx <- X.B.xmn <- X.B
  X.D.xmx <- X.D.xmn <- X.D
  X.B.xmx[,X.B.ind] <- X.D.xmx[,X.D.ind] <- x.max
  X.B.xmn[,X.B.ind] <- X.D.xmn[,X.D.ind] <- x.min
  BETA_xmx <- B0 + apply(B1 * X.B.xmx, 1, sum)
  DELTA_xmx <- D0 + apply(D1 * X.D.xmx, 1, sum)
  BETA_xmn <- B0 + apply(B1 * X.B.xmn, 1, sum)
  DELTA_xmn <- D0 + apply(D1 * X.D.xmn, 1, sum)
  PSI_xmx.pred <- PSI_xmn.pred <-
    matrix(NA, nrow = nsims, ncol = nyr)
  for(t in 1:nyr) {
    PSI_xmx.pred[,t] <- QSLpersonal::expit(BETA_xmx + DELTA_xmx * X.trend[t])
    PSI_xmn.pred[,t] <- QSLpersonal::expit(BETA_xmn + DELTA_xmn * X.trend[t])
  }
  out.plot$PSI.pred[which(out.plot$Level == xmx.label)] <-
    apply(PSI_xmx.pred, 2, median)
  out.plot$PSI.pred.lo[which(out.plot$Level == xmx.label)] <-
    apply(PSI_xmx.pred, 2, quantile, prob = 0.1, type = 8)
  out.plot$PSI.pred.hi[which(out.plot$Level == xmx.label)] <-
    apply(PSI_xmx.pred, 2, quantile, prob = 0.9, type = 8)
  out.plot$PSI.pred[which(out.plot$Level == xmn.label)] <-
    apply(PSI_xmn.pred, 2, median)
  out.plot$PSI.pred.lo[which(out.plot$Level == xmn.label)] <-
    apply(PSI_xmn.pred, 2, quantile, prob = 0.1, type = 8)
  out.plot$PSI.pred.hi[which(out.plot$Level == xmn.label)] <-
    apply(PSI_xmn.pred, 2, quantile, prob = 0.9, type = 8)

  return(out.plot)
}

## Sagebrush sparrow, Well_1km ##

spp.name <- "Sagebrush Sparrow"
spp <- "SABS"
ind.spp <- which(spp.list == spp)

var.nam <- "Well_1km"
X.B <- X.psi.pred.hi
X.B.ind <- which(dimnames(X.psi.pred.hi)[[2]] == var.nam)
X.D <- X.lambda.pred.hi
X.D.ind <- which(dimnames(X.lambda.pred.hi)[[2]] == var.nam)

x.min <- min(X.psi[,X.B.ind])
x.max <- max(X.psi[,X.B.ind])
xmn.label <- "0 wells per sqr km"
xmx.label <- "5 wells per sqr km"

dat.plot <- plot.table.fn(B0 = mod$mcmcOutput$beta0[,ind.spp],
                          B1 = mod$mcmcOutput$beta1[,ind.spp,],
                          X.B = X.B,
                          X.B.ind = X.B.ind,
                          D0 = mod$mcmcOutput$delta0[,ind.spp],
                          D1 = mod$mcmcOutput$delta1[,ind.spp,],
                          X.D = X.D,
                          X.D.ind = X.D.ind,
                          X.trend = X.trend,
                          years = years,
                          x.min = x.min,
                          x.max = x.max,
                          xmn.label = xmn.label,
                          xmx.label = xmx.label
                          ) %>%
  mutate(Level = factor(Level, levels = c(xmn.label, xmx.label)))

p.SABS <- ggplot(dat.plot, aes(x = Year, y = PSI.pred)) +
  geom_ribbon(aes(ymin = PSI.pred.lo, ymax = PSI.pred.hi, fill = Level),
              size = 0, alpha = 0.3) +
  geom_line(aes(color = Level), size = 1) +
  scale_color_manual(values = c("#000000", "#D55E00")) +
  scale_fill_manual(values = c("#000000", "#D55E00")) +
  scale_x_continuous(breaks = seq(2010, 2019, by = 2)) +
  ylim(0, 1) +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  guides(color=guide_legend(title = NULL),
         fill=guide_legend(title = NULL)) +
  ggtitle(spp.name) +
  xlab(NULL) + ylab("Point occupancy")  


## Brewer's sparrow, Well_1km ##

spp.name <- "Brewer's Sparrow"
spp <- "BRSP"
ind.spp <- which(spp.list == spp)

var.nam <- "Well_1km"
X.B <- X.psi.pred.hi
X.B.ind <- which(dimnames(X.psi.pred.hi)[[2]] == var.nam)
X.D <- X.lambda.pred.hi
X.D.ind <- which(dimnames(X.lambda.pred.hi)[[2]] == var.nam)

x.min <- min(X.psi[,X.B.ind])
x.max <- max(X.psi[,X.B.ind])
xmn.label <- "0 wells per sqr km"
xmx.label <- "5 wells per sqr km"

dat.plot <- plot.table.fn(B0 = mod$mcmcOutput$beta0[,ind.spp],
                          B1 = mod$mcmcOutput$beta1[,ind.spp,],
                          X.B = X.B,
                          X.B.ind = X.B.ind,
                          D0 = mod$mcmcOutput$delta0[,ind.spp],
                          D1 = mod$mcmcOutput$delta1[,ind.spp,],
                          X.D = X.D,
                          X.D.ind = X.D.ind,
                          X.trend = X.trend,
                          years = years,
                          x.min = x.min,
                          x.max = x.max,
                          xmn.label = xmn.label,
                          xmx.label = xmx.label
) %>%
  mutate(Level = factor(Level, levels = c(xmn.label, xmx.label)))

p.BRSP <- ggplot(dat.plot, aes(x = Year, y = PSI.pred)) +
  geom_ribbon(aes(ymin = PSI.pred.lo, ymax = PSI.pred.hi, fill = Level),
              size = 0, alpha = 0.3) +
  geom_line(aes(color = Level), size = 1) +
  scale_color_manual(values = c("#000000", "#D55E00")) +
  scale_fill_manual(values = c("#000000", "#D55E00")) +
  scale_x_continuous(breaks = seq(2010, 2019, by = 2)) +
  ylim(0, 1) +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  guides(color=guide_legend(title = NULL),
         fill=guide_legend(title = NULL)) +
  ggtitle(spp.name) +
  xlab(NULL) + ylab("Point occupancy")  

## Horned Lark, Annual herbaceous ##

spp.name <- "Horned Lark"
spp <- "HOLA"
ind.spp <- which(spp.list == spp)

var.nam <- "AHerb"
X.B <- X.psi.pred.hi
X.B.ind <- which(dimnames(X.psi.pred.hi)[[2]] == var.nam)
X.D <- X.lambda.pred.hi
X.D.ind <- which(dimnames(X.lambda.pred.hi)[[2]] == var.nam)

x.min <- min(X.psi[,X.B.ind])
x.max <- max(X.psi[,X.B.ind])
xmn.label <- "0% annual herbaceous cover"
xmx.label <- "17% annuall herbaceous cover"

dat.plot <- plot.table.fn(B0 = mod$mcmcOutput$beta0[,ind.spp],
                          B1 = mod$mcmcOutput$beta1[,ind.spp,],
                          X.B = X.B,
                          X.B.ind = X.B.ind,
                          D0 = mod$mcmcOutput$delta0[,ind.spp],
                          D1 = mod$mcmcOutput$delta1[,ind.spp,],
                          X.D = X.D,
                          X.D.ind = X.D.ind,
                          X.trend = X.trend,
                          years = years,
                          x.min = x.min,
                          x.max = x.max,
                          xmn.label = xmn.label,
                          xmx.label = xmx.label
) %>%
  mutate(Level = factor(Level, levels = c(xmn.label, xmx.label)))

p.HOLA <- ggplot(dat.plot, aes(x = Year, y = PSI.pred)) +
  geom_ribbon(aes(ymin = PSI.pred.lo, ymax = PSI.pred.hi, fill = Level),
              size = 0, alpha = 0.3) +
  geom_line(aes(color = Level), size = 1) +
  scale_color_manual(values = c("#000000", "#D55E00")) +
  scale_fill_manual(values = c("#000000", "#D55E00")) +
  scale_x_continuous(breaks = seq(2010, 2019, by = 2)) +
  ylim(0, 1) +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  guides(color=guide_legend(title = NULL),
         fill=guide_legend(title = NULL)) +
  ggtitle(spp.name) +
  xlab(NULL) + ylab("Point occupancy")  

p <- ggdraw() +
  draw_plot(p.BRSP, x = 0, y = 0.6833333, width = 1, height = 0.3166667) +
  draw_plot(p.SABS, x = 0, y = 0.3666667, width = 1, height = 0.3166667) +
  draw_plot(p.HOLA, x = 0, y = 0.05,      width = 1, height = 0.3166667) +
  draw_plot_label("Year", x = 0.5, y = 0.05, angle = 0, hjust = 0)

save_plot("Figure_spp_mechanisms.jpg", p, ncol = 0.7, nrow = 2.5, dpi = 600)
