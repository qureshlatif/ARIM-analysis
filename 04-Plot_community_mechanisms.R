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

#__________ Table compilation function __________#
plot.table.fn <- function(w, B0, B1, X.B, D0, D1, X.D,
                          b0, b1, X.b, X.b.ind,
                          d0, d1, X.d, X.d.ind, X.trend, years,
                          x.min, x.max, xmn.label, xmx.label) {
  # Setup #
  out.plot <- expand.grid(Level = c(xmx.label, xmn.label), Year = years,
                          stringsAsFactors = F) %>%
    mutate(SR.pred = 0,
           SR.pred.lo = 0,
           SR.pred.hi = 0)
  nsp <- dim(B0)[[2]]
  nyr <- length(years)
  
  BETA <- B0
  DELTA <- D0
  for(sp in 1:nsp) {
    BETA[,sp] <- BETA[,sp] + apply(B1[,sp,] * X.B, 1, sum)
    DELTA[,sp] <- DELTA[,sp] + apply(D1[,sp,] * X.D, 1, sum)
  }
  PSI.pred <- array(NA, dim = c(nsims, nsp, nyr))
  for(t in 1:nyr) {
    PSI.pred[,,t] <- QSLpersonal::expit(BETA + DELTA * X.trend[t]) * w
  }
  X.b.xmx <- X.b.xmn <- X.b
  X.d.xmx <- X.d.xmn <- X.d
  X.b.xmn[,X.b.ind] <- X.d.xmn[,X.d.ind] <- x.min
  X.b.xmx[,X.b.ind] <- X.d.xmx[,X.d.ind] <- x.max
  beta_xmx <- beta_xmn <- b0
  delta_xmx <- delta_xmn <- d0
  for(sp in 1:nsp) {
    beta_xmn[,sp] <- beta_xmn[,sp] + apply(b1[,sp,] * X.b.xmn, 1, sum)
    beta_xmx[,sp] <- beta_xmx[,sp] + apply(b1[,sp,] * X.b.xmx, 1, sum)
    delta_xmn[,sp] <- delta_xmn[,sp] + apply(d1[,sp,] * X.d.xmn, 1, sum)
    delta_xmx[,sp] <- delta_xmx[,sp] + apply(d1[,sp,] * X.d.xmx, 1, sum)
  }
  psi_xmn.pred <- psi_xmx.pred <- array(NA, dim = c(nsims, nsp, nyr))
  for(t in 1:nyr) {
    psi_xmn.pred[,,t] <- QSLpersonal::expit(beta_xmn + delta_xmn * X.trend[t])
    psi_xmx.pred[,,t] <- QSLpersonal::expit(beta_xmx + delta_xmx * X.trend[t])
  }
  out.plot$SR.pred[which(out.plot$Level == xmn.label)] <-
    (PSI.pred * psi_xmn.pred) %>% apply(c(1, 3), sum) %>% apply(2, median)
  out.plot$SR.pred.lo[which(out.plot$Level == xmn.label)] <-
    (PSI.pred * psi_xmn.pred) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.1, type = 8)
  out.plot$SR.pred.hi[which(out.plot$Level == xmn.label)] <-
    (PSI.pred * psi_xmn.pred) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.9, type = 8)
  
  out.plot$SR.pred[which(out.plot$Level == xmx.label)] <-
    (PSI.pred * psi_xmx.pred) %>% apply(c(1, 3), sum) %>% apply(2, median)
  out.plot$SR.pred.lo[which(out.plot$Level == xmx.label)] <-
    (PSI.pred * psi_xmx.pred) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.1, type = 8)
  out.plot$SR.pred.hi[which(out.plot$Level == xmx.label)] <-
    (PSI.pred * psi_xmx.pred) %>% apply(c(1, 3), sum) %>% apply(2, quantile, prob = 0.9, type = 8)
  
  return(out.plot)
}
#_____________________________________________#

## All species, Well_1km ##
gld <- "All"
gld.nam <- str_c(gld, " species")
ind.spp <- which(guild.mem[,gld])

var.nam <- "Well_1km"
X.b <- X.psi.pred.hi
X.b.ind <- which(dimnames(X.psi.pred.hi)[[2]] == var.nam)
X.d <- X.lambda.pred.hi
X.d.ind <- which(dimnames(X.lambda.pred.hi)[[2]] == var.nam)

x.min <- min(X.psi[,X.b.ind])
x.max <- max(X.psi[,X.b.ind])
xmn.label <- "0 wells/km^2"
xmx.label <- "5 wells/km^2"

dat.plot <- plot.table.fn(w = mod$mcmcOutput$w[,ind.spp],
                          B0 = mod$mcmcOutput$BETA0[,ind.spp],
                          B1 = mod$mcmcOutput$BETA1[,ind.spp,],
                          X.B = X.PSI.pred.hi,
                          D0 = mod$mcmcOutput$DELTA0[,ind.spp],
                          D1 = mod$mcmcOutput$DELTA1[,ind.spp,],
                          X.D = X.LAMBDA.pred.hi,
                          b0 = mod$mcmcOutput$beta0[,ind.spp],
                          b1 = mod$mcmcOutput$beta1[,ind.spp,],
                          X.b = X.b,
                          X.b.ind = X.b.ind,
                          d0 = mod$mcmcOutput$delta0[,ind.spp],
                          d1 = mod$mcmcOutput$delta1[,ind.spp,],
                          X.d = X.d,
                          X.d.ind = X.d.ind,
                          years = years,
                          X.trend = X.trend,
                          x.min = x.min, x.max = x.max,
                          xmn.label = xmn.label, xmx.label = xmx.label) %>%
  mutate(Level = factor(Level, levels = c(xmn.label, xmx.label)))

p.All_Well <- ggplot(dat.plot, aes(x = Year, y = SR.pred)) +
  geom_ribbon(aes(ymin = SR.pred.lo, ymax = SR.pred.hi, fill = Level),
              size = 0, alpha = 0.3) +
  geom_line(aes(color = Level), size = 1) +
  scale_color_manual(values = c("#000000", "#D55E00"), labels = c(expression(0~wells/km^2), expression(5~wells/km^2))) +
  scale_fill_manual(values = c("#000000", "#D55E00"), labels = c(expression(0~wells/km^2), expression(5~wells/km^2))) +
  ylim(0, max(dat.plot$SR.pred.hi)) +
  scale_x_continuous(breaks = seq(2010, 2019, by = 2)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  guides(color=guide_legend(title = NULL),
         fill=guide_legend(title = NULL)) +
  ggtitle(gld.nam) +
  xlab(NULL) + ylab("Fine-scale richness")

## Generalists, Well_1km ##
gld <- "Generalist"
gld.nam <- str_c(gld, " species")
ind.spp <- which(guild.mem[,gld])

var.nam <- "Well_1km"
X.b <- X.psi.pred.hi
X.b.ind <- which(dimnames(X.psi.pred.hi)[[2]] == var.nam)
X.d <- X.lambda.pred.hi
X.d.ind <- which(dimnames(X.lambda.pred.hi)[[2]] == var.nam)

x.min <- min(X.psi[,X.b.ind])
x.max <- max(X.psi[,X.b.ind])
xmn.label <- "0 wells per sqr km"
xmx.label <- "5 wells per sqr km"

dat.plot <- plot.table.fn(w = mod$mcmcOutput$w[,ind.spp],
                          B0 = mod$mcmcOutput$BETA0[,ind.spp],
                          B1 = mod$mcmcOutput$BETA1[,ind.spp,],
                          X.B = X.PSI.pred.hi,
                          D0 = mod$mcmcOutput$DELTA0[,ind.spp],
                          D1 = mod$mcmcOutput$DELTA1[,ind.spp,],
                          X.D = X.LAMBDA.pred.hi,
                          b0 = mod$mcmcOutput$beta0[,ind.spp],
                          b1 = mod$mcmcOutput$beta1[,ind.spp,],
                          X.b = X.b,
                          X.b.ind = X.b.ind,
                          d0 = mod$mcmcOutput$delta0[,ind.spp],
                          d1 = mod$mcmcOutput$delta1[,ind.spp,],
                          X.d = X.d,
                          X.d.ind = X.d.ind,
                          years = years,
                          X.trend = X.trend,
                          x.min = x.min, x.max = x.max,
                          xmn.label = xmn.label, xmx.label = xmx.label) %>%
  mutate(Level = factor(Level, levels = c(xmn.label, xmx.label)))

p.Gen_Well <- ggplot(dat.plot, aes(x = Year, y = SR.pred)) +
  geom_ribbon(aes(ymin = SR.pred.lo, ymax = SR.pred.hi, fill = Level),
              size = 0, alpha = 0.3) +
  geom_line(aes(color = Level), size = 1) +
  scale_color_manual(values = c("#000000", "#D55E00"), labels = c(expression(0~wells/km^2), expression(5~wells/km^2))) +
  scale_fill_manual(values = c("#000000", "#D55E00"), labels = c(expression(0~wells/km^2), expression(5~wells/km^2))) +
  scale_x_continuous(breaks = seq(2010, 2019, by = 2)) +
  ylim(0, max(dat.plot$SR.pred.hi)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  guides(color=guide_legend(title = NULL),
         fill=guide_legend(title = NULL)) +
  ggtitle(gld.nam) +
  xlab(NULL) + ylab("Fine-scale richness")

## Sagebrush, Well_1km ##
gld <- "Sagebrush"
gld.nam <- str_c(gld, " species")
ind.spp <- which(guild.mem[,gld])

var.nam <- "Well_1km"
X.b <- X.psi.pred.hi
X.b.ind <- which(dimnames(X.psi.pred.hi)[[2]] == var.nam)
X.d <- X.lambda.pred.hi
X.d.ind <- which(dimnames(X.lambda.pred.hi)[[2]] == var.nam)

x.min <- min(X.psi[,X.b.ind])
x.max <- max(X.psi[,X.b.ind])
xmn.label <- "0 wells per sqr km"
xmx.label <- "5 wells per sqr km"

dat.plot <- plot.table.fn(w = mod$mcmcOutput$w[,ind.spp],
                          B0 = mod$mcmcOutput$BETA0[,ind.spp],
                          B1 = mod$mcmcOutput$BETA1[,ind.spp,],
                          X.B = X.PSI.pred.hi,
                          D0 = mod$mcmcOutput$DELTA0[,ind.spp],
                          D1 = mod$mcmcOutput$DELTA1[,ind.spp,],
                          X.D = X.LAMBDA.pred.hi,
                          b0 = mod$mcmcOutput$beta0[,ind.spp],
                          b1 = mod$mcmcOutput$beta1[,ind.spp,],
                          X.b = X.b,
                          X.b.ind = X.b.ind,
                          d0 = mod$mcmcOutput$delta0[,ind.spp],
                          d1 = mod$mcmcOutput$delta1[,ind.spp,],
                          X.d = X.d,
                          X.d.ind = X.d.ind,
                          years = years,
                          X.trend = X.trend,
                          x.min = x.min, x.max = x.max,
                          xmn.label = xmn.label, xmx.label = xmx.label) %>%
  mutate(Level = factor(Level, levels = c(xmn.label, xmx.label)))

p.Sage_Well <- ggplot(dat.plot, aes(x = Year, y = SR.pred)) +
  geom_ribbon(aes(ymin = SR.pred.lo, ymax = SR.pred.hi, fill = Level),
              size = 0, alpha = 0.3) +
  geom_line(aes(color = Level), size = 1) +
  scale_color_manual(values = c("#000000", "#D55E00"), labels = c(expression(0~wells/km^2), expression(5~wells/km^2))) +
  scale_fill_manual(values = c("#000000", "#D55E00"), labels = c(expression(0~wells/km^2), expression(5~wells/km^2))) +
  scale_x_continuous(breaks = seq(2010, 2019, by = 2)) +
  ylim(0, max(dat.plot$SR.pred.hi)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  guides(color=guide_legend(title = NULL),
         fill=guide_legend(title = NULL)) +
  ggtitle(gld.nam) +
  xlab(NULL) + ylab("Fine-scale richness")

## Montane, Well_1km ##
gld <- "Montane"
gld.nam <- str_c(gld, " species")
ind.spp <- which(guild.mem[,gld])

var.nam <- "Well_1km"
X.b <- X.psi.pred.hi
X.b.ind <- which(dimnames(X.psi.pred.hi)[[2]] == var.nam)
X.d <- X.lambda.pred.hi
X.d.ind <- which(dimnames(X.lambda.pred.hi)[[2]] == var.nam)

x.min <- min(X.psi[,X.b.ind])
x.max <- max(X.psi[,X.b.ind])
xmn.label <- "0 wells per sqr km"
xmx.label <- "5 wells per sqr km"

dat.plot <- plot.table.fn(w = mod$mcmcOutput$w[,ind.spp],
                          B0 = mod$mcmcOutput$BETA0[,ind.spp],
                          B1 = mod$mcmcOutput$BETA1[,ind.spp,],
                          X.B = X.PSI.pred.hi,
                          D0 = mod$mcmcOutput$DELTA0[,ind.spp],
                          D1 = mod$mcmcOutput$DELTA1[,ind.spp,],
                          X.D = X.LAMBDA.pred.hi,
                          b0 = mod$mcmcOutput$beta0[,ind.spp],
                          b1 = mod$mcmcOutput$beta1[,ind.spp,],
                          X.b = X.b,
                          X.b.ind = X.b.ind,
                          d0 = mod$mcmcOutput$delta0[,ind.spp],
                          d1 = mod$mcmcOutput$delta1[,ind.spp,],
                          X.d = X.d,
                          X.d.ind = X.d.ind,
                          years = years,
                          X.trend = X.trend,
                          x.min = x.min, x.max = x.max,
                          xmn.label = xmn.label, xmx.label = xmx.label) %>%
  mutate(Level = factor(Level, levels = c(xmn.label, xmx.label)))

p.Mont_Well <- ggplot(dat.plot, aes(x = Year, y = SR.pred)) +
  geom_ribbon(aes(ymin = SR.pred.lo, ymax = SR.pred.hi, fill = Level),
              size = 0, alpha = 0.3) +
  geom_line(aes(color = Level), size = 1) +
  scale_color_manual(values = c("#000000", "#D55E00"), labels = c(expression(0~wells/km^2), expression(5~wells/km^2))) +
  scale_fill_manual(values = c("#000000", "#D55E00"), labels = c(expression(0~wells/km^2), expression(5~wells/km^2))) +
  scale_x_continuous(breaks = seq(2010, 2019, by = 2)) +
  ylim(0, max(dat.plot$SR.pred.hi)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  guides(color=guide_legend(title = NULL),
         fill=guide_legend(title = NULL)) +
  ggtitle(gld.nam) +
  xlab(NULL) + ylab("Fine-scale richness")

## Riparian, Well_1km ##
gld <- "Riparian"
gld.nam <- str_c(gld, " species")
ind.spp <- which(guild.mem[,gld])

var.nam <- "Well_1km"
X.b <- X.psi.pred.hi
X.b.ind <- which(dimnames(X.psi.pred.hi)[[2]] == var.nam)
X.d <- X.lambda.pred.hi
X.d.ind <- which(dimnames(X.lambda.pred.hi)[[2]] == var.nam)

x.min <- min(X.psi[,X.b.ind])
x.max <- max(X.psi[,X.b.ind])
xmn.label <- "0 wells per sqr km"
xmx.label <- "5 wells per sqr km"

dat.plot <- plot.table.fn(w = mod$mcmcOutput$w[,ind.spp],
                          B0 = mod$mcmcOutput$BETA0[,ind.spp],
                          B1 = mod$mcmcOutput$BETA1[,ind.spp,],
                          X.B = X.PSI.pred.hi,
                          D0 = mod$mcmcOutput$DELTA0[,ind.spp],
                          D1 = mod$mcmcOutput$DELTA1[,ind.spp,],
                          X.D = X.LAMBDA.pred.hi,
                          b0 = mod$mcmcOutput$beta0[,ind.spp],
                          b1 = mod$mcmcOutput$beta1[,ind.spp,],
                          X.b = X.b,
                          X.b.ind = X.b.ind,
                          d0 = mod$mcmcOutput$delta0[,ind.spp],
                          d1 = mod$mcmcOutput$delta1[,ind.spp,],
                          X.d = X.d,
                          X.d.ind = X.d.ind,
                          years = years,
                          X.trend = X.trend,
                          x.min = x.min, x.max = x.max,
                          xmn.label = xmn.label, xmx.label = xmx.label) %>%
  mutate(Level = factor(Level, levels = c(xmn.label, xmx.label)))

p.Ripar_Well <- ggplot(dat.plot, aes(x = Year, y = SR.pred)) +
  geom_ribbon(aes(ymin = SR.pred.lo, ymax = SR.pred.hi, fill = Level),
              size = 0, alpha = 0.3) +
  geom_line(aes(color = Level), size = 1) +
  scale_color_manual(values = c("#000000", "#D55E00"), labels = c(expression(0~wells/km^2), expression(5~wells/km^2))) +
  scale_fill_manual(values = c("#000000", "#D55E00"), labels = c(expression(0~wells/km^2), expression(5~wells/km^2))) +
  scale_x_continuous(breaks = seq(2010, 2019, by = 2)) +
  ylim(0, max(dat.plot$SR.pred.hi)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  guides(color=guide_legend(title = NULL),
         fill=guide_legend(title = NULL)) +
  ggtitle(gld.nam) +
  xlab(NULL) + ylab("Fine-scale richness")

## Grassland, Annual herbaceous cover ##
gld <- "Grassland"
gld.nam <- str_c(gld, " species")
ind.spp <- which(guild.mem[,gld])

var.nam <- "AHerb"
X.b <- X.psi.pred.hi
X.b.ind <- which(dimnames(X.psi.pred.hi)[[2]] == var.nam)
X.d <- X.lambda.pred.hi
X.d.ind <- which(dimnames(X.lambda.pred.hi)[[2]] == var.nam)

x.min <- min(X.psi[,X.b.ind])
x.max <- max(X.psi[,X.b.ind])
xmn.label <- "0% cheatgrass"
xmx.label <- "17% cheatgrass"

dat.plot <- plot.table.fn(w = mod$mcmcOutput$w[,ind.spp],
                          B0 = mod$mcmcOutput$BETA0[,ind.spp],
                          B1 = mod$mcmcOutput$BETA1[,ind.spp,],
                          X.B = X.PSI.pred.hi,
                          D0 = mod$mcmcOutput$DELTA0[,ind.spp],
                          D1 = mod$mcmcOutput$DELTA1[,ind.spp,],
                          X.D = X.LAMBDA.pred.hi,
                          b0 = mod$mcmcOutput$beta0[,ind.spp],
                          b1 = mod$mcmcOutput$beta1[,ind.spp,],
                          X.b = X.b,
                          X.b.ind = X.b.ind,
                          d0 = mod$mcmcOutput$delta0[,ind.spp],
                          d1 = mod$mcmcOutput$delta1[,ind.spp,],
                          X.d = X.d,
                          X.d.ind = X.d.ind,
                          years = years,
                          X.trend = X.trend,
                          x.min = x.min, x.max = x.max,
                          xmn.label = xmn.label, xmx.label = xmx.label) %>%
  mutate(Level = factor(Level, levels = c(xmn.label, xmx.label)))

p.Grass_AHerb <- ggplot(dat.plot, aes(x = Year, y = SR.pred)) +
  geom_ribbon(aes(ymin = SR.pred.lo, ymax = SR.pred.hi, fill = Level),
              size = 0, alpha = 0.3) +
  geom_line(aes(color = Level), size = 1) +
  scale_color_manual(values = c("#000000", "#D55E00")) +
  scale_fill_manual(values = c("#000000", "#D55E00")) +
  scale_x_continuous(breaks = seq(2010, 2019, by = 2)) +
  ylim(0, max(dat.plot$SR.pred.hi)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  guides(color=guide_legend(title = NULL),
         fill=guide_legend(title = NULL)) +
  ggtitle(gld.nam) +
  xlab(NULL) + ylab("Fine-scale richness")

# Put it all together #
p <- ggdraw() +
  draw_plot(p.All_Well,    x = 0,   y = 0.6833333, width = 0.5, height = 0.3166667) +
  draw_plot(p.Gen_Well,    x = 0.5, y = 0.6833333, width = 0.5, height = 0.3166667) +
  draw_plot(p.Sage_Well,   x = 0,   y = 0.3666667, width = 0.5, height = 0.3166667) +
  draw_plot(p.Grass_AHerb, x = 0.5, y = 0.3666667, width = 0.5, height = 0.3166667) +
  draw_plot(p.Mont_Well,   x = 0,   y = 0.05,      width = 0.5, height = 0.3166667) +
  draw_plot(p.Ripar_Well,  x = 0.5, y = 0.05,      width = 0.5, height = 0.3166667) +
  draw_plot_label("Year",  x = 0.5, y = 0.05, angle = 0, hjust = 0)

save_plot("Figure_guild_mechanisms.jpg", p, ncol = 1.5, nrow = 3, dpi = 600)

