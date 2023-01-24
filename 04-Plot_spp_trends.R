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


########################
# Focal species trends #
########################

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
    apply(PSI.yr, 2, quantile, prob = 0.1, type = 8)
  out.plot$PSI.yr.hi[which(out.plot$Development == "High")] <-
    apply(PSI.yr, 2, quantile, prob = 0.9, type = 8)
  out.plot$PSI.pred[which(out.plot$Development == "High")] <-
    apply(PSI.pred, 2, median)
  out.plot$PSI.pred.lo[which(out.plot$Development == "High")] <-
    apply(PSI.pred, 2, quantile, prob = 0.1, type = 8)
  out.plot$PSI.pred.hi[which(out.plot$Development == "High")] <-
    apply(PSI.pred, 2, quantile, prob = 0.9, type = 8)

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
    apply(PSI.yr, 2, quantile, prob = 0.1, type = 8)
  out.plot$PSI.yr.hi[which(out.plot$Development == "Low")] <-
    apply(PSI.yr, 2, quantile, prob = 0.9, type = 8)
  out.plot$PSI.pred[which(out.plot$Development == "Low")] <-
    apply(PSI.pred, 2, median)
  out.plot$PSI.pred.lo[which(out.plot$Development == "Low")] <-
    apply(PSI.pred, 2, quantile, prob = 0.1, type = 8)
  out.plot$PSI.pred.hi[which(out.plot$Development == "Low")] <-
    apply(PSI.pred, 2, quantile, prob = 0.9, type = 8)
  
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
    apply(PSI.yr, 2, quantile, prob = 0.1, type = 8)
  out.plot$PSI.yr.hi[which(out.plot$Development == "Background")] <-
    apply(PSI.yr, 2, quantile, prob = 0.9, type = 8)
  out.plot$PSI.pred[which(out.plot$Development == "Background")] <-
    apply(PSI.pred, 2, median)
  out.plot$PSI.pred.lo[which(out.plot$Development == "Background")] <-
    apply(PSI.pred, 2, quantile, prob = 0.1, type = 8)
  out.plot$PSI.pred.hi[which(out.plot$Development == "Background")] <-
    apply(PSI.pred, 2, quantile, prob = 0.9, type = 8)
  
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
      xlab(NULL) + ylab("Coarse-scale occupancy")  
  } else {
    p.grid <- p.grid +
      guides(color = "none", fill = "none", shape = "none") +
      xlab(NULL) + ylab("Coarse-scale occupancy")  
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
    guides(color = "none", fill = "none", shape = "none") +
    xlab(NULL) + ylab("Fine-scale occupancy")
  
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


######################
# All species trends #
######################

## Set up output table ##
cols <- c("Spp", "TREND_hi", "TREND_hi_lo", "TREND_hi_hi",
          "TREND_lo", "TREND_lo_lo", "TREND_lo_hi", "supp_Thi_lt_Tlo", "supp_Thi_2lt_Tlo",
          "TREND_bg", "TREND_bg_lo", "TREND_bg_hi", "supp_Thi_lt_Tbg", "supp_Thi_2lt_Tbg",
          "trend_hi", "trend_hi_lo", "trend_hi_hi",
          "trend_lo", "trend_lo_lo", "trend_lo_hi", "supp_thi_lt_tlo", "supp_thi_2lt_tlo",
          "trend_bg", "trend_bg_lo", "trend_bg_hi", "supp_thi_lt_tbg", "supp_thi_2lt_tbg")
out <- matrix(NA, nrow = length(spp.list), ncol = length(cols),
              dimnames = list(spp.list, cols)) %>%
  data.frame(stringsAsFactors = F) %>%
  mutate(Spp = as.character(Spp)) %>%
  mutate(Spp = spp.list)

## Stratum covariate values ##
source(str_c(scripts.loc, "Calculate_mech_path_covariate_values.R"))

## Regression parameters ##
BETA_hi <- BETA_lo <- BETA_bg <- mod$mcmcOutput$BETA0
DELTA_hi <- DELTA_lo <- DELTA_bg <- mod$mcmcOutput$DELTA0
beta_hi <- beta_lo <- beta_bg <- mod$mcmcOutput$beta0
delta_hi <- delta_lo <- delta_bg <- mod$mcmcOutput$delta0
for(sp in 1:n.spp) {
  BETA_hi[,sp] <- BETA_hi[,sp] + apply(mod$mcmcOutput$BETA1[,sp,] * X.PSI.pred.hi, 1, sum)
  BETA_lo[,sp] <- BETA_lo[,sp] + apply(mod$mcmcOutput$BETA1[,sp,] * X.PSI.pred.lo, 1, sum)
  BETA_bg[,sp] <- BETA_bg[,sp] + apply(mod$mcmcOutput$BETA1[,sp,] * X.PSI.pred.bg, 1, sum)
  
  DELTA_hi[,sp] <- DELTA_hi[,sp] + apply(mod$mcmcOutput$DELTA1[,sp,] * X.LAMBDA.pred.hi, 1, sum)
  DELTA_lo[,sp] <- DELTA_lo[,sp] + apply(mod$mcmcOutput$DELTA1[,sp,] * X.LAMBDA.pred.lo, 1, sum)
  DELTA_bg[,sp] <- DELTA_bg[,sp] + apply(mod$mcmcOutput$DELTA1[,sp,] * X.LAMBDA.pred.bg, 1, sum)
  
  beta_hi[,sp] <- beta_hi[,sp] + apply(mod$mcmcOutput$beta1[,sp,] * X.psi.pred.hi, 1, sum)
  beta_lo[,sp] <- beta_lo[,sp] + apply(mod$mcmcOutput$beta1[,sp,] * X.psi.pred.lo, 1, sum)
  beta_bg[,sp] <- beta_bg[,sp] + apply(mod$mcmcOutput$beta1[,sp,] * X.psi.pred.bg, 1, sum)
  
  delta_hi[,sp] <- delta_hi[,sp] + apply(mod$mcmcOutput$delta1[,sp,] * X.lambda.pred.hi, 1, sum)
  delta_lo[,sp] <- delta_lo[,sp] + apply(mod$mcmcOutput$delta1[,sp,] * X.lambda.pred.lo, 1, sum)
  delta_bg[,sp] <- delta_bg[,sp] + apply(mod$mcmcOutput$delta1[,sp,] * X.lambda.pred.bg, 1, sum)
}

## Grid-cell trend ##
PSI1 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[1])
PSI10 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[10])
TREND_hi <- log(PSI10 / PSI1)
out$TREND_hi <- apply(TREND_hi, 2, median)
out$TREND_hi_lo <- apply(TREND_hi, 2, function(x) quantile(x, prob = 0.1, type = 8))
out$TREND_hi_hi <- apply(TREND_hi, 2, function(x) quantile(x, prob = 0.9, type = 8))

PSI1 <- QSLpersonal::expit(BETA_lo + DELTA_lo * X.trend[1])
PSI10 <- QSLpersonal::expit(BETA_lo + DELTA_lo * X.trend[10])
TREND_lo <- log(PSI10 / PSI1)
out$TREND_lo <- apply(TREND_lo, 2, median)
out$TREND_lo_lo <- apply(TREND_lo, 2, function(x) quantile(x, prob = 0.1, type = 8))
out$TREND_lo_hi <- apply(TREND_lo, 2, function(x) quantile(x, prob = 0.9, type = 8))
DIFF <- exp(TREND_hi) - exp(TREND_lo)
out$supp_Thi_lt_Tlo <- DIFF %>%
  apply(2, function(x) sum(x < 0)/nsims) %>%
  (function(x) x >= 0.9 | x <= 0.1)
out$supp_Thi_2lt_Tlo <- DIFF %>%
  apply(2, function(x) sum(x < -0.02)/nsims) %>%
  (function(x) x >= 0.9)

PSI1 <- QSLpersonal::expit(BETA_bg + DELTA_bg * X.trend[1])
PSI10 <- QSLpersonal::expit(BETA_bg + DELTA_bg * X.trend[10])
TREND_bg <- log(PSI10 / PSI1)
out$TREND_bg <- apply(TREND_bg, 2, median)
out$TREND_bg_lo <- apply(TREND_bg, 2, function(x) quantile(x, prob = 0.1, type = 8))
out$TREND_bg_hi <- apply(TREND_bg, 2, function(x) quantile(x, prob = 0.9, type = 8))
DIFF <- exp(TREND_hi) - exp(TREND_bg)
out$supp_Thi_lt_Tbg <- DIFF %>%
  apply(2, function(x) sum(x < 0)/nsims) %>%
  (function(x) x >= 0.9 | x <= 0.1)
out$supp_Thi_2lt_Tbg <- DIFF %>%
  apply(2, function(x) sum(x < -0.02)/nsims) %>%
  (function(x) x >= 0.9)

## Point-scale trend ##
psi1 <- QSLpersonal::expit(beta_hi + delta_hi * X.trend[1])
psi10 <- QSLpersonal::expit(beta_hi + delta_hi * X.trend[10])
trend_hi <- log(psi10 / psi1)
out$trend_hi <- apply(trend_hi, 2, median)
out$trend_hi_lo <- apply(trend_hi, 2, function(x) quantile(x, prob = 0.1, type = 8))
out$trend_hi_hi <- apply(trend_hi, 2, function(x) quantile(x, prob = 0.9, type = 8))

psi1 <- QSLpersonal::expit(beta_lo + delta_lo * X.trend[1])
psi10 <- QSLpersonal::expit(beta_lo + delta_lo * X.trend[10])
trend_lo <- log(psi10 / psi1)
out$trend_lo <- apply(trend_lo, 2, median)
out$trend_lo_lo <- apply(trend_lo, 2, function(x) quantile(x, prob = 0.1, type = 8))
out$trend_lo_hi <- apply(trend_lo, 2, function(x) quantile(x, prob = 0.9, type = 8))
DIFF <- exp(trend_hi) - exp(trend_lo)
out$supp_thi_lt_tlo <- DIFF %>%
  apply(2, function(x) sum(x < 0)/nsims) %>%
  (function(x) x >= 0.9 | x <= 0.1)
out$supp_thi_2lt_tlo <- DIFF %>%
  apply(2, function(x) sum(x < -0.02)/nsims) %>%
  (function(x) x >= 0.9)

psi1 <- QSLpersonal::expit(beta_bg + delta_bg * X.trend[1])
psi10 <- QSLpersonal::expit(beta_bg + delta_bg * X.trend[10])
trend_bg <- log(psi10 / psi1)
out$trend_bg <- apply(trend_bg, 2, median)
out$trend_bg_lo <- apply(trend_bg, 2, function(x) quantile(x, prob = 0.1, type = 8))
out$trend_bg_hi <- apply(trend_bg, 2, function(x) quantile(x, prob = 0.9, type = 8))
DIFF <- exp(trend_hi) - exp(trend_bg)
out$supp_thi_lt_tbg <- DIFF %>%
  apply(2, function(x) sum(x < 0)/nsims) %>%
  (function(x) x >= 0.9 | x <= 0.1)
out$supp_thi_2lt_tbg <- DIFF %>%
  apply(2, function(x) sum(x < -0.02)/nsims) %>%
  (function(x) x >= 0.9)

## Filter down to supported species ##
guild.codes <- c("Gn", "Gr", "M", "R", "Sg", "Sh", "Wt", "Wd")
names(guild.codes) <- sort(unique(spp.out$Guild))
guild.codes <- guild.codes[spp.out$Guild]
names(guild.codes) <- spp.out$BirdCode

dat_plot <- out %>%
  filter(supp_Thi_lt_Tlo |
           supp_Thi_lt_Tbg |
           supp_thi_lt_tlo |
           supp_thi_lt_tbg) %>%
  mutate(index = row_number() %>% rev(),
         Spp = str_c(Spp, "(", guild.codes[Spp], ")"))
dat_plot_long <- dat_plot %>%
  select(Spp, index,
         TREND_hi, TREND_hi_lo, TREND_hi_hi,
         trend_hi, trend_hi_lo, trend_hi_hi) %>%
  rename(TREND = TREND_hi, TREND_lo = TREND_hi_lo, TREND_hi = TREND_hi_hi,
         trend = trend_hi, trend_lo = trend_hi_lo, trend_hi = trend_hi_hi) %>%
  mutate(Development = "High") %>%
  bind_rows(dat_plot %>%
              select(Spp, index,
                     TREND_lo, TREND_lo_lo, TREND_lo_hi,
                     trend_lo, trend_lo_lo, trend_lo_hi) %>%
              rename(TREND = TREND_lo, TREND_lo = TREND_lo_lo, TREND_hi = TREND_lo_hi,
                     trend = trend_lo, trend_lo = trend_lo_lo, trend_hi = trend_lo_hi) %>%
              mutate(Development = "Low")) %>%
  bind_rows(dat_plot %>%
              select(Spp, index,
                     TREND_bg, TREND_bg_lo, TREND_bg_hi,
                     trend_bg, trend_bg_lo, trend_bg_hi) %>%
              rename(TREND = TREND_bg, TREND_lo = TREND_bg_lo, TREND_hi = TREND_bg_hi,
                     trend = trend_bg, trend_lo = trend_bg_lo, trend_hi = trend_bg_hi) %>%
              mutate(Development = "Background")) %>%
  mutate(Development = factor(Development, levels = c("Background", "Low", "High")))

# Grid cell trend plot #
dat_supp_lo <- dat_plot %>%
  filter(supp_Thi_lt_Tlo) %>%
  mutate(x = index, y = TREND_lo_hi + 0.05) %>%
  select(x, y)
dat_supp_lo2 <- dat_plot %>%
  filter(supp_Thi_2lt_Tlo) %>%
  mutate(x = index, y = TREND_lo_hi + 0.1) %>%
  select(x, y)
dat_supp_bg <- dat_plot %>%
  filter(supp_Thi_lt_Tbg) %>%
  mutate(x = index - 0.35, y = TREND_bg_hi + 0.05) %>%
  select(x, y)
dat_supp_bg2 <- dat_plot %>%
  filter(supp_Thi_2lt_Tbg) %>%
  mutate(x = index - 0.35, y = TREND_bg_hi + 0.1) %>%
  select(x, y)

dodge <- position_dodge(width=0.7)
ymin <- min(dat_plot_long$TREND_lo)
ymax <- max(dat_plot_long$TREND_hi + 0.1)
breaks <- c(seq(ymin, 0, length.out = 5),
            seq(0, ymax, length.out = 5)[-1])
break.labs <- round((exp(breaks) - 1) * 100)
p.TREND <- ggplot(dat = dat_plot_long, aes(x = index, y = TREND)) +
  geom_errorbar(aes(ymin = TREND_lo, ymax = TREND_hi, color = Development),
                size=1, width=0, position = dodge) +
  geom_point(size = 2.5, aes(color = Development), position = dodge) +
  geom_text(data = dat_supp_lo, aes(x = x, y = y), label = "*", size = 5) +
  geom_text(data = dat_supp_lo2, aes(x = x, y = y), label = "*", size = 5) +
  geom_text(data = dat_supp_bg, aes(x = x, y = y), label = "*", size = 5) +
  geom_text(data = dat_supp_bg2, aes(x = x, y = y), label = "*", size = 5) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = log(0.75), linetype = "dashed") +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat_plot), labels = rev(dat_plot$Spp),
                     expand=c(0, 1)) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = breaks, labels = break.labs) +
  scale_color_manual(values = c("#000000", "#009E73", "#D55E00")) +
  ylab("Coarse-scale trend (%)") + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  theme(legend.position = c(0.01,1), legend.justification = c(0.01,1))

# Point scale trend plot #
dat_supp_lo <- dat_plot %>%
  filter(supp_thi_lt_tlo) %>%
  mutate(x = index, y = trend_lo_hi + 0.1) %>%
  select(x, y)
dat_supp_lo2 <- dat_plot %>%
  filter(supp_thi_2lt_tlo) %>%
  mutate(x = index, y = trend_lo_hi + 0.2) %>%
  select(x, y)
dat_supp_bg <- dat_plot %>%
  filter(supp_thi_lt_tbg) %>%
  mutate(x = index - 0.35, y = trend_bg_hi + 0.1) %>%
  select(x, y)
dat_supp_bg2 <- dat_plot %>%
  filter(supp_thi_2lt_tbg) %>%
  mutate(x = index - 0.35, y = trend_bg_hi + 0.2) %>%
  select(x, y)

dodge <- position_dodge(width=0.7)
ymin <- min(dat_plot_long$trend_lo)
ymax <- max(dat_plot_long$trend_hi + 0.2)
breaks <- c(seq(ymin, 0, length.out = 5),
            seq(0, ymax, length.out = 5)[-1])
break.labs <- round((exp(breaks) - 1) * 100)
p.trend <- ggplot(dat = dat_plot_long, aes(x = index, y = trend)) +
  geom_errorbar(aes(ymin = trend_lo, ymax = trend_hi, color = Development),
                size=1, width=0, position = dodge) +
  geom_point(size = 2.5, aes(color = Development), position = dodge) +
  geom_text(data = dat_supp_lo, aes(x = x, y = y), label = "*", size = 5) +
  geom_text(data = dat_supp_lo2, aes(x = x, y = y), label = "*", size = 5) +
  geom_text(data = dat_supp_bg, aes(x = x, y = y), label = "*", size = 5) +
  geom_text(data = dat_supp_bg2, aes(x = x, y = y), label = "*", size = 5) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = log(0.9), linetype = "dashed") +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat_plot), labels = rev(dat_plot$Spp),
                     expand=c(0, 1)) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = breaks, labels = break.labs) +
  scale_color_manual(values = c("#000000", "#009E73", "#D55E00")) +
  ylab("Fine-scale trend (%)") + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = "none")

p <- ggdraw() +
  draw_plot(p.TREND, x = 0.05, y = 0, width = 0.475, height = 1) +
  draw_plot(p.trend, x = 0.525, y = 0, width = 0.475, height = 1) +
  draw_plot_label("Species", x = 0, y = 0.5, angle = 90, hjust = 0, size = 30)

save_plot("Figure_spp_trends.jpg", p, ncol = 2, nrow = 3, dpi = 600)
