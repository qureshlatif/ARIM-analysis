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

# Tabulate parameter estimates
pars <- c(str_c("BETA.", dimnames(X.PSI)[[2]]),
          str_c("DELTA.", dimnames(X.LAMBDA)[[2]]),
          str_c("beta.", dimnames(X.psi)[[2]]),
          str_c("delta.", dimnames(X.lambda)[[2]]))
cols <- (c("", ".lo", ".hi") %>%
           expand.grid(pars, stringsAsFactors = F) %>%
           mutate(Var3 = str_c(Var2, Var1, sep = "")))$Var3
tbl_pars <- matrix(NA, nrow = length(spp.list), ncol = length(cols), dimnames = list(NULL, cols))

pars.ind <- which(str_detect(pars, "BETA"))
for(i in 1:length(pars.ind)) {
  parm <- mod$mcmcOutput$BETA1[,,i]
  tbl_pars[, pars[pars.ind[i]]] <- apply(parm, 2, median)
  tbl_pars[, str_c(pars[pars.ind[i]], ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.1, type = 8))
  tbl_pars[, str_c(pars[pars.ind[i]], ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.9, type = 8))
}

pars.ind <- which(str_detect(pars, "DELTA"))
for(i in 1:length(pars.ind)) {
  parm <- mod$mcmcOutput$DELTA1[,,i]
  tbl_pars[, pars[pars.ind[i]]] <- apply(parm, 2, median)
  tbl_pars[, str_c(pars[pars.ind[i]], ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.1, type = 8))
  tbl_pars[, str_c(pars[pars.ind[i]], ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.9, type = 8))
}

pars.ind <- which(str_detect(pars, "beta"))
for(i in 1:length(pars.ind)) {
  parm <- mod$mcmcOutput$beta1[,,i]
  tbl_pars[, pars[pars.ind[i]]] <- apply(parm, 2, median)
  tbl_pars[, str_c(pars[pars.ind[i]], ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.1, type = 8))
  tbl_pars[, str_c(pars[pars.ind[i]], ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.9, type = 8))
}

pars.ind <- which(str_detect(pars, "delta"))
for(i in 1:length(pars.ind)) {
  parm <- mod$mcmcOutput$delta1[,,i]
  tbl_pars[, pars[pars.ind[i]]] <- apply(parm, 2, median)
  tbl_pars[, str_c(pars[pars.ind[i]], ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.1, type = 8))
  tbl_pars[, str_c(pars[pars.ind[i]], ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.9, type = 8))
}

rm(parm)
spp.detected <- apply(Y.mat, 2, function(x) any(x > 0))
tbl_pars <- tbl_pars[which(spp.detected),]

#### Stratum effects on trend ####
pars.sub <- c("DELTA.Dev_lo", "DELTA.Dev_bg", "delta.Dev_lo", "delta.Dev_bg")
dat.plt <- tbl_pars %>% as_tibble() %>%
  select(contains(pars.sub, ignore.case = F))%>%
  mutate(Spp = spp.list[which(spp.detected)]) %>%
  mutate(index = row_number() %>% rev())

dat.plt.supp <- dat.plt %>%
  filter_at(vars(ends_with(".lo")), any_vars(. > 0)) %>%
  bind_rows(
    dat.plt %>%
      filter_at(vars(ends_with(".hi")), any_vars(. < 0))
  ) %>%
  distinct() %>%
  arrange(index %>% desc()) %>%
  mutate(index = row_number() %>% rev())

cols <- pars.sub %>% str_c(".supp")
dat.supp <- matrix("none", nrow = nrow(dat.plt.supp), ncol = length(cols),
                   dimnames = list(NULL, cols))
for(i in 1:length(cols)) {
  col.chck <- str_sub(cols[i], 1, -6)
  chck <- dat.plt.supp[, which(str_detect(names(dat.plt.supp), col.chck))]
  dat.supp[which(chck[, 2] > 0), cols[i]] <- "pos"
  dat.supp[which(chck[, 3] < 0), cols[i]] <- "neg"
}
dat.plt.supp <- dat.plt.supp %>%
  bind_cols(
    dat.supp %>% data.frame(stringsAsFactors = F)
  )

p.DDlo <- ggplot(dat = dat.plt.supp, aes(x = index, y = DELTA.Dev_lo)) +
  geom_errorbar(aes(ymin = DELTA.Dev_lo.lo, ymax = DELTA.Dev_lo.hi, color = DELTA.Dev_lo.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = DELTA.Dev_lo.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#000000", "#D55E00")) +
  ylab(expression(hat(Delta)["Low-development"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.DDbg <- ggplot(dat = dat.plt.supp, aes(x = index, y = DELTA.Dev_bg)) +
  geom_errorbar(aes(ymin = DELTA.Dev_bg.lo, ymax = DELTA.Dev_bg.hi, color = DELTA.Dev_bg.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = DELTA.Dev_bg.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#000000", "#D55E00")) +
  ylab(expression(hat(Delta)["Background"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.dDlo <- ggplot(dat = dat.plt.supp, aes(x = index, y = delta.Dev_lo)) +
  geom_errorbar(aes(ymin = delta.Dev_lo.lo, ymax = delta.Dev_lo.hi, color = delta.Dev_lo.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = delta.Dev_lo.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#000000", "#D55E00")) +
  ylab(expression(hat(delta)["Low-development"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.dDbg <- ggplot(dat = dat.plt.supp, aes(x = index, y = delta.Dev_bg)) +
  geom_errorbar(aes(ymin = delta.Dev_bg.lo, ymax = delta.Dev_bg.hi, color = delta.Dev_bg.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = delta.Dev_bg.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(delta)["Background"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=30)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p <- ggdraw() + 
  draw_plot(p.DDlo, x = 0.0500, y = 0, width = 0.2375, height = 1) +
  draw_plot(p.DDbg, x = 0.2875, y = 0, width = 0.2375, height = 1) +
  draw_plot(p.dDlo, x = 0.5250, y = 0, width = 0.2375, height = 1) +
  draw_plot(p.dDbg, x = 0.7625, y = 0, width = 0.2375, height = 1) +
  draw_plot_label("Species", x = 0, y = 0.5, size = 40, angle = 90, hjust = 0)

save_plot("Plot_TrendStrataEffects.jpg", p, ncol = 3, nrow = 3.5, dpi = 200)

#### Mechanistic covariate effects on occupancy ####
pars.sub <- c("BETA.Well_3km", "BETA.Road_1km",
              "beta.Well_1km", "beta.Road_125m", "beta.AHerb")
dat.plt <- tbl_pars %>% as_tibble() %>%
  select(contains(pars.sub, ignore.case = F))%>%
  mutate(Spp = spp.list[which(spp.detected)]) %>%
  mutate(index = row_number() %>% rev())

dat.plt.supp <- dat.plt %>%
  filter_at(vars(ends_with(".lo")), any_vars(. > 0)) %>%
  bind_rows(
    dat.plt %>%
      filter_at(vars(ends_with(".hi")), any_vars(. < 0))
  ) %>%
  distinct() %>%
  arrange(index %>% desc()) %>%
  mutate(index = row_number() %>% rev())

cols <- pars.sub %>% str_c(".supp")
dat.supp <- matrix("none", nrow = nrow(dat.plt.supp), ncol = length(cols),
                   dimnames = list(NULL, cols))
for(i in 1:length(cols)) {
  col.chck <- str_sub(cols[i], 1, -6)
  chck <- dat.plt.supp[, which(str_detect(names(dat.plt.supp), col.chck))]
  dat.supp[which(chck[, 2] > 0), cols[i]] <- "pos"
  dat.supp[which(chck[, 3] < 0), cols[i]] <- "neg"
}
dat.plt.supp <- dat.plt.supp %>%
  bind_cols(
    dat.supp %>% data.frame(stringsAsFactors = F)
  )

p.BWell <- ggplot(dat = dat.plt.supp, aes(x = index, y = BETA.Well_3km)) +
  geom_errorbar(aes(ymin = BETA.Well_3km.lo, ymax = BETA.Well_3km.hi, color = BETA.Well_3km.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = BETA.Well_3km.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(Beta)["Well density (900 ha)"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.BRoad <- ggplot(dat = dat.plt.supp, aes(x = index, y = BETA.Road_1km)) +
  geom_errorbar(aes(ymin = BETA.Road_1km.lo, ymax = BETA.Road_1km.hi, color = BETA.Road_1km.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = BETA.Road_1km.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(Beta)["Road density (100 ha)"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.bWell <- ggplot(dat = dat.plt.supp, aes(x = index, y = beta.Well_1km)) +
  geom_errorbar(aes(ymin = beta.Well_1km.lo, ymax = beta.Well_1km.hi, color = beta.Well_1km.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = beta.Well_1km.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(beta)["Well density (100 ha)"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.bRoad <- ggplot(dat = dat.plt.supp, aes(x = index, y = beta.Road_125m)) +
  geom_errorbar(aes(ymin = beta.Road_125m.lo, ymax = beta.Road_125m.hi, color = beta.Road_125m.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = beta.Road_125m.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(beta)["Road density (5 ha)"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.bAHerb <- ggplot(dat = dat.plt.supp, aes(x = index, y = beta.AHerb)) +
  geom_errorbar(aes(ymin = beta.AHerb.lo, ymax = beta.AHerb.hi, color = beta.AHerb.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = beta.AHerb.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(beta)["Cheatgrass"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p <- ggdraw() + 
  draw_plot(p.BWell,  x = 0.05, y = 0, width = 0.19, height = 1) +
  draw_plot(p.BRoad,  x = 0.24, y = 0, width = 0.19, height = 1) +
  draw_plot(p.bWell,  x = 0.43, y = 0, width = 0.19, height = 1) +
  draw_plot(p.bRoad,  x = 0.62, y = 0, width = 0.19, height = 1) +
  draw_plot(p.bAHerb, x = 0.81, y = 0, width = 0.19, height = 1) +
  draw_plot_label("Species", x = 0, y = 0.5, size = 40, angle = 90, hjust = 0)

save_plot("Plot_MechEffects_occupancy.jpg", p, ncol = 3, nrow = 3.5, dpi = 200)

#### Mechanistic covariate effects on trend ####
pars.sub <- c("DELTA.Well_3km", "delta.Well_1km", "delta.Road_125m", "delta.AHerb")
dat.plt <- tbl_pars %>% as_tibble() %>%
  select(contains(pars.sub, ignore.case = F))%>%
  mutate(Spp = spp.list[which(spp.detected)]) %>%
  mutate(index = row_number() %>% rev())

dat.plt.supp <- dat.plt %>%
  filter_at(vars(ends_with(".lo")), any_vars(. > 0)) %>%
  bind_rows(
    dat.plt %>%
      filter_at(vars(ends_with(".hi")), any_vars(. < 0))
  ) %>%
  distinct() %>%
  arrange(index %>% desc()) %>%
  mutate(index = row_number() %>% rev())

cols <- pars.sub %>% str_c(".supp")
dat.supp <- matrix("none", nrow = nrow(dat.plt.supp), ncol = length(cols),
                   dimnames = list(NULL, cols))
for(i in 1:length(cols)) {
  col.chck <- str_sub(cols[i], 1, -6)
  chck <- dat.plt.supp[, which(str_detect(names(dat.plt.supp), col.chck))]
  dat.supp[which(chck[, 2] > 0), cols[i]] <- "pos"
  dat.supp[which(chck[, 3] < 0), cols[i]] <- "neg"
}
dat.plt.supp <- dat.plt.supp %>%
  bind_cols(
    dat.supp %>% data.frame(stringsAsFactors = F)
  )

p.DWell <- ggplot(dat = dat.plt.supp, aes(x = index, y = DELTA.Well_3km)) +
  geom_errorbar(aes(ymin = DELTA.Well_3km.lo, ymax = DELTA.Well_3km.hi, color = DELTA.Well_3km.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = DELTA.Well_3km.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#000000", "#D55E00")) +
  ylab(expression(hat(Delta)["Well density (900 ha)"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.dWell <- ggplot(dat = dat.plt.supp, aes(x = index, y = delta.Well_1km)) +
  geom_errorbar(aes(ymin = delta.Well_1km.lo, ymax = delta.Well_1km.hi, color = delta.Well_1km.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = delta.Well_1km.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(delta)["Well density (100 ha)"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.dRoad <- ggplot(dat = dat.plt.supp, aes(x = index, y = delta.Road_125m)) +
  geom_errorbar(aes(ymin = delta.Road_125m.lo, ymax = delta.Road_125m.hi, color = delta.Road_125m.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = delta.Road_125m.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(delta)["Road density (5 ha)"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.dAHerb <- ggplot(dat = dat.plt.supp, aes(x = index, y = delta.AHerb)) +
  geom_errorbar(aes(ymin = delta.AHerb.lo, ymax = delta.AHerb.hi, color = delta.AHerb.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = delta.AHerb.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(delta)["Cheatgrass"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p <- ggdraw() + 
  draw_plot(p.DWell,  x = 0.0500, y = 0, width = 0.2375, height = 1) +
  draw_plot(p.dWell,  x = 0.2875, y = 0, width = 0.2375, height = 1) +
  draw_plot(p.dRoad,  x = 0.5250, y = 0, width = 0.2375, height = 1) +
  draw_plot(p.dAHerb, x = 0.7625, y = 0, width = 0.2375, height = 1) +
  draw_plot_label("Species", x = 0, y = 0.5, size = 40, angle = 90, hjust = 0)

save_plot("Plot_MechEffects_trend.jpg", p, ncol = 2.5, nrow = 3.5, dpi = 200)

#### Control covariate effects on occupancy ####
pars.sub <- c("BETA.PJ_area", "BETA.NDVI", "beta.TPI_min",
              "beta.Sage", "beta.Herb")
dat.plt <- tbl_pars %>% as_tibble() %>%
  select(contains(pars.sub, ignore.case = F))%>%
  mutate(Spp = spp.list[which(spp.detected)]) %>%
  mutate(index = row_number() %>% rev())

dat.plt.supp <- dat.plt %>%
  filter_at(vars(ends_with(".lo")), any_vars(. > 0)) %>%
  bind_rows(
    dat.plt %>%
      filter_at(vars(ends_with(".hi")), any_vars(. < 0))
  ) %>%
  distinct() %>%
  arrange(index %>% desc()) %>%
  mutate(index = row_number() %>% rev())

cols <- pars.sub %>% str_c(".supp")
dat.supp <- matrix("none", nrow = nrow(dat.plt.supp), ncol = length(cols),
                   dimnames = list(NULL, cols))
for(i in 1:length(cols)) {
  col.chck <- str_sub(cols[i], 1, -6)
  chck <- dat.plt.supp[, which(str_detect(names(dat.plt.supp), col.chck))]
  dat.supp[which(chck[, 2] > 0), cols[i]] <- "pos"
  dat.supp[which(chck[, 3] < 0), cols[i]] <- "neg"
}
dat.plt.supp <- dat.plt.supp %>%
  bind_cols(
    dat.supp %>% data.frame(stringsAsFactors = F)
  )

p.BPJA <- ggplot(dat = dat.plt.supp, aes(x = index, y = BETA.PJ_area)) +
  geom_errorbar(aes(ymin = BETA.PJ_area.lo, ymax = BETA.PJ_area.hi, color = BETA.PJ_area.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = BETA.PJ_area.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(Beta)["Pinyon-juniper forest"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.BNDVI <- ggplot(dat = dat.plt.supp, aes(x = index, y = BETA.NDVI)) +
  geom_errorbar(aes(ymin = BETA.NDVI.lo, ymax = BETA.NDVI.hi, color = BETA.NDVI.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = BETA.NDVI.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(Beta)["NDVI"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.bTPI <- ggplot(dat = dat.plt.supp, aes(x = index, y = beta.TPI_min)) +
  geom_errorbar(aes(ymin = beta.TPI_min.lo, ymax = beta.TPI_min.hi, color = beta.TPI_min.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = beta.TPI_min.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(beta)["Topographic position index"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.bSage <- ggplot(dat = dat.plt.supp, aes(x = index, y = beta.Sage)) +
  geom_errorbar(aes(ymin = beta.Sage.lo, ymax = beta.Sage.hi, color = beta.Sage.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = beta.Sage.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(beta)["Sagebrush"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p.bHerb <- ggplot(dat = dat.plt.supp, aes(x = index, y = beta.Herb)) +
  geom_errorbar(aes(ymin = beta.Herb.lo, ymax = beta.Herb.hi, color = beta.Herb.supp), size=1, width=0) +
  geom_point(size = 2.5, aes(color = beta.Herb.supp)) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.supp), labels = rev(dat.plt.supp$Spp),
                     expand=c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(beta)["Herbaceous cover"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  theme(axis.text.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=15)) +
  guides(color = F)

p <- ggdraw() + 
  draw_plot(p.BPJA,  x = 0.05, y = 0, width = 0.19, height = 1) +
  draw_plot(p.BNDVI, x = 0.24, y = 0, width = 0.19, height = 1) +
  draw_plot(p.bTPI,  x = 0.43, y = 0, width = 0.19, height = 1) +
  draw_plot(p.bSage, x = 0.62, y = 0, width = 0.19, height = 1) +
  draw_plot(p.bHerb, x = 0.81, y = 0, width = 0.19, height = 1) +
  draw_plot_label("Species", x = 0, y = 0.5, size = 40, angle = 90, hjust = 0)

save_plot("Plot_ControlEffects.jpg", p, ncol = 4, nrow = 4.5, dpi = 200)
