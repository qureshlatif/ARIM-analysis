library(mcmcOutput)
library(stringr)
library(dplyr)
library(R.utils)

setwd("C:/Users/Quresh.Latif/files/projects/BLM/ARIM")
load("Data_compiled.RData")

#__________ Script inputs _____________#
mod.nam <- "mod_path"
scripts.loc <- "ARIM-analysis/"
mod <- loadObject(mod.nam)
nsims <- dim(mod$mcmcOutput)[1]
out.vals <- c("est", "f")

# Data processing
reduce.data.aug <- F # If TRUE, reduce data augmentation to 10 additional species.
development <- F # Set to TRUE for running test model with only develop.spp, and FALSE to run the full model.
#develop.spp <- c("MODO", "MOPL", "HOLA", "WEME")
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
guild.mem[which(guild.mem["Sagebrush"]), "Shrubland"] <- TRUE

## Richness and trend estimates ##
cols <- c("SR_hi", "SR_lo", "SR_bg",
          "TREND_hi", "TREND_lo", "TREND_bg",
          "DIFFTot_lo",  "DIFFTot_bg",
          "sr_hi", "sr_lo", "sr_bg",
          "trend_hi", "trend_lo", "trend_bg",
          "diffTot_lo",  "diffTot_bg")
out <- matrix("", nrow = length(guilds), ncol = length(cols),
              dimnames = list(guilds, cols))

for(g in 1:length(guilds)) {
  spp.ind <- which(guild.mem[,g])
  
  # Grid-cell richness #
  BETA_hi <- mod$mcmcOutput$BETA0[, spp.ind]
  for(i in 1:length(spp.ind)) BETA_hi[, i] <- BETA_hi[, i] +
    apply(mod$mcmcOutput$BETA1[,spp.ind[i],] * X.PSI.pred.hi, 1, sum)
  out[g, "SR_hi"] <-
    QSLpersonal::BCI(
      (QSLpersonal::expit(BETA_hi) *
         mod$mcmcOutput$w[, spp.ind]) %>% apply(1, sum),
      flag.sig = F)
  
  BETA_lo <- mod$mcmcOutput$BETA0[, spp.ind]
  for(i in 1:length(spp.ind)) BETA_lo[, i] <- BETA_lo[, i] +
    apply(mod$mcmcOutput$BETA1[,spp.ind[i],] * X.PSI.pred.lo, 1, sum)
  out[g, "SR_lo"] <-
    QSLpersonal::BCI(
      (QSLpersonal::expit(BETA_lo) *
         mod$mcmcOutput$w[, spp.ind]) %>% apply(1, sum),
      flag.sig = F)
  
  BETA_bg <- mod$mcmcOutput$BETA0[, spp.ind]
  for(i in 1:length(spp.ind)) BETA_bg[, i] <- BETA_bg[, i] +
    apply(mod$mcmcOutput$BETA1[,spp.ind[i],] * X.PSI.pred.bg, 1, sum)
  out[g, "SR_bg"] <-
    QSLpersonal::BCI(
      (QSLpersonal::expit(BETA_bg) *
         mod$mcmcOutput$w[, spp.ind]) %>% apply(1, sum),
      flag.sig = F)

  # Grid-cell trend #
  L.hi <- mod$mcmcOutput$DELTA0[, spp.ind]
  for(i in 1:length(spp.ind)) L.hi[,i] <- L.hi[,i] +
    apply(mod$mcmcOutput$DELTA1[,spp.ind[i],] * X.LAMBDA.pred.hi, 1, sum)
  PSI1 <- QSLpersonal::expit(BETA_hi + L.hi * X.trend[1]) * mod$mcmcOutput$w[, spp.ind]
  PSI10 <- QSLpersonal::expit(BETA_hi + L.hi * X.trend[10]) * mod$mcmcOutput$w[, spp.ind]
  TREND_hi <- apply(PSI10, 1, sum) / apply(PSI1, 1, sum)
  out[g, "TREND_hi"] <- QSLpersonal::BCI(TREND_hi, flag.sig = F)
  if(quantile(TREND_hi, prob = 0.975, type = 8) < 1 |
     quantile(TREND_hi, prob = 0.025, type = 8) > 1)
    out[g, "TREND_hi"] <- str_c(out[g, "TREND_hi"], "*")

  L.lo <- mod$mcmcOutput$DELTA0[, spp.ind]
  for(i in 1:length(spp.ind)) L.lo[,i] <- L.lo[,i] +
    apply(mod$mcmcOutput$DELTA1[,spp.ind[i],] * X.LAMBDA.pred.lo, 1, sum)
  PSI1 <- QSLpersonal::expit(BETA_lo + L.lo * X.trend[1]) * mod$mcmcOutput$w[, spp.ind]
  PSI10 <- QSLpersonal::expit(BETA_lo + L.lo * X.trend[10]) * mod$mcmcOutput$w[, spp.ind]
  TREND_lo <- apply(PSI10, 1, sum) / apply(PSI1, 1, sum)
  out[g, "TREND_lo"] <- QSLpersonal::BCI(TREND_lo, flag.sig = F)
  if(quantile(TREND_lo, prob = 0.975, type = 8) < 1 |
     quantile(TREND_lo, prob = 0.025, type = 8) > 1)
    out[g, "TREND_lo"] <- str_c(out[g, "TREND_lo"], "*")
  
  L.bg <- mod$mcmcOutput$DELTA0[, spp.ind]
  for(i in 1:length(spp.ind)) L.bg[,i] <- L.bg[,i] +
    apply(mod$mcmcOutput$DELTA1[,spp.ind[i],] * X.LAMBDA.pred.bg, 1, sum)
  PSI1 <- QSLpersonal::expit(BETA_bg + L.bg * X.trend[1]) * mod$mcmcOutput$w[, spp.ind]
  PSI10 <- QSLpersonal::expit(BETA_bg + L.bg * X.trend[10]) * mod$mcmcOutput$w[, spp.ind]
  TREND_bg <- apply(PSI10, 1, sum) / apply(PSI1, 1, sum)
  out[g, "TREND_bg"] <- QSLpersonal::BCI(TREND_bg, flag.sig = F)
  if(quantile(TREND_bg, prob = 0.975, type = 8) < 1 |
     quantile(TREND_bg, prob = 0.025, type = 8) > 1)
    out[g, "TREND_bg"] <- str_c(out[g, "TREND_bg"], "*")
  
  # Grid-cell trend differences #
  out[g, "DIFFTot_lo"] <- QSLpersonal::BCI(TREND_hi - TREND_lo)
  out[g, "DIFFTot_bg"] <- QSLpersonal::BCI(TREND_hi - TREND_bg)

  # Point-scale richness #
  beta_hi <- mod$mcmcOutput$beta0[, spp.ind]
  for(i in 1:length(spp.ind)) beta_hi[, i] <- beta_hi[, i] +
    apply(mod$mcmcOutput$beta1[,spp.ind[i],] * X.psi.pred.hi, 1, sum)
  beta_hi <- beta_hi
  out[g, "sr_hi"] <-
    QSLpersonal::BCI(
      (QSLpersonal::expit(BETA_hi) *
        QSLpersonal::expit(beta_hi) *
         mod$mcmcOutput$w[, spp.ind]) %>% apply(1, sum),
      flag.sig = F)
  
  beta_lo <- mod$mcmcOutput$beta0[, spp.ind]
  for(i in 1:length(spp.ind)) beta_lo[, i] <- beta_lo[, i] +
    apply(mod$mcmcOutput$beta1[,spp.ind[i],] * X.psi.pred.lo, 1, sum)
  beta_lo <- beta_lo
  out[g, "sr_lo"] <-
    QSLpersonal::BCI(
      (QSLpersonal::expit(BETA_lo) *
         QSLpersonal::expit(beta_lo) *
         mod$mcmcOutput$w[, spp.ind]) %>% apply(1, sum),
      flag.sig = F)
  
  beta_bg <- mod$mcmcOutput$beta0[, spp.ind]
  for(i in 1:length(spp.ind)) beta_bg[, i] <- beta_bg[, i] +
    apply(mod$mcmcOutput$beta1[,spp.ind[i],] * X.psi.pred.bg, 1, sum)
  beta_bg <- beta_bg
  out[g, "sr_bg"] <-
    QSLpersonal::BCI(
      (QSLpersonal::expit(BETA_bg) *
         QSLpersonal::expit(beta_bg) *
         mod$mcmcOutput$w[, spp.ind]) %>% apply(1, sum),
      flag.sig = F)
  
  # Point-scale trend #
  l.hi <- mod$mcmcOutput$delta0[, spp.ind]
  for(i in 1:length(spp.ind)) l.hi[,i] <- l.hi[,i] +
    apply(mod$mcmcOutput$delta1[,spp.ind[i],] * X.lambda.pred.hi, 1, sum)
  psi1 <- QSLpersonal::expit(BETA_hi + L.hi * X.trend[1]) *
    QSLpersonal::expit(beta_hi + l.hi * X.trend[1]) *
    mod$mcmcOutput$w[, spp.ind]
  psi10 <- QSLpersonal::expit(BETA_hi + L.hi * X.trend[10]) *
    QSLpersonal::expit(beta_hi + l.hi * X.trend[10]) *
    mod$mcmcOutput$w[, spp.ind]
  trend_hi <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
  out[g, "trend_hi"] <- QSLpersonal::BCI(trend_hi, flag.sig = F)
  if(quantile(trend_hi, prob = 0.975, type = 8) < 1 |
     quantile(trend_hi, prob = 0.025, type = 8) > 1)
    out[g, "trend_hi"] <- str_c(out[g, "trend_hi"], "*")
  
  l.lo <- mod$mcmcOutput$delta0[, spp.ind]
  for(i in 1:length(spp.ind)) l.lo[,i] <- l.lo[,i] +
    apply(mod$mcmcOutput$delta1[,spp.ind[i],] * X.lambda.pred.lo, 1, sum)
  psi1 <- QSLpersonal::expit(BETA_lo + L.lo * X.trend[1]) *
    QSLpersonal::expit(beta_lo + l.lo * X.trend[1]) *
    mod$mcmcOutput$w[, spp.ind]
  psi10 <- QSLpersonal::expit(BETA_lo + L.lo * X.trend[10]) *
    QSLpersonal::expit(beta_lo + l.lo * X.trend[10]) *
    mod$mcmcOutput$w[, spp.ind]
  trend_lo <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
  out[g, "trend_lo"] <- QSLpersonal::BCI(trend_lo, flag.sig = F)
  if(quantile(trend_lo, prob = 0.975, type = 8) < 1 |
     quantile(trend_lo, prob = 0.025, type = 8) > 1)
    out[g, "trend_lo"] <- str_c(out[g, "trend_lo"], "*")
  
  l.bg <- mod$mcmcOutput$delta0[, spp.ind]
  for(i in 1:length(spp.ind)) l.bg[,i] <- l.bg[,i] +
    apply(mod$mcmcOutput$delta1[,spp.ind[i],] * X.lambda.pred.bg, 1, sum)
  psi1 <- QSLpersonal::expit(BETA_bg + L.bg * X.trend[1]) *
    QSLpersonal::expit(beta_bg + l.bg * X.trend[1]) *
    mod$mcmcOutput$w[, spp.ind]
  psi10 <- QSLpersonal::expit(BETA_bg + L.bg * X.trend[10]) *
    QSLpersonal::expit(beta_bg + l.bg * X.trend[10]) *
    mod$mcmcOutput$w[, spp.ind]
  trend_bg <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
  out[g, "trend_bg"] <- QSLpersonal::BCI(trend_bg, flag.sig = F)
  if(quantile(trend_bg, prob = 0.975, type = 8) < 1 |
     quantile(trend_bg, prob = 0.025, type = 8) > 1)
    out[g, "trend_bg"] <- str_c(out[g, "trend_bg"], "*")
  
  # Point-scale trend differences #
  out[g, "diffTot_lo"] <- QSLpersonal::BCI(trend_hi - trend_lo)
  out[g, "diffTot_bg"] <- QSLpersonal::BCI(trend_hi - trend_bg)
}

write.csv(out, "Summary_community_dev_effects.csv", row.names = T)

