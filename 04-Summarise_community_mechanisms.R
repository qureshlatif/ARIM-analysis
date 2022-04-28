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
guild.mem[which(guild.mem["Sagebrush"]), "Shrubland"] <- TRUE

## Compile evidence for mechanisms ##
cols <- c("Diff_lo_tot", "Diff_bg_tot",
          "Well_1km_supp",
          "Well_1km_lo_effect", "Well_1km_lo_contrib",
          "Well_1km_bg_effect", "Well_1km_bg_contrib",
          "Road_125m_supp",
          "Road_125m_lo_effect", "Road_125m_lo_contrib",
          "Road_125m_bg_effect", "Road_125m_bg_contrib",
          "AHerb_supp",
          "AHerb_lo_effect", "AHerb_lo_contrib",
          "AHerb_bg_effect", "AHerb_bg_contrib")
out <- matrix("", nrow = length(guilds), ncol = length(cols),
              dimnames = list(guilds, cols))

for(g in guilds) {
  spp.ind <- which(guild.mem[,g])
  
  # Differences in trend #
  BETA_hi <- mod$mcmcOutput$BETA0[, spp.ind]
  DELTA_hi <- mod$mcmcOutput$DELTA0[, spp.ind]
  beta_hi <- mod$mcmcOutput$beta0[, spp.ind]
  delta_hi <- mod$mcmcOutput$delta0[, spp.ind]
  for(i in 1:length(spp.ind)) {
    BETA_hi[, i] <- BETA_hi[, i] +
      apply(mod$mcmcOutput$BETA1[,spp.ind[i],] * X.PSI.pred.hi, 1, sum)    
    DELTA_hi[,i] <- DELTA_hi[,i] +
      apply(mod$mcmcOutput$DELTA1[,spp.ind[i],] * X.LAMBDA.pred.hi, 1, sum)
    beta_hi[, i] <- beta_hi[, i] +
      apply(mod$mcmcOutput$beta1[,spp.ind[i],] * X.psi.pred.hi, 1, sum)
    delta_hi[,i] <- delta_hi[,i] +
      apply(mod$mcmcOutput$delta1[,spp.ind[i],] * X.lambda.pred.hi, 1, sum)
  }
  psi1 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[1]) *
    QSLpersonal::expit(beta_hi + delta_hi * X.trend[1]) *
    mod$mcmcOutput$w[, spp.ind]
  psi10 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[10]) *
    QSLpersonal::expit(beta_hi + delta_hi * X.trend[10]) *
    mod$mcmcOutput$w[, spp.ind]
  trend_hi <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
  
  BETA_lo <- mod$mcmcOutput$BETA0[, spp.ind]
  DELTA_lo <- mod$mcmcOutput$DELTA0[, spp.ind]
  beta_lo <- mod$mcmcOutput$beta0[, spp.ind]
  delta_lo <- mod$mcmcOutput$delta0[, spp.ind]
  for(i in 1:length(spp.ind)) {
    BETA_lo[, i] <- BETA_lo[, i] +
      apply(mod$mcmcOutput$BETA1[,spp.ind[i],] * X.PSI.pred.lo, 1, sum)    
    DELTA_lo[,i] <- DELTA_lo[,i] +
      apply(mod$mcmcOutput$DELTA1[,spp.ind[i],] * X.LAMBDA.pred.lo, 1, sum)
    beta_lo[, i] <- beta_lo[, i] +
      apply(mod$mcmcOutput$beta1[,spp.ind[i],] * X.psi.pred.lo, 1, sum)
    delta_lo[,i] <- delta_lo[,i] +
      apply(mod$mcmcOutput$delta1[,spp.ind[i],] * X.lambda.pred.lo, 1, sum)
  }
  psi1 <- QSLpersonal::expit(BETA_lo + DELTA_lo * X.trend[1]) *
    QSLpersonal::expit(beta_lo + delta_lo * X.trend[1]) *
    mod$mcmcOutput$w[, spp.ind]
  psi10 <- QSLpersonal::expit(BETA_lo + DELTA_lo * X.trend[10]) *
    QSLpersonal::expit(beta_lo + delta_lo * X.trend[10]) *
    mod$mcmcOutput$w[, spp.ind]
  trend_lo <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
  diff_lo <- trend_hi - trend_lo
  out[g, "Diff_lo_tot"] <- QSLpersonal::BCI(diff_lo)
  Diff_lo_supp <- str_detect(out[g, "Diff_lo_tot"], "\\)\\*")
  
  BETA_bg <- mod$mcmcOutput$BETA0[, spp.ind]
  DELTA_bg <- mod$mcmcOutput$DELTA0[, spp.ind]
  beta_bg <- mod$mcmcOutput$beta0[, spp.ind]
  delta_bg <- mod$mcmcOutput$delta0[, spp.ind]
  for(i in 1:length(spp.ind)) {
    BETA_bg[, i] <- BETA_bg[, i] +
      apply(mod$mcmcOutput$BETA1[,spp.ind[i],] * X.PSI.pred.bg, 1, sum)    
    DELTA_bg[,i] <- DELTA_bg[,i] +
      apply(mod$mcmcOutput$DELTA1[,spp.ind[i],] * X.LAMBDA.pred.bg, 1, sum)
    beta_bg[, i] <- beta_bg[, i] +
      apply(mod$mcmcOutput$beta1[,spp.ind[i],] * X.psi.pred.bg, 1, sum)
    delta_bg[,i] <- delta_bg[,i] +
      apply(mod$mcmcOutput$delta1[,spp.ind[i],] * X.lambda.pred.bg, 1, sum)
  }
  psi1 <- QSLpersonal::expit(BETA_bg + DELTA_bg * X.trend[1]) *
    QSLpersonal::expit(beta_bg + delta_bg * X.trend[1]) *
    mod$mcmcOutput$w[, spp.ind]
  psi10 <- QSLpersonal::expit(BETA_bg + DELTA_bg * X.trend[10]) *
    QSLpersonal::expit(beta_bg + delta_bg * X.trend[10]) *
    mod$mcmcOutput$w[, spp.ind]
  trend_bg <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
  diff_bg <- trend_hi - trend_bg
  out[g, "Diff_bg_tot"] <- QSLpersonal::BCI(diff_bg)
  Diff_bg_supp <- str_detect(out[g, "Diff_bg_tot"], "\\)\\*")
  
  # Evidence for Well_1km mechanism #
    # Consistency of relationship with mechanism #
  ind.x.delta <- which(dimnames(X.lambda)[[2]] == "Well_1km")
  delta.x <- mod$mcmcOutput$delta1[,spp.ind,ind.x.delta] %>%
    apply(1, mean)
  if(quantile(delta.x, prob = 0.975, type = 8) < 0) {
    out[g, "Well_1km_supp"] <- "supported"
  } else {
    if(quantile(delta.x, prob = 0.025, type = 8) > 0) {
      out[g, "Well_1km_supp"] <- "contradicted"
    } else {
      out[g, "Well_1km_supp"] <- "unclear"
    }
  }
  
  if(out[g, "Well_1km_supp"] == "supported") {
    ind.beta.x <- which(dimnames(X.psi)[[2]] == "Well_1km")
    ind.delta.x <- which(dimnames(X.lambda)[[2]] == "Well_1km")
    if(Diff_lo_supp) {
      beta_hi.rmx <- mod$mcmcOutput$beta0[,spp.ind]
      delta_hi.rmx <- mod$mcmcOutput$delta0[,spp.ind]
      for(i in 1:length(spp.ind)) {
        beta_hi.rmx[,i] <- beta_hi.rmx[,i] +
          apply(mod$mcmcOutput$beta1[,spp.ind[i],-ind.beta.x] * X.psi.pred.hi[,-ind.beta.x], 1, sum) +
          mod$mcmcOutput$beta1[,spp.ind[i],ind.beta.x] * X.psi.pred.lo[,ind.beta.x]
        delta_hi.rmx[,i] <- delta_hi.rmx[,i] +
          apply(mod$mcmcOutput$delta1[,spp.ind[i],-ind.delta.x] * X.lambda.pred.hi[,-ind.delta.x], 1, sum) +
          mod$mcmcOutput$delta1[,spp.ind[i],ind.delta.x] * X.lambda.pred.lo[,ind.delta.x]
      }
      psi1 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[1]) *
        QSLpersonal::expit(beta_hi.rmx + delta_hi.rmx * X.trend[1]) *
        mod$mcmcOutput$w[, spp.ind]
      psi10 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[10]) *
        QSLpersonal::expit(beta_hi.rmx + delta_hi.rmx * X.trend[10]) *
        mod$mcmcOutput$w[, spp.ind]
      trend_hi.rmx <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
      diff_lo.rmx <- trend_hi.rmx - trend_lo
      out[g, "Well_1km_lo_effect"] <- QSLpersonal::BCI(diff_lo - diff_lo.rmx)
      out[g, "Well_1km_lo_contrib"] <- QSLpersonal::BCI(((diff_lo - diff_lo.rmx) / diff_lo) * 100)
    }
    if(Diff_bg_supp) {
      beta_hi.rmx <- mod$mcmcOutput$beta0[,spp.ind]
      delta_hi.rmx <- mod$mcmcOutput$delta0[,spp.ind]
      for(i in 1:length(spp.ind)) {
        beta_hi.rmx[,i] <- beta_hi.rmx[,i] +
          apply(mod$mcmcOutput$beta1[,spp.ind[i],-ind.beta.x] * X.psi.pred.hi[,-ind.beta.x], 1, sum) +
          mod$mcmcOutput$beta1[,spp.ind[i],ind.beta.x] * X.psi.pred.bg[,ind.beta.x]
        delta_hi.rmx[,i] <- delta_hi.rmx[,i] +
          apply(mod$mcmcOutput$delta1[,spp.ind[i],-ind.delta.x] * X.lambda.pred.hi[,-ind.delta.x], 1, sum) +
          mod$mcmcOutput$delta1[,spp.ind[i],ind.delta.x] * X.lambda.pred.bg[,ind.delta.x]
      }
      psi1 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[1]) *
        QSLpersonal::expit(beta_hi.rmx + delta_hi.rmx * X.trend[1]) *
        mod$mcmcOutput$w[, spp.ind]
      psi10 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[10]) *
        QSLpersonal::expit(beta_hi.rmx + delta_hi.rmx * X.trend[10]) *
        mod$mcmcOutput$w[, spp.ind]
      trend_hi.rmx <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
      diff_bg.rmx <- trend_hi.rmx - trend_bg
      out[g, "Well_1km_bg_effect"] <- QSLpersonal::BCI(diff_bg - diff_bg.rmx)
      out[g, "Well_1km_bg_contrib"] <- QSLpersonal::BCI(((diff_bg - diff_bg.rmx) / diff_bg) * 100)
    }
  }

  # Evidence for Road_125m mechanism #
    # Consistency of relationship with mechanism #
  ind.x.delta <- which(dimnames(X.lambda)[[2]] == "Road_125m")
  delta.x <- mod$mcmcOutput$delta1[,spp.ind,ind.x.delta] %>%
    apply(1, mean)
  if(quantile(delta.x, prob = 0.975, type = 8) < 0) {
    out[g, "Road_125m_supp"] <- "supported"
  } else {
    if(quantile(delta.x, prob = 0.025, type = 8) > 0) {
      out[g, "Road_125m_supp"] <- "contradicted"
    } else {
      out[g, "Road_125m_supp"] <- "unclear"
    }
  }
  
  # Evidence for AHerb mechanism #
    # Consistency of relationship with mechanism #
  ind.x.delta <- which(dimnames(X.lambda)[[2]] == "AHerb")
  delta.x <- mod$mcmcOutput$delta1[,spp.ind,ind.x.delta] %>%
    apply(1, mean)
  if(quantile(delta.x, prob = 0.975, type = 8) < 0) {
    out[g, "AHerb_supp"] <- "supported"
  } else {
    if(quantile(delta.x, prob = 0.025, type = 8) > 0) {
      out[g, "AHerb_supp"] <- "contradicted"
    } else {
      out[g, "AHerb_supp"] <- "unclear"
    }
  }
}

write.csv(out, "Community_mechanisms.csv", row.names = T)
