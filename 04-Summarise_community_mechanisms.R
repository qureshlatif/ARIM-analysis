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
guild.mem[which(guild.mem[,"Sagebrush"]), "Shrubland"] <- TRUE

## Compile evidence for mechanisms ##
cols <- c("DIFF_lo", "DIFF_plo",
          "Well_3km_contrib_Dlo", "Well_3km_pcont_Dlo",
          "DIFF_bg", "DIFF_pbg",
          "Well_3km_contrib_Dbg", "Well_3km_pcont_Dbg",
          "Diff_lo", "Diff_plo",
          "Well_1km_contrib_dlo", "Well_1km_pcont_dlo",
          "Diff_bg", "Diff_pbg",
          "Well_1km_contrib_dbg", "Well_1km_pcont_dbg",
          "Road_125m_contrib_dbg", "Road_125m_pcont_dbg",
          "AHerb_contrib_dbg", "AHerb_pcont_dbg")
out <- matrix("", nrow = length(guilds), ncol = length(cols),
              dimnames = list(guilds, cols))

for(g in guilds) {
  spp.ind <- which(guild.mem[,g])
  
  ## Coefficients ##
  BETA_hi <- mod$mcmcOutput$BETA0[, spp.ind]
  DELTA_hi <- mod$mcmcOutput$DELTA0[, spp.ind]
  beta_hi <- mod$mcmcOutput$beta0[, spp.ind]
  delta_hi <- mod$mcmcOutput$delta0[, spp.ind]

  BETA_lo <- mod$mcmcOutput$BETA0[, spp.ind]
  DELTA_lo <- mod$mcmcOutput$DELTA0[, spp.ind]
  beta_lo <- mod$mcmcOutput$beta0[, spp.ind]
  delta_lo <- mod$mcmcOutput$delta0[, spp.ind]
  
  BETA_bg <- mod$mcmcOutput$BETA0[, spp.ind]
  DELTA_bg <- mod$mcmcOutput$DELTA0[, spp.ind]
  beta_bg <- mod$mcmcOutput$beta0[, spp.ind]
  delta_bg <- mod$mcmcOutput$delta0[, spp.ind]
  
  for(i in 1:length(spp.ind)) {
    BETA_hi[, i] <- BETA_hi[, i] +
      apply(mod$mcmcOutput$BETA1[,spp.ind[i],] * X.PSI.pred.hi, 1, sum)    
    DELTA_hi[,i] <- DELTA_hi[,i] +
      apply(mod$mcmcOutput$DELTA1[,spp.ind[i],] * X.LAMBDA.pred.hi, 1, sum)
    beta_hi[, i] <- beta_hi[, i] +
      apply(mod$mcmcOutput$beta1[,spp.ind[i],] * X.psi.pred.hi, 1, sum)
    delta_hi[,i] <- delta_hi[,i] +
      apply(mod$mcmcOutput$delta1[,spp.ind[i],] * X.lambda.pred.hi, 1, sum)
    
    BETA_lo[, i] <- BETA_lo[, i] +
      apply(mod$mcmcOutput$BETA1[,spp.ind[i],] * X.PSI.pred.lo, 1, sum)    
    DELTA_lo[,i] <- DELTA_lo[,i] +
      apply(mod$mcmcOutput$DELTA1[,spp.ind[i],] * X.LAMBDA.pred.lo, 1, sum)
    beta_lo[, i] <- beta_lo[, i] +
      apply(mod$mcmcOutput$beta1[,spp.ind[i],] * X.psi.pred.lo, 1, sum)
    delta_lo[,i] <- delta_lo[,i] +
      apply(mod$mcmcOutput$delta1[,spp.ind[i],] * X.lambda.pred.lo, 1, sum)
    
    BETA_bg[, i] <- BETA_bg[, i] +
      apply(mod$mcmcOutput$BETA1[,spp.ind[i],] * X.PSI.pred.bg, 1, sum)    
    DELTA_bg[,i] <- DELTA_bg[,i] +
      apply(mod$mcmcOutput$DELTA1[,spp.ind[i],] * X.LAMBDA.pred.bg, 1, sum)
    beta_bg[, i] <- beta_bg[, i] +
      apply(mod$mcmcOutput$beta1[,spp.ind[i],] * X.psi.pred.bg, 1, sum)
    delta_bg[,i] <- delta_bg[,i] +
      apply(mod$mcmcOutput$delta1[,spp.ind[i],] * X.lambda.pred.bg, 1, sum)
  }
  
  ## Grid scale ##
  # Differences in trend #
  PSI1 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[1]) *
    mod$mcmcOutput$w[, spp.ind]
  PSI10 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[10]) *
    mod$mcmcOutput$w[, spp.ind]
  TREND_hi <- apply(PSI10, 1, sum) / apply(PSI1, 1, sum)
  
  PSI1 <- QSLpersonal::expit(BETA_lo + DELTA_lo * X.trend[1]) *
    mod$mcmcOutput$w[, spp.ind]
  PSI10 <- QSLpersonal::expit(BETA_lo + DELTA_lo * X.trend[10]) *
    mod$mcmcOutput$w[, spp.ind]
  TREND_lo <- apply(PSI10, 1, sum) / apply(PSI1, 1, sum)
  DIFF_lo <- TREND_hi - TREND_lo
  out[g, "DIFF_lo"] <- QSLpersonal::BCI(DIFF_lo, BCIpercent = 80)
  out[g, "DIFF_plo"] <- round(sum(DIFF_lo < 0) / nsims, digits = 2)

  PSI1 <- QSLpersonal::expit(BETA_bg + DELTA_bg * X.trend[1]) *
    mod$mcmcOutput$w[, spp.ind]
  PSI10 <- QSLpersonal::expit(BETA_bg + DELTA_bg * X.trend[10]) *
    mod$mcmcOutput$w[, spp.ind]
  TREND_bg <- apply(PSI10, 1, sum) / apply(PSI1, 1, sum)
  DIFF_bg <- TREND_hi - TREND_bg
  out[g, "DIFF_bg"] <- QSLpersonal::BCI(DIFF_bg, BCIpercent = 80)
  out[g, "DIFF_pbg"] <- round(sum(DIFF_bg < 0) / nsims, digits = 2)

  # Evidence for mechanisms #
  if(out[g, "DIFF_plo"] >= 0.9) {
    # Well density contribution to difference from low development
    ind.bx <- which(dimnames(X.PSI)[[2]] == "Well_3km")
    ind.dx <- which(dimnames(X.LAMBDA)[[2]] == "Well_3km")
    bhi.rmx <- mod$mcmcOutput$BETA0[,spp.ind]
    dhi.rmx <- mod$mcmcOutput$DELTA0[,spp.ind]
    for(i in 1:length(spp.ind)) {
      bhi.rmx[,i] <- bhi.rmx[,i] + apply(mod$mcmcOutput$BETA1[,spp.ind[i],-ind.bx] *
                                   X.PSI.pred.hi[,-ind.bx], 1, sum) +
        mod$mcmcOutput$BETA1[,spp.ind[i],ind.bx] * X.PSI.pred.lo[,ind.bx]
      dhi.rmx[,i] <- dhi.rmx[,i] + apply(mod$mcmcOutput$DELTA1[,spp.ind[i],-ind.dx] *
                                   X.LAMBDA.pred.hi[,-ind.dx], 1, sum) +
        mod$mcmcOutput$DELTA1[,spp.ind[i],ind.dx] * X.LAMBDA.pred.lo[,ind.dx]
    }
    psi1 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[1]) *
      mod$mcmcOutput$w[, spp.ind]
    psi10 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[10]) *
      mod$mcmcOutput$w[, spp.ind]
    trend_hi.rmx <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
    diff_lo.rmx <- trend_hi.rmx - TREND_lo
    out[g, "Well_3km_contrib_Dlo"] <- QSLpersonal::BCI(((DIFF_lo - diff_lo.rmx) / DIFF_lo) * 100, BCIpercent = 80)
    out[g, "Well_3km_pcont_Dlo"] <- round(sum(((DIFF_lo - diff_lo.rmx) / DIFF_lo) > 0) / nsims, digits = 2)
  } else {
    out[g, "Well_3km_contrib_Dlo"] <- "not considered"
    out[g, "Well_3km_pcont_Dlo"] <-  "not considered"
  }
  if(out[g, "DIFF_pbg"] >= 0.9) {
    # Well density contribution to difference from background
    ind.bx <- which(dimnames(X.PSI)[[2]] == "Well_3km")
    ind.dx <- which(dimnames(X.LAMBDA)[[2]] == "Well_3km")
    bhi.rmx <- mod$mcmcOutput$BETA0[,spp.ind]
    dhi.rmx <- mod$mcmcOutput$DELTA0[,spp.ind]
    for(i in 1:length(spp.ind)) {
      bhi.rmx[,i] <- bhi.rmx[,i] + apply(mod$mcmcOutput$BETA1[,spp.ind[i],-ind.bx] *
                                   X.PSI.pred.hi[,-ind.bx], 1, sum) +
        mod$mcmcOutput$BETA1[,spp.ind[i],ind.bx] * X.PSI.pred.bg[,ind.bx]
      dhi.rmx[,i] <- dhi.rmx[,i] + apply(mod$mcmcOutput$DELTA1[,spp.ind[i],-ind.dx] *
                                   X.LAMBDA.pred.hi[,-ind.dx], 1, sum) +
        mod$mcmcOutput$DELTA1[,spp.ind[i],ind.dx] * X.LAMBDA.pred.bg[,ind.dx]
    }
    psi1 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[1]) *
      mod$mcmcOutput$w[, spp.ind]
    psi10 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[10]) *
      mod$mcmcOutput$w[, spp.ind]
    trend_hi.rmx <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
    diff_bg.rmx <- trend_hi.rmx - TREND_bg
    out[g, "Well_3km_contrib_Dbg"] <- QSLpersonal::BCI(((DIFF_bg - diff_bg.rmx) / DIFF_bg) * 100, BCIpercent = 80)
    out[g, "Well_3km_pcont_Dbg"] <- round(sum(((DIFF_bg - diff_bg.rmx) / DIFF_bg) > 0) / nsims, digits = 2)
  } else {
    out[g, "Well_3km_contrib_Dbg"] <- "not considered"
    out[g, "Well_3km_pcont_Dbg"] <-  "not considered"
  }
  
  ## Point scale ##
  # Differences in trend #
  psi1 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[1]) *
    QSLpersonal::expit(beta_hi + delta_hi * X.trend[1]) *
    mod$mcmcOutput$w[, spp.ind]
  psi10 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[10]) *
    QSLpersonal::expit(beta_hi + delta_hi * X.trend[10]) *
    mod$mcmcOutput$w[, spp.ind]
  trend_hi <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
  
  psi1 <- QSLpersonal::expit(BETA_lo + DELTA_lo * X.trend[1]) *
    QSLpersonal::expit(beta_lo + delta_lo * X.trend[1]) *
    mod$mcmcOutput$w[, spp.ind]
  psi10 <- QSLpersonal::expit(BETA_lo + DELTA_lo * X.trend[10]) *
    QSLpersonal::expit(beta_lo + delta_lo * X.trend[10]) *
    mod$mcmcOutput$w[, spp.ind]
  trend_lo <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
  diff_lo <- trend_hi - trend_lo
  out[g, "Diff_lo"] <- QSLpersonal::BCI(diff_lo, BCIpercent = 80)
  out[g, "Diff_plo"] <- round(sum(diff_lo < 0) / nsims, digits = 2)

  psi1 <- QSLpersonal::expit(BETA_bg + DELTA_bg * X.trend[1]) *
    QSLpersonal::expit(beta_bg + delta_bg * X.trend[1]) *
    mod$mcmcOutput$w[, spp.ind]
  psi10 <- QSLpersonal::expit(BETA_bg + DELTA_bg * X.trend[10]) *
    QSLpersonal::expit(beta_bg + delta_bg * X.trend[10]) *
    mod$mcmcOutput$w[, spp.ind]
  trend_bg <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
  diff_bg <- trend_hi - trend_bg
  out[g, "Diff_bg"] <- QSLpersonal::BCI(diff_bg, BCIpercent = 80)
  out[g, "Diff_pbg"] <- round(sum(diff_bg < 0) / nsims, digits = 2)

  # Evidence for mechanisms #
  if(out[g, "Diff_plo"] >= 0.9) {
    # Well density contribution to difference from low development
    ind.bx <- which(dimnames(X.psi)[[2]] == "Well_1km")
    ind.dx <- which(dimnames(X.lambda)[[2]] == "Well_1km")
    bhi.rmx <- mod$mcmcOutput$beta0[,spp.ind]
    dhi.rmx <- mod$mcmcOutput$delta0[,spp.ind]
    for(i in 1:length(spp.ind)) {
      bhi.rmx[,i] <- bhi.rmx[,i] + apply(mod$mcmcOutput$beta1[,spp.ind[i],-ind.bx] *
                                   X.psi.pred.hi[,-ind.bx], 1, sum) +
        mod$mcmcOutput$beta1[,spp.ind[i],ind.bx] * X.psi.pred.lo[,ind.bx]
      dhi.rmx[,i] <- dhi.rmx[,i] + apply(mod$mcmcOutput$delta1[,spp.ind[i],-ind.dx] *
                                   X.lambda.pred.hi[,-ind.dx], 1, sum) +
        mod$mcmcOutput$delta1[,spp.ind[i],ind.dx] * X.lambda.pred.lo[,ind.dx]
    }
    psi1 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[1]) *
      QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[1]) *
      mod$mcmcOutput$w[, spp.ind]
    psi10 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[10]) *
      QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[10]) *
      mod$mcmcOutput$w[, spp.ind]
    trend_hi.rmx <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
    diff_lo.rmx <- trend_hi.rmx - trend_lo
    out[g, "Well_1km_contrib_dlo"] <- QSLpersonal::BCI(((diff_lo - diff_lo.rmx) / diff_lo) * 100, BCIpercent = 80)
    out[g, "Well_1km_pcont_dlo"] <- round(sum(((diff_lo - diff_lo.rmx) / diff_lo) > 0) / nsims, digits = 2)
  } else {
    out[g, "Well_1km_contrib_dlo"] <- "not considered"
    out[g, "Well_1km_pcont_dlo"] <-  "not considered"
  }
  if(out[g, "Diff_pbg"] >= 0.9) {
    # Well density contribution to difference from background
    ind.bx <- which(dimnames(X.psi)[[2]] == "Well_1km")
    ind.dx <- which(dimnames(X.lambda)[[2]] == "Well_1km")
    bhi.rmx <- mod$mcmcOutput$beta0[,spp.ind]
    dhi.rmx <- mod$mcmcOutput$delta0[,spp.ind]
    for(i in 1:length(spp.ind)) {
      bhi.rmx[,i] <- bhi.rmx[,i] + apply(mod$mcmcOutput$beta1[,spp.ind[i],-ind.bx] *
                                           X.psi.pred.hi[,-ind.bx], 1, sum) +
        mod$mcmcOutput$beta1[,spp.ind[i],ind.bx] * X.psi.pred.bg[,ind.bx]
      dhi.rmx[,i] <- dhi.rmx[,i] + apply(mod$mcmcOutput$delta1[,spp.ind[i],-ind.dx] *
                                           X.lambda.pred.hi[,-ind.dx], 1, sum) +
        mod$mcmcOutput$delta1[,spp.ind[i],ind.dx] * X.lambda.pred.bg[,ind.dx]
    }
    psi1 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[1]) *
      QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[1]) *
      mod$mcmcOutput$w[, spp.ind]
    psi10 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[10]) *
      QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[10]) *
      mod$mcmcOutput$w[, spp.ind]
    trend_hi.rmx <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
    diff_bg.rmx <- trend_hi.rmx - trend_bg
    out[g, "Well_1km_contrib_dbg"] <- QSLpersonal::BCI(((diff_bg - diff_bg.rmx) / diff_bg) * 100, BCIpercent = 80)
    out[g, "Well_1km_pcont_dbg"] <- round(sum(((diff_bg - diff_bg.rmx) / diff_bg) > 0) / nsims, digits = 2)

    # Road density contribution to difference from background
    ind.bx <- which(dimnames(X.psi)[[2]] == "Road_125m")
    ind.dx <- which(dimnames(X.lambda)[[2]] == "Road_125m")
    bhi.rmx <- mod$mcmcOutput$beta0[,spp.ind]
    dhi.rmx <- mod$mcmcOutput$delta0[,spp.ind]
    for(i in 1:length(spp.ind)) {
      bhi.rmx[,i] <- bhi.rmx[,i] + apply(mod$mcmcOutput$beta1[,spp.ind[i],-ind.bx] *
                                           X.psi.pred.hi[,-ind.bx], 1, sum) +
        mod$mcmcOutput$beta1[,spp.ind[i],ind.bx] * X.psi.pred.bg[,ind.bx]
      dhi.rmx[,i] <- dhi.rmx[,i] + apply(mod$mcmcOutput$delta1[,spp.ind[i],-ind.dx] *
                                           X.lambda.pred.hi[,-ind.dx], 1, sum) +
        mod$mcmcOutput$delta1[,spp.ind[i],ind.dx] * X.lambda.pred.bg[,ind.dx]
    }
    psi1 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[1]) *
      QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[1]) *
      mod$mcmcOutput$w[, spp.ind]
    psi10 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[10]) *
      QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[10]) *
      mod$mcmcOutput$w[, spp.ind]
    trend_hi.rmx <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
    diff_bg.rmx <- trend_hi.rmx - trend_bg
    out[g, "Road_125m_contrib_dbg"] <- QSLpersonal::BCI(((diff_bg - diff_bg.rmx) / diff_bg) * 100, BCIpercent = 80)
    out[g, "Road_125m_pcont_dbg"] <- round(sum(((diff_bg - diff_bg.rmx) / diff_bg) > 0) / nsims, digits = 2)
    
    # Well density contribution to difference from background
    ind.bx <- which(dimnames(X.psi)[[2]] == "AHerb")
    ind.dx <- which(dimnames(X.lambda)[[2]] == "AHerb")
    bhi.rmx <- mod$mcmcOutput$beta0[,spp.ind]
    dhi.rmx <- mod$mcmcOutput$delta0[,spp.ind]
    for(i in 1:length(spp.ind)) {
      bhi.rmx[,i] <- bhi.rmx[,i] + apply(mod$mcmcOutput$beta1[,spp.ind[i],-ind.bx] *
                                           X.psi.pred.hi[,-ind.bx], 1, sum) +
        mod$mcmcOutput$beta1[,spp.ind[i],ind.bx] * X.psi.pred.bg[,ind.bx]
      dhi.rmx[,i] <- dhi.rmx[,i] + apply(mod$mcmcOutput$delta1[,spp.ind[i],-ind.dx] *
                                           X.lambda.pred.hi[,-ind.dx], 1, sum) +
        mod$mcmcOutput$delta1[,spp.ind[i],ind.dx] * X.lambda.pred.bg[,ind.dx]
    }
    psi1 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[1]) *
      QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[1]) *
      mod$mcmcOutput$w[, spp.ind]
    psi10 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[10]) *
      QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[10]) *
      mod$mcmcOutput$w[, spp.ind]
    trend_hi.rmx <- apply(psi10, 1, sum) / apply(psi1, 1, sum)
    diff_bg.rmx <- trend_hi.rmx - trend_bg
    out[g, "AHerb_contrib_dbg"] <- QSLpersonal::BCI(((diff_bg - diff_bg.rmx) / diff_bg) * 100, BCIpercent = 80)
    out[g, "AHerb_pcont_dbg"] <- round(sum(((diff_bg - diff_bg.rmx) / diff_bg) > 0) / nsims, digits = 2)
  } else {
    out[g, "Well_1km_contrib_dbg"] <- "not considered"
    out[g, "Well_1km_pcont_dbg"] <-  "not considered"
    out[g, "Road_125m_contrib_dbg"] <- "not considered"
    out[g, "Road_125m_pcont_dbg"] <-  "not considered"
    out[g, "AHerb_contrib_dbg"] <- "not considered"
    out[g, "AHerb_pcont_dbg"] <-  "not considered"
  }
}

write.csv(out, "Community_mechanisms.csv", row.names = T)
