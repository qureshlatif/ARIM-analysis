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

spp.mech <- loadObject("Spp_mechanisms")
cols <- c("DIFF_lo", "DIFF_plo",
          "Well_3km_contrib_Dlo", "Well_3km_pcont_Dlo",
          "DIFF_bg", "DIFF_pbg",
          "Well_3km_contrib_Dbg", "Well_3km_pcont_Dbg",
          "Diff_lo", "Diff_plo",
          "Well_1km_contrib_dlo", "Well_1km_pcont_dlo",
          "AHerb_contrib_dlo", "AHerb_pcont_dlo",
          "Diff_bg", "Diff_pbg",
          "Well_1km_contrib_dbg", "Well_1km_pcont_dbg",
          "Road_125m_contrib_dbg", "Road_125m_pcont_dbg",
          "AHerb_contrib_dbg", "AHerb_pcont_dbg")
out <- matrix("", nrow = length(spp.mech), ncol = length(cols),
              dimnames = list(spp.mech, cols))

for(spp in spp.mech) {
  spp.ind <- which(spp.list == spp)
  
  ## Grid scale ##
  # Differences in trend #
  BETA_hi <- mod$mcmcOutput$BETA0[,spp.ind] +
    apply(mod$mcmcOutput$BETA1[,spp.ind,] * X.PSI.pred.hi, 1, sum)
  DELTA_hi <- mod$mcmcOutput$DELTA0[,spp.ind] +
    apply(mod$mcmcOutput$DELTA1[,spp.ind,] * X.LAMBDA.pred.hi, 1, sum)
  PSI1 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[1])
  PSI10 <- QSLpersonal::expit(BETA_hi + DELTA_hi * X.trend[10])
  TREND_hi <- (PSI10 / (1 - PSI10)) / (PSI1 / (1 - PSI1))
  
  BETA_lo <- mod$mcmcOutput$BETA0[,spp.ind] +
    apply(mod$mcmcOutput$BETA1[,spp.ind,] * X.PSI.pred.lo, 1, sum)
  DELTA_lo <- mod$mcmcOutput$DELTA0[,spp.ind] +
    apply(mod$mcmcOutput$DELTA1[,spp.ind,] * X.LAMBDA.pred.lo, 1, sum)
  PSI1 <- QSLpersonal::expit(BETA_lo + DELTA_lo * X.trend[1])
  PSI10 <- QSLpersonal::expit(BETA_lo + DELTA_lo * X.trend[10])
  TREND_lo <- (PSI10 / (1 - PSI10)) / (PSI1 / (1 - PSI1))
  DIFF_lo <- TREND_hi - TREND_lo
  out[spp, "DIFF_lo"] <- QSLpersonal::BCI(DIFF_lo, BCIpercent = 90)
  out[spp, "DIFF_plo"] <- round(sum(DIFF_lo < 0) / nsims, digits = 2)
  
  BETA_bg <- mod$mcmcOutput$BETA0[,spp.ind] +
    apply(mod$mcmcOutput$BETA1[,spp.ind,] * X.PSI.pred.bg, 1, sum)
  DELTA_bg <- mod$mcmcOutput$DELTA0[,spp.ind] +
    apply(mod$mcmcOutput$DELTA1[,spp.ind,] * X.LAMBDA.pred.bg, 1, sum)
  PSI1 <- QSLpersonal::expit(BETA_bg + DELTA_bg * X.trend[1])
  PSI10 <- QSLpersonal::expit(BETA_bg + DELTA_bg * X.trend[10])
  TREND_bg <- (PSI10 / (1 - PSI10)) / (PSI1 / (1 - PSI1))
  DIFF_bg <- TREND_hi - TREND_bg
  out[spp, "DIFF_bg"] <- QSLpersonal::BCI(DIFF_bg, BCIpercent = 90)
  out[spp, "DIFF_pbg"] <- round(sum(DIFF_bg < 0) / nsims, digits = 2)
  
  # Evidence for mechanisms #
  if(out[spp, "DIFF_plo"] >= 0.9) {
    # Well density contribution to difference from low development
    ind.bx <- which(dimnames(X.PSI)[[2]] == "Well_3km")
    ind.dx <- which(dimnames(X.LAMBDA)[[2]] == "Well_3km")
    bhi.rmx <- mod$mcmcOutput$BETA0[,spp.ind] +
      apply(mod$mcmcOutput$BETA1[,spp.ind,-ind.bx] *
              X.PSI.pred.hi[,-ind.bx], 1, sum) +
      mod$mcmcOutput$BETA1[,spp.ind,ind.bx] * X.PSI.pred.lo[,ind.bx]
    dhi.rmx <- mod$mcmcOutput$DELTA0[,spp.ind] +
      apply(mod$mcmcOutput$DELTA1[,spp.ind,-ind.dx] *
              X.LAMBDA.pred.hi[,-ind.dx], 1, sum) +
      mod$mcmcOutput$DELTA1[,spp.ind,ind.dx] * X.LAMBDA.pred.lo[,ind.dx]
    psi1 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[1])
    psi10 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[10])
    trend_hi.rmx <- (psi10 / (1 - psi10)) / (psi1 / (1 - psi1))
    diff_lo.rmx <- trend_hi.rmx - TREND_lo
    out[spp, "Well_3km_contrib_Dlo"] <- QSLpersonal::BCI(((DIFF_lo - diff_lo.rmx) / DIFF_lo) * 100, BCIpercent = 90)
    out[spp, "Well_3km_pcont_Dlo"] <- round(sum(((DIFF_lo - diff_lo.rmx) / DIFF_lo) > 0) / nsims, digits = 2)
  } else {
    out[spp, "Well_3km_contrib_Dlo"] <- "not considered"
    out[spp, "Well_3km_pcont_Dlo"] <-  "not considered"
  }
  if(out[spp, "DIFF_pbg"] >= 0.9) {
    # Well density contribution to difference from background
    ind.bx <- which(dimnames(X.PSI)[[2]] == "Well_3km")
    ind.dx <- which(dimnames(X.LAMBDA)[[2]] == "Well_3km")
    bhi.rmx <- mod$mcmcOutput$BETA0[,spp.ind] +
      apply(mod$mcmcOutput$BETA1[,spp.ind,-ind.bx] *
              X.PSI.pred.hi[,-ind.bx], 1, sum) +
      mod$mcmcOutput$BETA1[,spp.ind,ind.bx] * X.PSI.pred.bg[,ind.bx]
    dhi.rmx <- mod$mcmcOutput$DELTA0[,spp.ind] +
      apply(mod$mcmcOutput$DELTA1[,spp.ind,-ind.dx] *
              X.LAMBDA.pred.hi[,-ind.dx], 1, sum) +
      mod$mcmcOutput$DELTA1[,spp.ind,ind.dx] * X.LAMBDA.pred.bg[,ind.dx]
    psi1 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[1])
    psi10 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[10])
    trend_hi.rmx <- (psi10 / (1 - psi10)) / (psi1 / (1 - psi1))
    diff_bg.rmx <- trend_hi.rmx - TREND_bg
    out[spp, "Well_3km_contrib_Dbg"] <- QSLpersonal::BCI(((DIFF_bg - diff_bg.rmx) / DIFF_bg) * 100, BCIpercent = 90)
    out[spp, "Well_3km_pcont_Dbg"] <- round(sum(((DIFF_bg - diff_bg.rmx) / DIFF_bg) > 0) / nsims, digits = 2)
  } else {
    out[spp, "Well_3km_contrib_Dbg"] <- "not considered"
    out[spp, "Well_3km_pcont_Dbg"] <-  "not considered"
  }

  ## Point scale ##
  # Differences in trend #
  beta_hi <- mod$mcmcOutput$beta0[,spp.ind] +
    apply(mod$mcmcOutput$beta1[,spp.ind,] * X.psi.pred.hi, 1, sum)
  delta_hi <- mod$mcmcOutput$delta0[,spp.ind] +
    apply(mod$mcmcOutput$delta1[,spp.ind,] * X.lambda.pred.hi, 1, sum)
  psi1 <- QSLpersonal::expit(beta_hi + delta_hi * X.trend[1])
  psi10 <- QSLpersonal::expit(beta_hi + delta_hi * X.trend[10])
  trend_hi <- (psi10 / (1 - psi10)) / (psi1 / (1 - psi1))
  
  beta_lo <- mod$mcmcOutput$beta0[,spp.ind] +
    apply(mod$mcmcOutput$beta1[,spp.ind,] * X.psi.pred.lo, 1, sum)
  delta_lo <- mod$mcmcOutput$delta0[,spp.ind] +
    apply(mod$mcmcOutput$delta1[,spp.ind,] * X.lambda.pred.lo, 1, sum)
  psi1 <- QSLpersonal::expit(beta_lo + delta_lo * X.trend[1])
  psi10 <- QSLpersonal::expit(beta_lo + delta_lo * X.trend[10])
  trend_lo <- (psi10 / (1 - psi10)) / (psi1 / (1 - psi1))
  diff_lo <- trend_hi - trend_lo
  out[spp, "Diff_lo"] <- QSLpersonal::BCI(diff_lo, BCIpercent = 90)
  out[spp, "Diff_plo"] <- round(sum(diff_lo < 0) / nsims, digits = 2)
  
  beta_bg <- mod$mcmcOutput$beta0[,spp.ind] +
    apply(mod$mcmcOutput$beta1[,spp.ind,] * X.psi.pred.bg, 1, sum)
  delta_bg <- mod$mcmcOutput$delta0[,spp.ind] +
    apply(mod$mcmcOutput$delta1[,spp.ind,] * X.lambda.pred.bg, 1, sum)
  psi1 <- QSLpersonal::expit(beta_bg + delta_bg * X.trend[1])
  psi10 <- QSLpersonal::expit(beta_bg + delta_bg * X.trend[10])
  trend_bg <- (psi10 / (1 - psi10)) / (psi1 / (1 - psi1))
  diff_bg <- trend_hi - trend_bg
  out[spp, "Diff_bg"] <- QSLpersonal::BCI(diff_bg, BCIpercent = 90)
  out[spp, "Diff_pbg"] <- round(sum(diff_bg < 0) / nsims, digits = 2)
  
  # Evidence for mechanisms #
  if(out[spp, "Diff_plo"] >= 0.9) {
    # Well density contribution to difference from low development
    ind.bx <- which(dimnames(X.psi)[[2]] == "Well_1km")
    ind.dx <- which(dimnames(X.lambda)[[2]] == "Well_1km")
    bhi.rmx <- mod$mcmcOutput$beta0[,spp.ind] +
      apply(mod$mcmcOutput$beta1[,spp.ind,-ind.bx] *
              X.psi.pred.hi[,-ind.bx], 1, sum) +
      mod$mcmcOutput$beta1[,spp.ind,ind.bx] * X.psi.pred.lo[,ind.bx]
    dhi.rmx <- mod$mcmcOutput$delta0[,spp.ind] +
      apply(mod$mcmcOutput$delta1[,spp.ind,-ind.dx] *
              X.lambda.pred.hi[,-ind.dx], 1, sum) +
      mod$mcmcOutput$delta1[,spp.ind,ind.dx] * X.lambda.pred.lo[,ind.dx]
    psi1 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[1])
    psi10 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[10])
    trend_hi.rmx <- (psi10 / (1 - psi10)) / (psi1 / (1 - psi1))
    diff_lo.rmx <- trend_hi.rmx - trend_lo
    out[spp, "Well_1km_contrib_dlo"] <- QSLpersonal::BCI(((diff_lo - diff_lo.rmx) / diff_lo) * 100, BCIpercent = 90)
    out[spp, "Well_1km_pcont_dlo"] <- round(sum(((diff_lo - diff_lo.rmx) / diff_lo) > 0) / nsims, digits = 2)
    
    # Annual herbaceous contribution to difference from low development
    ind.bx <- which(dimnames(X.psi)[[2]] == "AHerb")
    ind.dx <- which(dimnames(X.lambda)[[2]] == "AHerb")
    bhi.rmx <- mod$mcmcOutput$beta0[,spp.ind] +
      apply(mod$mcmcOutput$beta1[,spp.ind,-ind.bx] *
              X.psi.pred.hi[,-ind.bx], 1, sum) +
      mod$mcmcOutput$beta1[,spp.ind,ind.bx] * X.psi.pred.lo[,ind.bx]
    dhi.rmx <- mod$mcmcOutput$delta0[,spp.ind] +
      apply(mod$mcmcOutput$delta1[,spp.ind,-ind.dx] *
              X.lambda.pred.hi[,-ind.dx], 1, sum) +
      mod$mcmcOutput$delta1[,spp.ind,ind.dx] * X.lambda.pred.lo[,ind.dx]
    psi1 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[1])
    psi10 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[10])
    trend_hi.rmx <- (psi10 / (1 - psi10)) / (psi1 / (1 - psi1))
    diff_lo.rmx <- trend_hi.rmx - trend_lo
    out[spp, "AHerb_contrib_dlo"] <- QSLpersonal::BCI(((diff_lo - diff_lo.rmx) / diff_lo) * 100, BCIpercent = 90)
    out[spp, "AHerb_pcont_dlo"] <- round(sum(((diff_lo - diff_lo.rmx) / diff_lo) > 0) / nsims, digits = 2)
    
    
  } else {
    out[spp, "Well_1km_contrib_dlo"] <- "not considered"
    out[spp, "Well_1km_pcont_dlo"] <-  "not considered"
    out[spp, "AHerb_contrib_dlo"] <- "not considered"
    out[spp, "AHerb_pcont_dlo"] <-  "not considered"
  }
  if(out[spp, "Diff_pbg"] >= 0.9) {
    # Well density contribution to difference from background
    ind.bx <- which(dimnames(X.psi)[[2]] == "Well_1km")
    ind.dx <- which(dimnames(X.lambda)[[2]] == "Well_1km")
    bhi.rmx <- mod$mcmcOutput$beta0[,spp.ind] +
      apply(mod$mcmcOutput$beta1[,spp.ind,-ind.bx] *
              X.psi.pred.hi[,-ind.bx], 1, sum) +
      mod$mcmcOutput$beta1[,spp.ind,ind.bx] * X.psi.pred.bg[,ind.bx]
    dhi.rmx <- mod$mcmcOutput$delta0[,spp.ind] +
      apply(mod$mcmcOutput$delta1[,spp.ind,-ind.dx] *
              X.lambda.pred.hi[,-ind.dx], 1, sum) +
      mod$mcmcOutput$delta1[,spp.ind,ind.dx] * X.lambda.pred.bg[,ind.dx]
    psi1 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[1])
    psi10 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[10])
    trend_hi.rmx <- (psi10 / (1 - psi10)) / (psi1 / (1 - psi1))
    diff_bg.rmx <- trend_hi.rmx - trend_bg
    out[spp, "Well_1km_contrib_dbg"] <- QSLpersonal::BCI(((diff_bg - diff_bg.rmx) / diff_bg) * 100, BCIpercent = 90)
    out[spp, "Well_1km_pcont_dbg"] <- round(sum(((diff_bg - diff_bg.rmx) / diff_bg) > 0) / nsims, digits = 2)
    
    # Road density contribution to difference from background
    ind.bx <- which(dimnames(X.psi)[[2]] == "Road_125m")
    ind.dx <- which(dimnames(X.lambda)[[2]] == "Road_125m")
    bhi.rmx <- mod$mcmcOutput$beta0[,spp.ind] +
      apply(mod$mcmcOutput$beta1[,spp.ind,-ind.bx] *
              X.psi.pred.hi[,-ind.bx], 1, sum) +
      mod$mcmcOutput$beta1[,spp.ind,ind.bx] * X.psi.pred.bg[,ind.bx]
    dhi.rmx <- mod$mcmcOutput$delta0[,spp.ind] +
      apply(mod$mcmcOutput$delta1[,spp.ind,-ind.dx] *
              X.lambda.pred.hi[,-ind.dx], 1, sum) +
      mod$mcmcOutput$delta1[,spp.ind,ind.dx] * X.lambda.pred.bg[,ind.dx]
    psi1 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[1])
    psi10 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[10])
    trend_hi.rmx <- (psi10 / (1 - psi10)) / (psi1 / (1 - psi1))
    diff_bg.rmx <- trend_hi.rmx - trend_bg
    out[spp, "Road_125m_contrib_dbg"] <- QSLpersonal::BCI(((diff_bg - diff_bg.rmx) / diff_bg) * 100, BCIpercent = 90)
    out[spp, "Road_125m_pcont_dbg"] <- round(sum(((diff_bg - diff_bg.rmx) / diff_bg) > 0) / nsims, digits = 2)
    
    # Annual herbaceous contribution to difference from background
    ind.bx <- which(dimnames(X.psi)[[2]] == "AHerb")
    ind.dx <- which(dimnames(X.lambda)[[2]] == "AHerb")
    bhi.rmx <- mod$mcmcOutput$beta0[,spp.ind] +
      apply(mod$mcmcOutput$beta1[,spp.ind,-ind.bx] *
              X.psi.pred.hi[,-ind.bx], 1, sum) +
      mod$mcmcOutput$beta1[,spp.ind,ind.bx] * X.psi.pred.bg[,ind.bx]
    dhi.rmx <- mod$mcmcOutput$delta0[,spp.ind] +
      apply(mod$mcmcOutput$delta1[,spp.ind,-ind.dx] *
              X.lambda.pred.hi[,-ind.dx], 1, sum) +
      mod$mcmcOutput$delta1[,spp.ind,ind.dx] * X.lambda.pred.bg[,ind.dx]
    psi1 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[1])
    psi10 <- QSLpersonal::expit(bhi.rmx + dhi.rmx * X.trend[10])
    trend_hi.rmx <- (psi10 / (1 - psi10)) / (psi1 / (1 - psi1))
    diff_bg.rmx <- trend_hi.rmx - trend_bg
    out[spp, "AHerb_contrib_dbg"] <- QSLpersonal::BCI(((diff_bg - diff_bg.rmx) / diff_bg) * 100, BCIpercent = 90)
    out[spp, "AHerb_pcont_dbg"] <- round(sum(((diff_bg - diff_bg.rmx) / diff_bg) > 0) / nsims, digits = 2)
    
  } else {
    out[spp, "Well_1km_contrib_dbg"] <- "not considered"
    out[spp, "Well_1km_pcont_dbg"] <-  "not considered"
    out[spp, "Road_125m_contrib_dbg"] <- "not considered"
    out[spp, "Road_125m_pcont_dbg"] <-  "not considered"
    out[spp, "AHerb_contrib_dbg"] <- "not considered"
    out[spp, "AHerb_pcont_dbg"] <-  "not considered"
  }
}

write.csv(out, "Spp_neg_trend_mechanisms.csv", row.names = T)
