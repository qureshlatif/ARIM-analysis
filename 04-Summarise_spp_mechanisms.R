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

spp.mech <- c("CONI", "KILL", "HOLA", "ROWR", "HOWR", "SATH", "AMRO", "BRSP",
              "SABS", "WEME", "BHCO")
cols <- c("Diff_lo_tot", "Diff_bg_tot",
          "Well_1km_supp",
          "Well_1km_lo_effect", "Well_1km_lo_contrib",
          "Well_1km_bg_effect", "Well_1km_bg_contrib",
          "Road_125m_supp",
          "Road_125m_bg_effect", "Road_125m_bg_contrib",
          "AHerb_supp",
          "AHerb_bg_effect", "AHerb_bg_contrib")
out <- matrix("", nrow = length(spp.mech), ncol = length(cols),
              dimnames = list(spp.mech, cols))

for(spp in spp.mech) {
  spp.ind <- which(spp.list == spp)

  # Total differences in trend #
  beta_hi <- mod$mcmcOutput$beta0[,spp.ind] +
    apply(mod$mcmcOutput$beta1[,spp.ind,] * X.psi.pred.hi, 1, sum)
  delta_hi <- mod$mcmcOutput$delta0[,spp.ind] +
    apply(mod$mcmcOutput$delta1[,spp.ind,] * X.lambda.pred.hi, 1, sum)
  psi1 <- QSLpersonal::expit(beta_hi + delta_hi * X.trend[1])
  psi10 <- QSLpersonal::expit(beta_hi + delta_hi * X.trend[10])
  trend_hi <- psi10 / psi1
  
  beta_lo <- mod$mcmcOutput$beta0[,spp.ind] +
    apply(mod$mcmcOutput$beta1[,spp.ind,] * X.psi.pred.lo, 1, sum)
  delta_lo <- mod$mcmcOutput$delta0[,spp.ind] +
    apply(mod$mcmcOutput$delta1[,spp.ind,] * X.lambda.pred.lo, 1, sum)
  psi1 <- QSLpersonal::expit(beta_lo + delta_lo * X.trend[1])
  psi10 <- QSLpersonal::expit(beta_lo + delta_lo * X.trend[10])
  trend_lo <- psi10 / psi1
  diff_lo <- trend_hi - trend_lo
  out[spp, "Diff_lo_tot"] <- QSLpersonal::BCI(diff_lo)
  Diff_lo_supp <- str_detect(out[spp, "Diff_lo_tot"], "\\)\\*")
  
  beta_bg <- mod$mcmcOutput$beta0[,spp.ind] +
    apply(mod$mcmcOutput$beta1[,spp.ind,] * X.psi.pred.bg, 1, sum)
  delta_bg <- mod$mcmcOutput$delta0[,spp.ind] +
    apply(mod$mcmcOutput$delta1[,spp.ind,] * X.lambda.pred.bg, 1, sum)
  psi1 <- QSLpersonal::expit(beta_bg + delta_bg * X.trend[1])
  psi10 <- QSLpersonal::expit(beta_bg + delta_bg * X.trend[10])
  trend_bg <- psi10 / psi1
  diff_bg <- trend_hi - trend_bg
  out[spp, "Diff_bg_tot"] <- QSLpersonal::BCI(diff_bg)
  Diff_bg_supp <- str_detect(out[spp, "Diff_bg_tot"], "\\)\\*")
  
  # Evidence for Well_1km mechanism #
    # Consistency of relationship with mechanism #
  ind.x.delta <- which(dimnames(X.lambda)[[2]] == "Well_1km")
  delta.x <- mod$mcmcOutput$delta1[,spp.ind,ind.x.delta]
  if(quantile(delta.x, prob = 0.975, type = 8) < 0) {
    out[spp, "Well_1km_supp"] <- "supported"
  } else {
    if(quantile(delta.x, prob = 0.025, type = 8) > 0) {
      out[spp, "Well_1km_supp"] <- "contradicted"
    } else {
      out[spp, "Well_1km_supp"] <- "unclear"
    }
  }
  
  if(out[spp, "Well_1km_supp"] == "supported") {
    if(Diff_lo_supp) {
      ind.beta.x <- which(dimnames(X.psi)[[2]] == "Well_1km")
      ind.delta.x <- which(dimnames(X.lambda)[[2]] == "Well_1km")
      beta_hi.rmx <- mod$mcmcOutput$beta0[,spp.ind] +
        apply(mod$mcmcOutput$beta1[,spp.ind,-ind.beta.x] * X.psi.pred.hi[,-ind.beta.x], 1, sum) +
        mod$mcmcOutput$beta1[,spp.ind,ind.beta.x] * X.psi.pred.lo[,ind.beta.x]
      delta_hi.rmx <- mod$mcmcOutput$delta0[,spp.ind] +
        apply(mod$mcmcOutput$delta1[,spp.ind,-ind.delta.x] * X.lambda.pred.hi[,-ind.delta.x], 1, sum) +
        mod$mcmcOutput$delta1[,spp.ind,ind.delta.x] * X.lambda.pred.lo[,ind.delta.x]
      psi1 <- QSLpersonal::expit(beta_hi.rmx + delta_hi.rmx * X.trend[1])
      psi10 <- QSLpersonal::expit(beta_hi.rmx + delta_hi.rmx * X.trend[10])
      trend_hi.rmx <- psi10 / psi1
      diff_lo.rmx <- trend_hi.rmx - trend_lo
      out[spp, "Well_1km_lo_effect"] <- QSLpersonal::BCI(diff_lo - diff_lo.rmx)
      out[spp, "Well_1km_lo_contrib"] <- QSLpersonal::BCI(((diff_lo - diff_lo.rmx) / diff_lo) * 100)
    }
    if(Diff_bg_supp) {
      ind.beta.x <- which(dimnames(X.psi)[[2]] == "Well_1km")
      ind.delta.x <- which(dimnames(X.lambda)[[2]] == "Well_1km")
      beta_hi.rmx <- mod$mcmcOutput$beta0[,spp.ind] +
        apply(mod$mcmcOutput$beta1[,spp.ind,-ind.beta.x] * X.psi.pred.hi[,-ind.beta.x], 1, sum) +
        mod$mcmcOutput$beta1[,spp.ind,ind.beta.x] * X.psi.pred.bg[,ind.beta.x]
      delta_hi.rmx <- mod$mcmcOutput$delta0[,spp.ind] +
        apply(mod$mcmcOutput$delta1[,spp.ind,-ind.delta.x] * X.lambda.pred.hi[,-ind.delta.x], 1, sum) +
        mod$mcmcOutput$delta1[,spp.ind,ind.delta.x] * X.lambda.pred.bg[,ind.delta.x]
      psi1 <- QSLpersonal::expit(beta_hi.rmx + delta_hi.rmx * X.trend[1])
      psi10 <- QSLpersonal::expit(beta_hi.rmx + delta_hi.rmx * X.trend[10])
      trend_hi.rmx <- psi10 / psi1
      diff_bg.rmx <- trend_hi.rmx - trend_bg
      out[spp, "Well_1km_bg_effect"] <- QSLpersonal::BCI(diff_bg - diff_bg.rmx)
      out[spp, "Well_1km_bg_contrib"] <- QSLpersonal::BCI(((diff_bg - diff_bg.rmx) / diff_bg) * 100)
    }
  }

  # Evidence for Road_125m mechanism #
    # Consistency of relationship with mechanism #
  ind.x.delta <- which(dimnames(X.lambda)[[2]] == "Road_125m")
  delta.x <- mod$mcmcOutput$delta1[,spp.ind,ind.x.delta]
  if(quantile(delta.x, prob = 0.975, type = 8) < 0) {
    out[spp, "Road_125m_supp"] <- "supported"
  } else {
    if(quantile(delta.x, prob = 0.025, type = 8) > 0) {
      out[spp, "Road_125m_supp"] <- "contradicted"
    } else {
      out[spp, "Road_125m_supp"] <- "unclear"
    }
  }
  
  if(out[spp, "Road_125m_supp"] == "supported") {
    if(Diff_bg_supp) {
      ind.beta.x <- which(dimnames(X.psi)[[2]] == "Road_125m")
      ind.delta.x <- which(dimnames(X.lambda)[[2]] == "Road_125m")
      beta_hi.rmx <- mod$mcmcOutput$beta0[,spp.ind] +
        apply(mod$mcmcOutput$beta1[,spp.ind,-ind.beta.x] * X.psi.pred.hi[,-ind.beta.x], 1, sum) +
        mod$mcmcOutput$beta1[,spp.ind,ind.beta.x] * X.psi.pred.bg[,ind.beta.x]
      delta_hi.rmx <- mod$mcmcOutput$delta0[,spp.ind] +
        apply(mod$mcmcOutput$delta1[,spp.ind,-ind.delta.x] * X.lambda.pred.hi[,-ind.delta.x], 1, sum) +
        mod$mcmcOutput$delta1[,spp.ind,ind.delta.x] * X.lambda.pred.bg[,ind.delta.x]
      psi1 <- QSLpersonal::expit(beta_hi.rmx + delta_hi.rmx * X.trend[1])
      psi10 <- QSLpersonal::expit(beta_hi.rmx + delta_hi.rmx * X.trend[10])
      trend_hi.rmx <- psi10 / psi1
      diff_bg.rmx <- trend_hi.rmx - trend_bg
      out[spp, "Road_125m_bg_effect"] <- QSLpersonal::BCI(diff_bg - diff_bg.rmx)
      out[spp, "Road_125m_bg_contrib"] <- QSLpersonal::BCI(((diff_bg - diff_bg.rmx) / diff_bg) * 100)
    }
  }

  # Evidence for AHerb mechanism #
    # Consistency of relationship with mechanism #
  ind.x.delta <- which(dimnames(X.lambda)[[2]] == "AHerb")
  delta.x <- mod$mcmcOutput$delta1[,spp.ind,ind.x.delta]
  if(quantile(delta.x, prob = 0.975, type = 8) < 0) {
    out[spp, "AHerb_supp"] <- "supported"
  } else {
    if(quantile(delta.x, prob = 0.025, type = 8) > 0) {
      out[spp, "AHerb_supp"] <- "contradicted"
    } else {
      out[spp, "AHerb_supp"] <- "unclear"
    }
  }
  
  if(out[spp, "AHerb_supp"] == "supported") {
    if(Diff_bg_supp) {
      ind.beta.x <- which(dimnames(X.psi)[[2]] == "AHerb")
      ind.delta.x <- which(dimnames(X.lambda)[[2]] == "AHerb")
      beta_hi.rmx <- mod$mcmcOutput$beta0[,spp.ind] +
        apply(mod$mcmcOutput$beta1[,spp.ind,-ind.beta.x] * X.psi.pred.hi[,-ind.beta.x], 1, sum) +
        mod$mcmcOutput$beta1[,spp.ind,ind.beta.x] * X.psi.pred.bg[,ind.beta.x]
      delta_hi.rmx <- mod$mcmcOutput$delta0[,spp.ind] +
        apply(mod$mcmcOutput$delta1[,spp.ind,-ind.delta.x] * X.lambda.pred.hi[,-ind.delta.x], 1, sum) +
        mod$mcmcOutput$delta1[,spp.ind,ind.delta.x] * X.lambda.pred.bg[,ind.delta.x]
      psi1 <- QSLpersonal::expit(beta_hi.rmx + delta_hi.rmx * X.trend[1])
      psi10 <- QSLpersonal::expit(beta_hi.rmx + delta_hi.rmx * X.trend[10])
      trend_hi.rmx <- psi10 / psi1
      diff_bg.rmx <- trend_hi.rmx - trend_bg
      out[spp, "AHerb_bg_effect"] <- QSLpersonal::BCI(diff_bg - diff_bg.rmx)
      out[spp, "AHerb_bg_contrib"] <- QSLpersonal::BCI(((diff_bg - diff_bg.rmx) / diff_bg) * 100)
    }
  }
}

write.csv(out, "Spp_neg_trend_mechanisms.csv", row.names = T)
