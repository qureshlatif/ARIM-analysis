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
source(str_c(scripts.loc, "Data_processing.R"))
#______________________________________#

## Set up output table ##
cols <- c("PSI_hi", "PSI_lo", "PSI_bg", # Needs to include differences in mechanistic covariate conditions.
          "PSIRatioTot_lo",  "PSIRatioTot_bg", "PSIRatioPart_lo", "PSIRatioPart_bg",
          "TREND_hi", "TREND_lo", "TREND_bg",  # Needs to include differences in mechanistic covariate conditions.
          "DIFFTot_lo",  "DIFFTot_bg", "DIFFPart_lo", "DIFFPart_bg",
          "psi_hi", "psi_lo", "psi_bg", # Needs to include differences in mechanistic covariate conditions.
          "psiRatioTot_lo",  "psiRatioTot_bg", "psiRatioPart_lo", "psiRatioPart_bg",
          "trend_hi", "trend_lo", "trend_bg",  # Needs to include differences in mechanistic covariate conditions.
          "diffTot_lo",  "diffTot_bg", "diffPart_lo", "diffPart_bg")
out <- matrix("", nrow = length(spp.list), ncol = length(cols),
              dimnames = list(spp.list, cols))

## Stratum covariate values ##
source(str_c(scripts.loc, "Calculate_mech_path_covariate_values.R"))

## Grid-cell occupancy ##
BETA_hi <- mod$mcmcOutput$BETA0
for(sp in 1:n.spp) BETA_hi[,sp] <- BETA_hi[,sp] +
  apply(mod$mcmcOutput$BETA1[,sp,] * X.PSI.pred.hi, 1, sum)
out[, "PSI_hi"] <-
  QSLpersonal::expit(BETA_hi) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)

BETA_lo <- mod$mcmcOutput$BETA0
for(sp in 1:n.spp) BETA_lo[,sp] <- BETA_lo[,sp] +
  apply(mod$mcmcOutput$BETA1[,sp,] * X.PSI.pred.lo, 1, sum)
out[, "PSI_lo"] <-
  QSLpersonal::expit(BETA_lo) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)

BETA_bg <- mod$mcmcOutput$BETA0
for(sp in 1:n.spp) BETA_bg[,sp] <- BETA_bg[,sp] +
  apply(mod$mcmcOutput$BETA1[,sp,] * X.PSI.pred.bg, 1, sum)
out[, "PSI_bg"] <-
  QSLpersonal::expit(BETA_bg) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)

## Grid-cell occupancy effects ##
ind.strat <- which(dimnames(X.PSI)[[2]] %in% c("Dev_lo", "Dev_bg"))
BETA_hi_partial <- mod$mcmcOutput$BETA0
for(sp in 1:n.spp) BETA_hi_partial[,sp] <- BETA_hi_partial[,sp] +
  apply(mod$mcmcOutput$BETA1[,sp,ind.strat] * X.PSI.pred.hi[,ind.strat], 1, sum)

Ratio <- QSLpersonal::expit(BETA_hi) /
  QSLpersonal::expit(BETA_lo)
out[, "PSIRatioTot_lo"] <- (Ratio) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)
ind.sig <- which(apply(Ratio, 2, quantile, prob = 0.025, type = 8) > 1 |
                   apply(Ratio, 2, quantile, prob = 0.975, type = 8) < 1)
if(length(ind.sig) > 0) out[ind.sig, "PSIRatioTot_lo"] <-
  str_c(out[ind.sig, "PSIRatioTot_lo"], "*")

BETA_lo_partial <- mod$mcmcOutput$BETA0
for(sp in 1:n.spp) BETA_lo_partial[,sp] <- BETA_lo_partial[,sp] +
  apply(mod$mcmcOutput$BETA1[,sp,ind.strat] * X.PSI.pred.lo[,ind.strat], 1, sum)
Ratio <- QSLpersonal::expit(BETA_hi_partial) /
  QSLpersonal::expit(BETA_lo_partial)
out[, "PSIRatioPart_lo"] <- (Ratio) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)
ind.sig <- which(apply(Ratio, 2, quantile, prob = 0.025, type = 8) > 1 |
                   apply(Ratio, 2, quantile, prob = 0.975, type = 8) < 1)
if(length(ind.sig) > 0) out[ind.sig, "PSIRatioPart_lo"] <-
  str_c(out[ind.sig, "PSIRatioPart_lo"], "*")

Ratio <- QSLpersonal::expit(BETA_hi) /
  QSLpersonal::expit(BETA_bg)
out[, "PSIRatioTot_bg"] <- (Ratio) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)
ind.sig <- which(apply(Ratio, 2, quantile, prob = 0.025, type = 8) > 1 |
                   apply(Ratio, 2, quantile, prob = 0.975, type = 8) < 1)
if(length(ind.sig) > 0) out[ind.sig, "PSIRatioTot_bg"] <-
  str_c(out[ind.sig, "PSIRatioTot_bg"], "*")

BETA_bg_partial <- mod$mcmcOutput$BETA0
for(sp in 1:n.spp) BETA_bg_partial[,sp] <- BETA_bg_partial[,sp] +
  apply(mod$mcmcOutput$BETA1[,sp,ind.strat] * X.PSI.pred.bg[,ind.strat], 1, sum)
Ratio <- QSLpersonal::expit(BETA_hi_partial) /
  QSLpersonal::expit(BETA_bg_partial)
out[, "PSIRatioPart_bg"] <- (Ratio) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)
ind.sig <- which(apply(Ratio, 2, quantile, prob = 0.025, type = 8) > 1 |
                   apply(Ratio, 2, quantile, prob = 0.975, type = 8) < 1)
if(length(ind.sig) > 0) out[ind.sig, "PSIRatioPart_bg"] <-
  str_c(out[ind.sig, "PSIRatioPart_bg"], "*")

## Grid-cell trend ##
L <- mod$mcmcOutput$DELTA0
for(sp in 1:n.spp) L[,sp] <- L[,sp] +
  apply(mod$mcmcOutput$DELTA1[,sp,] * X.LAMBDA.pred.hi, 1, sum)
PSI1 <- QSLpersonal::expit(BETA_hi + L * X.trend[1])
PSI10 <- QSLpersonal::expit(BETA_hi + L * X.trend[10])
TREND_hi <- PSI10 / PSI1
out[, "TREND_hi"] <- (TREND_hi) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)
ind.flag <- which(apply(TREND_hi, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 1 |
                    apply(TREND_hi, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 1)
if(length(ind.flag) > 0) out[ind.flag, "TREND_hi"] <- str_c(out[ind.flag, "TREND_hi"], "*")
ind.flag <- which(apply(TREND_hi, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0.75)
if(length(ind.flag) > 0) out[ind.flag, "TREND_hi"] <- str_c(out[ind.flag, "TREND_hi"], "*")

L <- mod$mcmcOutput$DELTA0
for(sp in 1:n.spp) L[,sp] <- L[,sp] +
  apply(mod$mcmcOutput$DELTA1[,sp,] * X.LAMBDA.pred.lo, 1, sum)
PSI1 <- QSLpersonal::expit(BETA_lo + L * X.trend[1])
PSI10 <- QSLpersonal::expit(BETA_lo + L * X.trend[10])
TREND_lo <- PSI10 / PSI1
out[, "TREND_lo"] <-(TREND_lo) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)
ind.flag <- which(apply(TREND_lo, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 1 |
                    apply(TREND_lo, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 1)
if(length(ind.flag) > 0) out[ind.flag, "TREND_lo"] <- str_c(out[ind.flag, "TREND_lo"], "*")

L <- mod$mcmcOutput$DELTA0
for(sp in 1:n.spp) L[,sp] <- L[,sp] +
  apply(mod$mcmcOutput$DELTA1[,sp,] * X.LAMBDA.pred.bg, 1, sum)
PSI1 <- QSLpersonal::expit(BETA_bg + L * X.trend[1])
PSI10 <- QSLpersonal::expit(BETA_bg + L * X.trend[10])
TREND_bg <- PSI10 / PSI1
out[, "TREND_bg"] <- (TREND_bg) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)
ind.flag <- which(apply(TREND_bg, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 1 |
                    apply(TREND_bg, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 1)
if(length(ind.flag) > 0) out[ind.flag, "TREND_bg"] <- str_c(out[ind.flag, "TREND_bg"], "*")

## Grid-cell trend differences ##
ind.strat <- which(dimnames(X.LAMBDA)[[2]] %in% c("Dev_lo", "Dev_bg"))

L <- mod$mcmcOutput$DELTA0
for(sp in 1:n.spp) L[,sp] <- L[,sp] +
  apply(mod$mcmcOutput$DELTA1[,sp,ind.strat] * X.LAMBDA.pred.hi[,ind.strat], 1, sum)
PSI1 <- QSLpersonal::expit(BETA_hi_partial + L * X.trend[1])
PSI10 <- QSLpersonal::expit(BETA_hi_partial + L * X.trend[10])
TREND_hi_Partial <- PSI10 / PSI1

out[, "DIFFTot_lo"] <- (TREND_hi - TREND_lo) %>%
  apply(2, QSLpersonal::BCI)
L <- mod$mcmcOutput$DELTA0
for(sp in 1:n.spp) L[,sp] <- L[,sp] +
  apply(mod$mcmcOutput$DELTA1[,sp,ind.strat] * X.LAMBDA.pred.lo[,ind.strat], 1, sum)
PSI1 <- QSLpersonal::expit(BETA_lo_partial + L * X.trend[1])
PSI10 <- QSLpersonal::expit(BETA_lo_partial + L * X.trend[10])
TREND_lo_Partial <- PSI10 / PSI1
out[, "DIFFPart_lo"] <- (TREND_hi_Partial - TREND_lo_Partial) %>%
  apply(2, QSLpersonal::BCI)

out[, "DIFFTot_bg"] <- (TREND_hi - TREND_bg) %>%
  apply(2, QSLpersonal::BCI)
L <- mod$mcmcOutput$DELTA0
for(sp in 1:n.spp) L[,sp] <- L[,sp] +
  apply(mod$mcmcOutput$DELTA1[,sp,ind.strat] * X.LAMBDA.pred.bg[,ind.strat], 1, sum)
PSI1 <- QSLpersonal::expit(BETA_bg_partial + L * X.trend[1])
PSI10 <- QSLpersonal::expit(BETA_bg_partial + L * X.trend[10])
TREND_bg_Partial <- PSI10 / PSI1
out[, "DIFFPart_bg"] <- (TREND_hi_Partial - TREND_bg_Partial) %>%
  apply(2, QSLpersonal::BCI)

## Point-cell occupancy ##
beta_hi <- mod$mcmcOutput$beta0
for(sp in 1:n.spp) beta_hi[,sp] <- beta_hi[,sp] +
  apply(mod$mcmcOutput$beta1[,sp,] * X.psi.pred.hi, 1, sum)
out[, "psi_hi"] <-
  QSLpersonal::expit(beta_hi) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)

beta_lo <- mod$mcmcOutput$beta0
for(sp in 1:n.spp) beta_lo[,sp] <- beta_lo[,sp] +
  apply(mod$mcmcOutput$beta1[,sp,] * X.psi.pred.lo, 1, sum)
out[, "psi_lo"] <-
  QSLpersonal::expit(beta_lo) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)

beta_bg <- mod$mcmcOutput$beta0
for(sp in 1:n.spp) beta_bg[,sp] <- beta_bg[,sp] +
  apply(mod$mcmcOutput$beta1[,sp,] * X.psi.pred.bg, 1, sum)
out[, "psi_bg"] <-
  QSLpersonal::expit(beta_bg) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)

## Point-scale occupancy effects ##
ind.strat <- which(dimnames(X.psi)[[2]] %in% c("Dev_lo", "Dev_bg"))
beta_hi_partial <- mod$mcmcOutput$beta0
for(sp in 1:n.spp) beta_hi_partial[,sp] <- beta_hi_partial[,sp] +
  apply(mod$mcmcOutput$beta1[,sp,ind.strat] * X.psi.pred.hi[,ind.strat], 1, sum)

Ratio <- QSLpersonal::expit(beta_hi) /
  QSLpersonal::expit(beta_lo)
out[, "psiRatioTot_lo"] <- (Ratio) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)
ind.sig <- which(apply(Ratio, 2, quantile, prob = 0.025, type = 8) > 1 |
                   apply(Ratio, 2, quantile, prob = 0.975, type = 8) < 1)
if(length(ind.sig) > 0) out[ind.sig, "psiRatioTot_lo"] <-
  str_c(out[ind.sig, "psiRatioTot_lo"], "*")

beta_lo_partial <- mod$mcmcOutput$beta0
for(sp in 1:n.spp) beta_lo_partial[,sp] <- beta_lo_partial[,sp] +
  apply(mod$mcmcOutput$beta1[,sp,ind.strat] * X.psi.pred.lo[,ind.strat], 1, sum)
Ratio <- QSLpersonal::expit(beta_hi_partial) /
  QSLpersonal::expit(beta_lo_partial)
out[, "psiRatioPart_lo"] <- (Ratio) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)
ind.sig <- which(apply(Ratio, 2, quantile, prob = 0.025, type = 8) > 1 |
                   apply(Ratio, 2, quantile, prob = 0.975, type = 8) < 1)
if(length(ind.sig) > 0) out[ind.sig, "psiRatioPart_lo"] <-
  str_c(out[ind.sig, "psiRatioPart_lo"], "*")

Ratio <- QSLpersonal::expit(beta_hi) /
  QSLpersonal::expit(beta_bg)
out[, "psiRatioTot_bg"] <- (Ratio) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)
ind.sig <- which(apply(Ratio, 2, quantile, prob = 0.025, type = 8) > 1 |
                   apply(Ratio, 2, quantile, prob = 0.975, type = 8) < 1)
if(length(ind.sig) > 0) out[ind.sig, "psiRatioTot_bg"] <-
  str_c(out[ind.sig, "psiRatioTot_bg"], "*")

beta_bg_partial <- mod$mcmcOutput$beta0
for(sp in 1:n.spp) beta_bg_partial[,sp] <- beta_bg_partial[,sp] +
  apply(mod$mcmcOutput$beta1[,sp,ind.strat] * X.psi.pred.bg[,ind.strat], 1, sum)
Ratio <- QSLpersonal::expit(beta_hi_partial) /
  QSLpersonal::expit(beta_bg_partial)
out[, "psiRatioPart_bg"] <- (Ratio) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)
ind.sig <- which(apply(Ratio, 2, quantile, prob = 0.025, type = 8) > 1 |
                   apply(Ratio, 2, quantile, prob = 0.975, type = 8) < 1)
if(length(ind.sig) > 0) out[ind.sig, "psiRatioPart_bg"] <-
  str_c(out[ind.sig, "psiRatioPart_bg"], "*")

## Point-scale trend ##
L <- mod$mcmcOutput$delta0
for(sp in 1:n.spp) L[,sp] <- L[,sp] +
  apply(mod$mcmcOutput$delta1[,sp,] * X.lambda.pred.hi, 1, sum)
psi1 <- QSLpersonal::expit(beta_hi + L * X.trend[1])
psi10 <- QSLpersonal::expit(beta_hi + L * X.trend[10])
trend_hi <- psi10 / psi1
out[, "trend_hi"] <- (trend_hi) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)
ind.flag <- which(apply(trend_hi, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 1 |
                    apply(trend_hi, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 1)
if(length(ind.flag) > 0) out[ind.flag, "trend_hi"] <- str_c(out[ind.flag, "trend_hi"], "*")
ind.flag <- which(apply(trend_hi, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 0.9)
if(length(ind.flag) > 0) out[ind.flag, "trend_hi"] <- str_c(out[ind.flag, "trend_hi"], "*")

L <- mod$mcmcOutput$delta0
for(sp in 1:n.spp) L[,sp] <- L[,sp] +
  apply(mod$mcmcOutput$delta1[,sp,] * X.lambda.pred.lo, 1, sum)
psi1 <- QSLpersonal::expit(beta_lo + L * X.trend[1])
psi10 <- QSLpersonal::expit(beta_lo + L * X.trend[10])
trend_lo <- psi10 / psi1
out[, "trend_lo"] <-(trend_lo) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)
ind.flag <- which(apply(trend_lo, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 1 |
                    apply(trend_lo, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 1)
if(length(ind.flag) > 0) out[ind.flag, "trend_lo"] <- str_c(out[ind.flag, "trend_lo"], "*")

L <- mod$mcmcOutput$delta0
for(sp in 1:n.spp) L[,sp] <- L[,sp] +
  apply(mod$mcmcOutput$delta1[,sp,] * X.lambda.pred.bg, 1, sum)
psi1 <- QSLpersonal::expit(beta_bg + L * X.trend[1])
psi10 <- QSLpersonal::expit(beta_bg + L * X.trend[10])
trend_bg <- psi10 / psi1
out[, "trend_bg"] <- (trend_bg) %>%
  apply(2, QSLpersonal::BCI, flag.sig = F)
ind.flag <- which(apply(trend_bg, 2, function(x) quantile(x, prob = 0.975, type = 8)) < 1 |
                    apply(trend_bg, 2, function(x) quantile(x, prob = 0.025, type = 8)) > 1)
if(length(ind.flag) > 0) out[ind.flag, "trend_bg"] <- str_c(out[ind.flag, "trend_bg"], "*")

## Point-scale trend differences ##
ind.strat <- which(dimnames(X.lambda)[[2]] %in% c("Dev_lo", "Dev_bg"))

L <- mod$mcmcOutput$delta0
for(sp in 1:n.spp) L[,sp] <- L[,sp] +
  apply(mod$mcmcOutput$delta1[,sp,ind.strat] * X.lambda.pred.hi[,ind.strat], 1, sum)
psi1 <- QSLpersonal::expit(beta_hi_partial + L * X.trend[1])
psi10 <- QSLpersonal::expit(beta_hi_partial + L * X.trend[10])
trend_hi_Partial <- psi10 / psi1

out[, "diffTot_lo"] <- (trend_hi - trend_lo) %>%
  apply(2, QSLpersonal::BCI)
L <- mod$mcmcOutput$delta0
for(sp in 1:n.spp) L[,sp] <- L[,sp] +
  apply(mod$mcmcOutput$delta1[,sp,ind.strat] * X.lambda.pred.lo[,ind.strat], 1, sum)
psi1 <- QSLpersonal::expit(beta_lo_partial + L * X.trend[1])
psi10 <- QSLpersonal::expit(beta_lo_partial + L * X.trend[10])
trend_lo_Partial <- psi10 / psi1
out[, "diffPart_lo"] <- (trend_hi_Partial - trend_lo_Partial) %>%
  apply(2, QSLpersonal::BCI)

out[, "diffTot_bg"] <- (trend_hi - trend_bg) %>%
  apply(2, QSLpersonal::BCI)
L <- mod$mcmcOutput$delta0
for(sp in 1:n.spp) L[,sp] <- L[,sp] +
  apply(mod$mcmcOutput$delta1[,sp,ind.strat] * X.lambda.pred.bg[,ind.strat], 1, sum)
psi1 <- QSLpersonal::expit(beta_bg_partial + L * X.trend[1])
psi10 <- QSLpersonal::expit(beta_bg_partial + L * X.trend[10])
trend_bg_Partial <- psi10 / psi1
out[, "diffPart_bg"] <- (trend_hi_Partial - trend_bg_Partial) %>%
  apply(2, QSLpersonal::BCI)

write.csv(out, "Summary_spp_development_effects.csv", row.names = T)

# # Species with any supported stratum effects or trends (partial or total)
# spp.list[which(apply(out, 1, function(x) any(str_detect(x, "\\)\\*"))))]

# # Species with supported total stratum effects or trends (not partial)
# spp.list[which(apply(out[,c("PSIRatioTot_lo",  "PSIRatioTot_bg",
#                             "TREND_hi", "TREND_lo", "TREND_bg",
#                             "DIFFTot_lo",  "DIFFTot_bg",
#                             "psiRatioTot_lo",  "psiRatioTot_bg",
#                             "trend_hi", "trend_lo", "trend_bg",
#                             "diffTot_lo",  "diffTot_bg")],
#                      1, function(x) any(str_detect(x, "\\)\\*"))))]
# 
# # Species with supported total stratum effects or trends at grid level (not partial)
# spp.list[which(apply(out[,c("PSIRatioTot_lo",  "PSIRatioTot_bg",
#                             "TREND_hi", "TREND_lo", "TREND_bg",
#                             "DIFFTot_lo",  "DIFFTot_bg")],
#                      1, function(x) any(str_detect(x, "\\)\\*"))))]

# Species with supported total stratum-level trends or stratum effects on trends (not partial and not on occupancy)
spp.list[which(apply(out[,c("TREND_hi", "TREND_lo",
                            "DIFFTot_lo",  "DIFFTot_bg",
                            "trend_hi", "trend_lo",
                            "diffTot_lo",  "diffTot_bg")],
                     1, function(x) any(str_detect(x, "\\)\\*"))))]
out[which(apply(out[,c("TREND_hi", "TREND_lo",
                   "DIFFTot_lo",  "DIFFTot_bg",
                   "trend_hi", "trend_lo",
                   "diffTot_lo",  "diffTot_bg")],
            1, function(x) any(str_detect(x, "\\)\\*")))),
    c("trend_hi", "trend_lo",
      "diffTot_lo",  "diffTot_bg")] %>%
  View()
