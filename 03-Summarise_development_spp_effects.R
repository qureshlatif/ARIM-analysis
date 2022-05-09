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
cols <- c("TREND_hi", "p_Tlt100_hi", "p_Tlt75_hi",
          "TREND_lo", "p_Thi_lt_Tlo", "p_Thi_2lt_Tlo",
          "TREND_bg", "p_Thi_lt_Tbg", "p_Thi_2lt_Tbg",
          "trend_hi", "p_tlt100_hi", "p_tlt90_hi",
          "trend_lo", "p_thi_lt_tlo", "p_thi_2lt_tlo",
          "trend_bg", "p_thi_lt_tbg", "p_thi_2lt_tbg")
out <- matrix("", nrow = length(spp.list), ncol = length(cols),
              dimnames = list(spp.list, cols))

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
TREND_hi <- PSI10 / PSI1
out[, "TREND_hi"] <- TREND_hi %>%
  apply(2, QSLpersonal::BCI, flag.sig = F, BCIpercent = 80)
out[,"p_Tlt100_hi"] <- TREND_hi %>%
  apply(2, function(x) sum(x < 1)/nsims) %>%
  round(digits = 2)
out[,"p_Tlt75_hi"] <- TREND_hi %>%
  apply(2, function(x) sum(x < 0.75)/nsims) %>%
  round(digits = 2)

PSI1 <- QSLpersonal::expit(BETA_lo + DELTA_lo * X.trend[1])
PSI10 <- QSLpersonal::expit(BETA_lo + DELTA_lo * X.trend[10])
TREND_lo <- PSI10 / PSI1
out[, "TREND_lo"] <- TREND_lo %>%
  apply(2, QSLpersonal::BCI, flag.sig = F, BCIpercent = 80)
DIFF <- TREND_hi - TREND_lo
out[,"p_Thi_lt_Tlo"] <- DIFF %>%
  apply(2, function(x) sum(x < 0)/nsims) %>%
  round(digits = 2)
out[,"p_Thi_2lt_Tlo"] <- DIFF %>%
  apply(2, function(x) sum(x < -0.02)/nsims) %>%
  round(digits = 2)

PSI1 <- QSLpersonal::expit(BETA_bg + DELTA_bg * X.trend[1])
PSI10 <- QSLpersonal::expit(BETA_bg + DELTA_bg * X.trend[10])
TREND_bg <- PSI10 / PSI1
out[, "TREND_bg"] <- TREND_bg %>%
  apply(2, QSLpersonal::BCI, flag.sig = F, BCIpercent = 80)
DIFF <- TREND_hi - TREND_bg
out[,"p_Thi_lt_Tbg"] <- DIFF %>%
  apply(2, function(x) sum(x < 0)/nsims) %>%
  round(digits = 2)
out[,"p_Thi_2lt_Tbg"] <- DIFF %>%
  apply(2, function(x) sum(x < -0.02)/nsims) %>%
  round(digits = 2)

## Point-scale trend ##
psi1 <- QSLpersonal::expit(beta_hi + delta_hi * X.trend[1])
psi10 <- QSLpersonal::expit(beta_hi + delta_hi * X.trend[10])
trend_hi <- psi10 / psi1
out[, "trend_hi"] <- trend_hi %>%
  apply(2, QSLpersonal::BCI, flag.sig = F, BCIpercent = 80)
out[,"p_tlt100_hi"] <- trend_hi %>%
  apply(2, function(x) sum(x < 1)/nsims) %>%
  round(digits = 2)
out[,"p_tlt90_hi"] <- trend_hi %>%
  apply(2, function(x) sum(x < 0.9)/nsims) %>%
  round(digits = 2)

psi1 <- QSLpersonal::expit(beta_lo + delta_lo * X.trend[1])
psi10 <- QSLpersonal::expit(beta_lo + delta_lo * X.trend[10])
trend_lo <- psi10 / psi1
out[, "trend_lo"] <- trend_lo %>%
  apply(2, QSLpersonal::BCI, flag.sig = F, BCIpercent = 80)
DIFF <- trend_hi - trend_lo
out[,"p_thi_lt_tlo"] <- DIFF %>%
  apply(2, function(x) sum(x < 0)/nsims) %>%
  round(digits = 2)
out[,"p_thi_2lt_tlo"] <- DIFF %>%
  apply(2, function(x) sum(x < -0.02)/nsims) %>%
  round(digits = 2)

psi1 <- QSLpersonal::expit(beta_bg + delta_bg * X.trend[1])
psi10 <- QSLpersonal::expit(beta_bg + delta_bg * X.trend[10])
trend_bg <- psi10 / psi1
out[, "trend_bg"] <- trend_bg %>%
  apply(2, QSLpersonal::BCI, flag.sig = F, BCIpercent = 80)
DIFF <- trend_hi - trend_bg
out[,"p_thi_lt_tbg"] <- DIFF %>%
  apply(2, function(x) sum(x < 0)/nsims) %>%
  round(digits = 2)
out[,"p_thi_2lt_tbg"] <- DIFF %>%
  apply(2, function(x) sum(x < -0.02)/nsims) %>%
  round(digits = 2)

write.csv(out, "Summary_spp_development_effects.csv", row.names = T)

## Species lists ##
out <- out %>% data.frame(stringsAsFactors = F) %>%
  mutate(Spp = row.names(out),
         p_Tlt100_hi = as.numeric(p_Tlt100_hi),
         p_Tlt75_hi = as.numeric(p_Tlt75_hi),
         p_Thi_lt_Tlo = as.numeric(p_Thi_lt_Tlo),
         p_Thi_2lt_Tlo = as.numeric(p_Thi_2lt_Tlo),
         p_Thi_lt_Tbg = as.numeric(p_Thi_lt_Tbg),
         p_Thi_2lt_Tbg = as.numeric(p_Thi_2lt_Tbg),
         p_tlt100_hi = as.numeric(p_tlt100_hi),
         p_tlt90_hi = as.numeric(p_tlt90_hi),
         p_thi_lt_tlo = as.numeric(p_thi_lt_tlo),
         p_thi_2lt_tlo = as.numeric(p_thi_2lt_tlo),
         p_thi_lt_tbg = as.numeric(p_thi_lt_tbg),
         p_thi_2lt_tbg = as.numeric(p_thi_2lt_tbg))

# Species with supported negative high-development trend
out %>%
  filter(p_Tlt100_hi >= 0.9 |
           p_tlt100_hi >= 0.9) %>%
  View()

# Species with supported stratum effects on trend
out %>%
  filter(p_Thi_lt_Tlo >= 0.9 | p_Thi_lt_Tlo <= 0.1 |
           p_Thi_lt_Tbg >= 0.9 | p_Thi_lt_Tbg <= 0.1 |
           p_thi_lt_tlo >= 0.9 | p_thi_lt_tlo <= 0.1 |
           p_thi_lt_tbg >= 0.9 | p_thi_lt_tbg <= 0.1) %>%
  left_join(spp.out %>% select(BirdCode, Guild),
            by = c("Spp" = "BirdCode")) %>%
  write.csv("Summary_spp_development_effects.csv", row.names = T)
out %>%
  filter(p_Thi_lt_Tlo >= 0.9 |
           p_Thi_lt_Tbg >= 0.9 |
           p_thi_lt_tlo >= 0.9 |
           p_thi_lt_tbg >= 0.9) %>%
  pull(Spp) %>%
  saveObject("Spp_mechanisms")

# Species with strongly supported negative high-dev effect on trend
out %>%
  filter(p_Thi_lt_Tlo >= 0.95 |
           p_Thi_lt_Tbg >= 0.95 |
           p_thi_lt_tlo >= 0.95 |
           p_thi_lt_tbg >= 0.95) %>%
  View()
out %>%
  filter(p_Thi_lt_Tlo >= 0.95 |
           p_Thi_lt_Tbg >= 0.95 |
           p_thi_lt_tlo >= 0.95 |
           p_thi_lt_tbg >= 0.95) %>%
  pull(Spp) %>%
  saveObject("Spp_plot_trends")

# Species with supported negative high-development trend and negative high-dev effect on trend
out %>%
  filter((p_Tlt100_hi >= 0.9 &
            (p_Thi_lt_Tlo >= 0.9 |
               p_Thi_lt_Tbg >= 0.9)) |
           (p_tlt100_hi >= 0.9 &
              (p_thi_lt_tlo >= 0.9 |
                 p_thi_lt_tbg >= 0.9))) %>%
  View()

# Species with considerable evidence for triggers
out %>%
  filter((p_Tlt75_hi >= 0.9 &
            (p_Thi_2lt_Tlo >= 0.9 |
              p_Thi_2lt_Tbg >= 0.9)) |
           (p_tlt90_hi >= 0.9 &
              (p_thi_2lt_tlo >= 0.9 |
                 p_thi_2lt_tbg >= 0.9))) %>%
  View()
