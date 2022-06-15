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
out.vals <- c("est", "f")
#______________________________________#

mod$summary %>% write.csv(str_c("Params_all.csv"), row.names = F) # Save full table of parameter estimates
mod$summary %>%
  filter((!str_detect(Parameter, "BETA0\\[") &
           !str_detect(Parameter, "BETA1\\[") &
           !str_detect(Parameter, "dev.BETA\\[") &
           !str_detect(Parameter, "DELTA0\\[") &
           !str_detect(Parameter, "DELTA1\\[") &
           !str_detect(Parameter, "beta0\\[") &
           !str_detect(Parameter, "beta1\\[") &
           !str_detect(Parameter, "dev.beta\\[") &
           !str_detect(Parameter, "delta0\\[") &
           !str_detect(Parameter, "delta1\\[") &
           !str_detect(Parameter, "zeta0\\[") &
           !str_detect(Parameter, "zeta1\\[") &
           !str_detect(Parameter, "dev.zeta\\[") &
           !str_detect(Parameter, "w\\[")) |
           str_detect(Parameter, "sigma")) %>%
  write.csv(str_c("HyperParams.csv"), row.names = F)

# Compile data to get covariate names #
source(str_c(scripts.loc, "Data_processing.R"))

params <- c("PSI0", str_c("BETA.", dimnames(X.PSI)[[2]]),
            "LAMBDA0", str_c("DELTA.", dimnames(X.LAMBDA)[[2]]),
            "psi0", str_c("beta.", dimnames(X.psi)[[2]]),
            "lambda0", str_c("delta.", dimnames(X.lambda)[[2]]),
            "pStar", str_c("zeta.", dimnames(X.zeta)[[2]]))
cols <- (expand.grid(out.vals, params, stringsAsFactors = F) %>%
  select(Var2, Var1) %>%
  mutate(Var3 = str_c(Var2, Var1, sep = ".")))$Var3
cols <- cols[-which(str_detect(cols, ".f") & (str_detect(cols, "PSI") |
                                                str_detect(cols, "LAMBDA") |
                                                str_detect(cols, "psi") |
                                                str_detect(cols, "lambda") |
                                                str_detect(cols, "pStar")))]
out <- matrix(NA, nrow = length(spp.list), ncol = length(cols),
              dimnames = list(spp.list, cols))

parm <- QSLpersonal::expit(mod$mcmcOutput$BETA0)
out[, "PSI0.est"] <- str_c(
  apply(parm, 2, median) %>% round(digits = 2),
  "(",
  apply(parm, 2, function(x) quantile(x, prob = 0.1, type = 8)) %>% round(digits = 2),
  ",",
  apply(parm, 2, function(x) quantile(x, prob = 0.9, type = 8)) %>% round(digits = 2),
  ")")

parm <- mod$mcmcOutput$DELTA0
out[, "LAMBDA0.est"] <- str_c(
  apply(parm, 2, median) %>% round(digits = 2),
  "(",
  apply(parm, 2, function(x) quantile(x, prob = 0.1, type = 8)) %>% round(digits = 2),
  ",",
  apply(parm, 2, function(x) quantile(x, prob = 0.9, type = 8)) %>% round(digits = 2),
  ")")

parm <- QSLpersonal::expit(mod$mcmcOutput$beta0)
out[, "psi0.est"] <- str_c(
  apply(parm, 2, median) %>% round(digits = 2),
  "(",
  apply(parm, 2, function(x) quantile(x, prob = 0.1, type = 8)) %>% round(digits = 2),
  ",",
  apply(parm, 2, function(x) quantile(x, prob = 0.9, type = 8)) %>% round(digits = 2),
  ")")

parm <- mod$mcmcOutput$delta0
out[, "lambda0.est"] <- str_c(
  apply(parm, 2, median) %>% round(digits = 2),
  "(",
  apply(parm, 2, function(x) quantile(x, prob = 0.1, type = 8)) %>% round(digits = 2),
  ",",
  apply(parm, 2, function(x) quantile(x, prob = 0.9, type = 8)) %>% round(digits = 2),
  ")")

parm <- QSLpersonal::expit(mod$mcmcOutput$zeta0) %>% (function(x) 1 - (1-x)^3)
out[, "pStar.est"] <- str_c(
  apply(parm, 2, median) %>% round(digits = 2),
  "(",
  apply(parm, 2, function(x) quantile(x, prob = 0.1, type = 8)) %>% round(digits = 2),
  ",",
  apply(parm, 2, function(x) quantile(x, prob = 0.9, type = 8)) %>% round(digits = 2),
  ")")

pars <- params[which(str_detect(params, "BETA."))]
for(i in 1:length(pars)) {
  parm <- mod$mcmcOutput$BETA1[,,i]
  out[, str_c(pars[i], ".est")] <- str_c(
    apply(parm, 2, median) %>% round(digits = 2),
    "(",
    apply(parm, 2, function(x) quantile(x, prob = 0.1, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.9, type = 8)) %>% round(digits = 2),
    ")")
  out[, str_c(pars[i], ".f")] <-
    apply(parm, 2, function(x) sum(x > 0) / length(x)) %>%
    round(digits = 2)
}

pars <- params[which(str_detect(params, "DELTA."))]
for(i in 1:length(pars)) {
  parm <- mod$mcmcOutput$DELTA1[,,i]
  out[, str_c(pars[i], ".est")] <- str_c(
    apply(parm, 2, median) %>% round(digits = 2),
    "(",
    apply(parm, 2, function(x) quantile(x, prob = 0.1, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.9, type = 8)) %>% round(digits = 2),
    ")")
  out[, str_c(pars[i], ".f")] <-
    apply(parm, 2, function(x) sum(x > 0) / length(x)) %>%
    round(digits = 2)
}

pars <- params[which(str_detect(params, "beta."))]
for(i in 1:length(pars)) {
  parm <- mod$mcmcOutput$beta1[,,i]
  out[, str_c(pars[i], ".est")] <- str_c(
    apply(parm, 2, median) %>% round(digits = 2),
    "(",
    apply(parm, 2, function(x) quantile(x, prob = 0.1, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.9, type = 8)) %>% round(digits = 2),
    ")")
  out[, str_c(pars[i], ".f")] <-
    apply(parm, 2, function(x) sum(x > 0) / length(x)) %>%
    round(digits = 2)
}

pars <- params[which(str_detect(params, "delta."))]
for(i in 1:length(pars)) {
  parm <- mod$mcmcOutput$delta1[,,i]
  out[, str_c(pars[i], ".est")] <- str_c(
    apply(parm, 2, median) %>% round(digits = 2),
    "(",
    apply(parm, 2, function(x) quantile(x, prob = 0.1, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.9, type = 8)) %>% round(digits = 2),
    ")")
  out[, str_c(pars[i], ".f")] <-
    apply(parm, 2, function(x) sum(x > 0) / length(x)) %>%
    round(digits = 2)
}

pars <- params[which(str_detect(params, "zeta."))]
for(i in 1:length(pars)) {
  parm <- mod$mcmcOutput$zeta1[,,i]
  out[, str_c(pars[i], ".est")] <- str_c(
    apply(parm, 2, median) %>% round(digits = 2),
    "(",
    apply(parm, 2, function(x) quantile(x, prob = 0.1, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.9, type = 8)) %>% round(digits = 2),
    ")")
  out[, str_c(pars[i], ".f")] <-
    apply(parm, 2, function(x) sum(x > 0) / length(x)) %>%
    round(digits = 2)
}

write.csv(out, "Parameter_est.csv", row.names = T)

## Detectability table for appendix ##
out.det <- out %>% data.frame(stringsAsFactors = F) %>%
  mutate(Spp = row.names(out)) %>%
  mutate(zeta.ShrubCov.est = ifelse(zeta.ShrubCov.f < 0.1 | zeta.ShrubCov.f > 0.9,
                                    str_c(zeta.ShrubCov.est, "*"), zeta.ShrubCov.est),
         zeta.DOY.est = ifelse(zeta.DOY.f < 0.1 | zeta.DOY.f > 0.9,
                                    str_c(zeta.DOY.est, "*"), zeta.DOY.est),
         zeta.Time_ssr.est = ifelse(zeta.Time_ssr.f < 0.1 | zeta.Time_ssr.f > 0.9,
                                    str_c(zeta.Time_ssr.est, "*"), zeta.Time_ssr.est),
         zeta.Dev_lo.est = ifelse(zeta.Dev_lo.f < 0.1 | zeta.Dev_lo.f > 0.9,
                                    str_c(zeta.Dev_lo.est, "*"), zeta.Dev_lo.est),
         zeta.Dev_bg.est = ifelse(zeta.Dev_bg.f < 0.1 | zeta.Dev_bg.f > 0.9,
                                    str_c(zeta.Dev_bg.est, "*"), zeta.Dev_bg.est),
         zeta.DOY2.est = ifelse(zeta.DOY2.f < 0.1 | zeta.DOY2.f > 0.9,
                                    str_c(zeta.DOY2.est, "*"), zeta.DOY2.est)) %>%
  select(Spp, pStar.est, zeta.ShrubCov.est, zeta.DOY.est,
         zeta.DOY2.est, zeta.Time_ssr.est, zeta.Dev_lo.est,
         zeta.Dev_bg.est)
write.csv(out.det, "Detectability_MS.csv", row.names = F)
out.det %>%
  filter(str_sub(zeta.ShrubCov.est, -1, -1) == "*" |
           str_sub(zeta.DOY.est, -1, -1) == "*" |
           str_sub(zeta.DOY2.est, -1, -1) == "*" |
           str_sub(zeta.Time_ssr.est, -1, -1) == "*" |
           str_sub(zeta.Dev_lo.est, -1, -1) == "*" |
           str_sub(zeta.Dev_bg.est, -1, -1) == "*") %>%
  nrow

## Hyper-parameter table for appendix ##
hpars <- c("omega", "rho.bB", "rho.zb",
           
           "BETA0.mu", "sigma.BETA0", "sigma.B0",
           "BETA1.mu[1]", "sigma.BETA1[1]",
           "BETA1.mu[2]", "sigma.BETA1[2]",
           "BETA1.mu[3]", "sigma.BETA1[3]",
           "BETA1.mu[4]", "sigma.BETA1[4]",
           "BETA1.mu[5]", "sigma.BETA1[5]",
           "BETA1.mu[6]", "sigma.BETA1[6]",
           
           "DELTA0.mu", "sigma.DELTA0",
           "DELTA1.mu[1]", "sigma.DELTA1[1]",
           "DELTA1.mu[2]", "sigma.DELTA1[2]",
           "DELTA1.mu[3]", "sigma.DELTA1[3]",
           
           "beta0.mu", "sigma.beta0", "sigma.b0",
           "beta1.mu[1]", "sigma.beta1[1]",
           "beta1.mu[2]", "sigma.beta1[2]",
           "beta1.mu[3]", "sigma.beta1[3]",
           "beta1.mu[4]", "sigma.beta1[4]",
           "beta1.mu[5]", "sigma.beta1[5]",
           "beta1.mu[6]", "sigma.beta1[6]",
           "beta1.mu[7]", "sigma.beta1[7]",
           "beta1.mu[8]", "sigma.beta1[8]",
           
           "delta0.mu", "sigma.delta0",
           "delta1.mu[1]", "sigma.delta1[1]",
           "delta1.mu[2]", "sigma.delta1[2]",
           "delta1.mu[3]", "sigma.delta1[3]",
           "delta1.mu[4]", "sigma.delta1[4]",
           "delta1.mu[5]", "sigma.delta1[5]",
           
           "zeta0.mu", "sigma.zeta0", "sigma.z0",
           "zeta1.mu[1]", "sigma.zeta1[1]",
           "zeta1.mu[2]", "sigma.zeta1[2]",
           "zeta1.mu[3]", "sigma.zeta1[3]",
           "zeta1.mu[4]", "sigma.zeta1[4]",
           "zeta1.mu[5]", "sigma.zeta1[5]",
           "zeta1.mu[6]", "sigma.zeta1[6]",
           
           "ALPHA0.Well_3km", "ALPHA.Dev_lo.Well_3km",
           "ALPHA.Dev_bg.Well_3km", "r.Well_3km",
           
           "ALPHA0.Well_1km", "ALPHA.Dev_lo.Well_1km",
           "ALPHA.Dev_bg.Well_1km", "r.Well_1km",
           
           "alpha0.Road_125m", "alpha.Dev_lo.Road_125m",
           "alpha.Dev_bg.Road_125m", "shape.Road_125m",
           
           "alpha0.AHerb", "alpha.Dev_lo.AHerb",
           "alpha.Dev_bg.AHerb", "alpha.Well_1km.AHerb",
           "alpha.Road_125m.AHerb", "phi.AHerb")

sum.vals <- mod$mcmcOutput[,hpars] %>%
  apply(2, function(x) QSLpersonal::BCI(x, BCIpercent = 80, flag.sig = F))
tab.vals <- cbind(Hpar = names(sum.vals), Est = sum.vals)
write.csv(tab.vals, "HyperPars_MS.csv", row.names = F)
