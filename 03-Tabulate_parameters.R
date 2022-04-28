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
  filter(!str_detect(Parameter, "BETA0\\[") &
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
           !str_detect(Parameter, "w\\[")) %>%
  write.csv(str_c("HyperParams.csv"), row.names = F)

# Compile data to get covariate names #
reduce.data.aug <- F # If TRUE, reduce data augmentation to 10 additional species.
development <- F # Set to TRUE for running test model with only develop.spp, and FALSE to run the full model.
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
  apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
  ",",
  apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
  ")")

parm <- mod$mcmcOutput$DELTA0
out[, "LAMBDA0.est"] <- str_c(
  apply(parm, 2, median) %>% round(digits = 2),
  "(",
  apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
  ",",
  apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
  ")")

parm <- QSLpersonal::expit(mod$mcmcOutput$beta0)
out[, "psi0.est"] <- str_c(
  apply(parm, 2, median) %>% round(digits = 2),
  "(",
  apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
  ",",
  apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
  ")")

parm <- mod$mcmcOutput$delta0
out[, "lambda0.est"] <- str_c(
  apply(parm, 2, median) %>% round(digits = 2),
  "(",
  apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
  ",",
  apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
  ")")

parm <- QSLpersonal::expit(mod$mcmcOutput$zeta0) %>% (function(x) 1 - (1-x)^3)
out[, "pStar.est"] <- str_c(
  apply(parm, 2, median) %>% round(digits = 2),
  "(",
  apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
  ",",
  apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
  ")")

pars <- params[which(str_detect(params, "BETA."))]
for(i in 1:length(pars)) {
  parm <- mod$mcmcOutput$BETA1[,,i]
  out[, str_c(pars[i], ".est")] <- str_c(
    apply(parm, 2, median) %>% round(digits = 2),
    "(",
    apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
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
    apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
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
    apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
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
    apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
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
    apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
    ")")
  out[, str_c(pars[i], ".f")] <-
    apply(parm, 2, function(x) sum(x > 0) / length(x)) %>%
    round(digits = 2)
}

write.csv(out, "Parameter_est.csv", row.names = T)

# # Long format for supplementary materials (Not sure yet if I'll do this.) #
# 
# beta_var_nams <- dimnames(X.beta)[[3]]
# beta_var_nams[3] <- "PAROpn"
# alpha_var_nams <- dimnames(X.alpha)[[2]]
# zeta_var_nams <- dimnames(X.zeta)[[2]]
# 
# tab_sum <- mod$summary %>% select(Parameter, median:u95) %>%
#   mutate(Species = ifelse(Parameter == "omega" |
#                             str_sub(Parameter, 1, 3) == "rho" |
#                             str_detect(Parameter, "mu") |
#                             str_detect(Parameter, "sigma"), "Community", ""))
# tab_sum$Species[which(tab_sum$Species == "")] <- spp.list
# 
# # Community parameters #
# tab_sum <- tab_sum %>% filter(Parameter == "omega") %>%
#   bind_rows(tab_sum %>% filter(Parameter == "rho.ab")) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "rho.za")) %>%
#   #psi Intercepts
#   bind_rows(tab_sum %>% filter(Parameter == "beta0.mu")) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.beta0")) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "beta0[") %>%
#               mutate(Parameter = "beta0")) %>%
#   #Psi covariate relationships
#     #Canopy gap extent
#   bind_rows(tab_sum %>% filter(Parameter == "beta1.mu[1]") %>%
#               mutate(Parameter = str_c("beta.", beta_var_nams[1], ".mu"))) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.beta1[1]") %>%
#               mutate(Parameter = str_c("sigma.beta.", beta_var_nams[1]))) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "beta1[" & str_sub(Parameter, -2, -1) == "1]") %>%
#               mutate(Parameter = str_c("beta.", beta_var_nams[1]))) %>%
#     #Open forest extent
#   bind_rows(tab_sum %>% filter(Parameter == "beta1.mu[2]") %>%
#               mutate(Parameter = str_c("beta.", beta_var_nams[2], ".mu"))) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.beta1[2]") %>%
#               mutate(Parameter = str_c("sigma.beta.", beta_var_nams[2]))) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "beta1[" & str_sub(Parameter, -2, -1) == "2]") %>%
#               mutate(Parameter = str_c("beta.", beta_var_nams[2]))) %>%
#     #Perimeter-area ratio for open forest
#   bind_rows(tab_sum %>% filter(Parameter == "beta1.mu[3]") %>%
#               mutate(Parameter = str_c("beta.", beta_var_nams[3], ".mu"))) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.beta1[3]") %>%
#               mutate(Parameter = str_c("sigma.beta.", beta_var_nams[3]))) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "beta1[" & str_sub(Parameter, -2, -1) == "3]") %>%
#               mutate(Parameter = str_c("beta.", beta_var_nams[3]))) %>%
#     #Latitude
#   bind_rows(tab_sum %>% filter(Parameter == "beta1.mu[4]") %>%
#               mutate(Parameter = str_c("beta.", beta_var_nams[4], ".mu"))) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.beta1[4]") %>%
#               mutate(Parameter = str_c("sigma.beta.", beta_var_nams[4]))) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "beta1[" & str_sub(Parameter, -2, -1) == "4]") %>%
#               mutate(Parameter = str_c("beta.", beta_var_nams[4]))) %>%
#     #Heat load
#   bind_rows(tab_sum %>% filter(Parameter == "beta1.mu[5]") %>%
#               mutate(Parameter = str_c("beta.", beta_var_nams[5], ".mu"))) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.beta1[5]") %>%
#               mutate(Parameter = str_c("sigma.beta.", beta_var_nams[5]))) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "beta1[" & str_sub(Parameter, -2, -1) == "5]") %>%
#               mutate(Parameter = str_c("beta.", beta_var_nams[5]))) %>%
#     #TWI
#   bind_rows(tab_sum %>% filter(Parameter == "beta1.mu[6]") %>%
#               mutate(Parameter = str_c("beta.", beta_var_nams[6], ".mu"))) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.beta1[6]") %>%
#               mutate(Parameter = str_c("sigma.beta.", beta_var_nams[6]))) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "beta1[" & str_sub(Parameter, -2, -1) == "6]") %>%
#               mutate(Parameter = str_c("beta.", beta_var_nams[6]))) %>%
#     #Life zone (1 = lower montane)
#   bind_rows(tab_sum %>% filter(Parameter == "beta1.mu[7]") %>%
#               mutate(Parameter = str_c("beta.", beta_var_nams[7], ".mu"))) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.beta1[7]") %>%
#               mutate(Parameter = str_c("sigma.beta.", beta_var_nams[7]))) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "beta1[" & str_sub(Parameter, -2, -1) == "7]") %>%
#               mutate(Parameter = str_c("beta.", beta_var_nams[7]))) %>%
#   #Year-specific deviations from mean logit(psi)
#     #Year 1
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "dev.b0" & str_sub(Parameter, -2, -1) == "1]") %>%
#               mutate(Parameter = "dev.b0.t1")) %>%
#     #Year 2
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "dev.b0" & str_sub(Parameter, -2, -1) == "2]") %>%
#               mutate(Parameter = "dev.b0.t2")) %>%
#     #Year 3
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "dev.b0" & str_sub(Parameter, -2, -1) == "3]") %>%
#               mutate(Parameter = "dev.b0.t3")) %>%
#   #theta intercepts
#   bind_rows(tab_sum %>% filter(Parameter == "alpha0.mu")) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.alpha0")) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "alpha0[") %>%
#               mutate(Parameter = "alpha0")) %>%
#   #theta covariate relationships
#     #Canopy cover
#   bind_rows(tab_sum %>% filter(Parameter == "alpha1.mu[1]") %>%
#               mutate(Parameter = str_c("alpha.", alpha_var_nams[1], ".mu"))) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.alpha1[1]") %>%
#               mutate(Parameter = str_c("sigma.alpha.", alpha_var_nams[1]))) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "alpha1[" & str_sub(Parameter, -2, -1) == "1]") %>%
#               mutate(Parameter = str_c("alpha.", alpha_var_nams[1]))) %>%
#     #Canopy cover squared
#   bind_rows(tab_sum %>% filter(Parameter == "alpha1.mu[2]") %>%
#               mutate(Parameter = str_c("alpha.", alpha_var_nams[2], ".mu"))) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.alpha1[2]") %>%
#               mutate(Parameter = str_c("sigma.alpha.", alpha_var_nams[2]))) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "alpha1[" & str_sub(Parameter, -2, -1) == "2]") %>%
#               mutate(Parameter = str_c("alpha.", alpha_var_nams[2]))) %>%
#   #p intercepts
#   bind_rows(tab_sum %>% filter(Parameter == "zeta0.mu"))  %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.zeta0"))  %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "zeta0[") %>%
#               mutate(Parameter = "zeta0")) %>%
#   #p covariate relationships
#     #Canopy cover
#   bind_rows(tab_sum %>% filter(Parameter == "zeta1.mu[1]") %>%
#               mutate(Parameter = str_c("zeta.", zeta_var_nams[1], ".mu"))) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.zeta1[1]") %>%
#               mutate(Parameter = str_c("sigma.zeta.", zeta_var_nams[1]))) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "zeta1[" & str_sub(Parameter, -2, -1) == "1]") %>%
#               mutate(Parameter = str_c("zeta.", zeta_var_nams[1]))) %>%
#     #Day of year
#   bind_rows(tab_sum %>% filter(Parameter == "zeta1.mu[2]") %>%
#               mutate(Parameter = str_c("zeta.", zeta_var_nams[2], ".mu"))) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.zeta1[2]") %>%
#               mutate(Parameter = str_c("sigma.zeta.", zeta_var_nams[2]))) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "zeta1[" & str_sub(Parameter, -2, -1) == "2]") %>%
#               mutate(Parameter = str_c("zeta.", zeta_var_nams[2]))) %>%
#     #Day of year squared
#   bind_rows(tab_sum %>% filter(Parameter == "zeta1.mu[4]") %>%
#               mutate(Parameter = str_c("zeta.", zeta_var_nams[4], ".mu"))) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.zeta1[4]") %>%
#               mutate(Parameter = str_c("sigma.zeta.", zeta_var_nams[4]))) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "zeta1[" & str_sub(Parameter, -2, -1) == "4]") %>%
#               mutate(Parameter = str_c("zeta.", zeta_var_nams[4]))) %>%
#     #Time since sunrise
#   bind_rows(tab_sum %>% filter(Parameter == "zeta1.mu[3]") %>%
#               mutate(Parameter = str_c("zeta.", zeta_var_nams[3], ".mu"))) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.zeta1[3]") %>%
#               mutate(Parameter = str_c("sigma.zeta.", zeta_var_nams[3]))) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "zeta1[" & str_sub(Parameter, -2, -1) == "3]") %>%
#               mutate(Parameter = str_c("zeta.", zeta_var_nams[3]))) %>%
#     #Time since sunrise squared
#   bind_rows(tab_sum %>% filter(Parameter == "zeta1.mu[5]") %>%
#               mutate(Parameter = str_c("zeta.", zeta_var_nams[5], ".mu"))) %>%
#   bind_rows(tab_sum %>% filter(Parameter == "sigma.zeta1[5]") %>%
#               mutate(Parameter = str_c("sigma.zeta.", zeta_var_nams[5]))) %>%
#   bind_rows(tab_sum %>% filter(str_sub(Parameter, 1, 6) == "zeta1[" & str_sub(Parameter, -2, -1) == "5]") %>%
#               mutate(Parameter = str_c("zeta.", zeta_var_nams[5]))) %>%
#   select(Parameter, Species, median:u95)
# 
# write.csv(tab_sum, "Parameter_est_AppndS3.csv", row.names = F)