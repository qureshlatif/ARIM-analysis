library(nimble)
library(stringr)
library(tidyr)
library(dplyr)

setwd("C:/Users/Quresh.Latif/files/projects/BLM/ARIM")
load("Data_compiled.RData")

#_____ Script inputs _____#
scripts.loc <- "ARIM-analysis/"
model.file <- str_c(scripts.loc, "model_path.nimble")
mod.nam <- "mod_path"
reduce.data.aug <- F # If TRUE, reduce data augmentation to 10 additional species.
development <- F # Set to TRUE for running test model with only develop.spp, and FALSE to run the full model.
develop.spp <- c("MODO", "MOPL", "HOLA", "WEME")

# MCMC values
nc <- 2 # number of chains
nb <- 10000 # burn in
ni <- 50000 # number of iterations
nt <- 10 # thinning
#_________________________#

# Compile data #
source(str_c(scripts.loc, "Data_processing.R"))

# Data objects to send to NIMBLE
# Parameters to set up variance-covariance prior for PSI, psi, and p
#R <- matrix(c(5, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
#Sig_df <- 4

data.nams <- c("Y.sum", "Y.int", "X.PSI.raw", "X.psi.raw")

constant.nams <- c("n.grdyr", "yearID.grdyr", #"gridID", "yearID", "pointID", "Y.point", "Y.grid", "Y.spp", "TPeriod", 
                   "n.year", "n.point", "K", "n.spp", #"n.grid", "n.pntyr", "gridXyrID",  "R", "Sig_df",
                   "n.int", "grdyrID.int", "sppID.int",
                   
                   "X.PSI", "p.PSI",
                   # "ind.PSI.offset", "ind.PSI.no_offset",
                   
                   "X.psi", "X.psi.dyn", "p.psi.init", "p.psi.dyn", "X.psi.dyn",
                   #"ind.psi.init.offset", "ind.psi.init.no_offset",
                   #"ind.psi.dyn.offset", "ind.psi.dyn.no_offset",
                   
                   "X.zeta", "p.zeta")#,

#"guildMem", "n.guild")
if(mod.nam %in% c("mod_path", "mod_interm_paths")) {
  constant.nams <- c(constant.nams,
                     
                     "ind.PSI.Dev_bg", "ind.PSI.Dev_lo", "ind.PSI.Well_3km",
                     "ind.PSI.Well_1km", "ind.PSI.Road_1km",
                     
                     "ind.psi.Road_125m", "ind.psi.Well_1km", "ind.psi.AHerb") # "ind.D_Road_125m", "ind.PA_Road_125m", "ind.Well_125m"
}

# Stuff to save from NIMBLE
parameters <- c(# Bird community parameters
  "omega", "rho.zb", "rho.bB",# "rho.zB",
  
  "BETA0.mu", "sigma.BETA0", "sigma.B0", "BETA1.mu", "sigma.BETA1",
  "DELTA0.mu", "sigma.DELTA0", "DELTA1.mu", "sigma.DELTA1",
  "beta0.mu", "sigma.beta0", "sigma.b0", "beta1.mu", "sigma.beta1",
  "delta0.mu", "sigma.delta0", "delta1.mu", "sigma.delta1",
  "zeta0.mu", "sigma.zeta0", "sigma.z0", "zeta1.mu", "sigma.zeta1",
  
  "BETA0", "dev.BETA", "BETA1",# "BETA1.offset",
  "DELTA0", "DELTA1",# "DELTA1.offset",
  "beta0", "dev.beta", "beta1",# "beta1.offset",
  "delta0", "delta1",# "delta1.offset",
  "zeta0", "dev.zeta", "zeta1",
  
  # Interm path parameters #
  "ALPHA0.Well_3km", "ALPHA.Dev_lo.Well_3km",
  "ALPHA.Dev_bg.Well_3km", "r.Well_3km",
  
  "ALPHA0.Well_1km", "ALPHA.Dev_lo.Well_1km",
  "ALPHA.Dev_bg.Well_1km", "r.Well_1km",
  
  "ALPHA0.Road_1km", "ALPHA.Dev_lo.Road_1km",
  "ALPHA.Dev_bg.Road_1km", "shape.Road_1km",
  
  "alpha0.Road_125m",
  "alpha.Dev_lo.Road_125m", "alpha.Dev_bg.Road_125m",
  "shape.Road_125m",
  
  "alpha0.AHerb", "alpha.Well_1km.AHerb",
  "alpha.Road_1km.AHerb")

# Function for setting initial values in JAGS
inits <- function()
  list(z=z.init, Z=Z.init, w=w.init,
       tvar.sigma.zeta0 = rnorm(1),
       tvar.sigma.z0 = rnorm(1),
       tvar.zeta1 = rnorm(p.zeta),
       
       tvar.sigma.BETA0 = rnorm(1),
       tvar.sigma.B0 = rnorm(1),
       tvar.BETA1 = rnorm(p.PSI),
       
       tvar.sigma.DELTA0 = rnorm(1),
       tvar.DELTA1 = rnorm(p.PSI),
       
       tvar.sigma.beta0 = rnorm(1),
       tvar.sigma.b0 = rnorm(1),
       tvar.beta1 = rnorm(p.psi.init),
       
       tvar.sigma.delta0 = rnorm(1),
       tvar.delta1 = rnorm(p.psi.dyn),
       
       # Needed for full initialization #
       BETA0 = BETA0.init, #rnorm(n.spp),
       BETA1 = matrix(0, nrow = n.spp, ncol = p.PSI), #rnorm(n.spp * p.PSI)
       #BETA1.base = matrix(0, nrow = n.spp, ncol = p.PSI), #rnorm(n.spp * p.PSI)
       #BETA1.offset = matrix(0, nrow = (n.guild-1), ncol = p.PSI),
       DELTA0 = rep(0, n.spp), #rnorm(n.spp),
       DELTA1 = matrix(0, nrow = n.spp, ncol = p.PSI), #rnorm(n.spp * p.PSI)
       #DELTA1.base = matrix(0, nrow = n.spp, ncol = p.PSI), #rnorm(n.spp * p.PSI)
       #DELTA1.offset = matrix(0, nrow = (n.guild-1), ncol = p.PSI),
       beta0 = beta0.init, #rnorm(n.spp),
       beta1 = matrix(0, nrow = n.spp, ncol = p.psi.init), #rnorm(n.spp * p.psi.init)
       #beta1.base = matrix(0, nrow = n.spp, ncol = p.psi.init), #rnorm(n.spp * p.psi.init)
       #beta1.offset = matrix(0, nrow = (n.guild-1), ncol = p.psi.init),
       delta0 = rep(0, n.spp), #rnorm(n.spp),
       delta1 = matrix(0, nrow = n.spp, ncol = p.psi.dyn), # rnorm(n.spp * p.psi.dyn)
       #delta1.base = matrix(0, nrow = n.spp, ncol = p.psi.dyn), # rnorm(n.spp * p.psi.dyn)
       #delta1.offset = matrix(0, nrow = (n.guild-1), ncol = p.psi.dyn),
       zeta0 = rep(2, n.spp),
       zeta1 = matrix(0, nrow = n.spp, ncol = p.zeta),
       
       dev.zeta = matrix(0, nrow = n.spp, ncol = n.year),
       dev.beta = matrix(0, nrow = n.spp, ncol = n.year),
       dev.BETA = matrix(0, nrow = n.spp, ncol = n.year),
       
       rho.zb = runif(1, -1, 1),
       rho.bB = runif(1, -1, 1),
       
       zeta0.mu = 2,
       zeta1.mu = rep(0, p.zeta), #rnorm(p.zeta),
       beta0.mu = mean(beta0.init),
       beta1.mu = rep(0, p.psi.init), #rnorm(p.psi.init),
       delta0.mu = 0, #rnorm(1),
       delta1.mu = rep(0, p.psi.dyn), #rnorm(p.psi.dyn),
       BETA0.mu = mean(BETA0.init),
       BETA1.mu = rep(0, p.PSI), #rnorm(p.PSI),
       DELTA0.mu = 0, #rnorm(1),
       DELTA1.mu = rep(0, p.PSI), #rnorm(p.PSI),
       
       omega = runif(1, 0.8201439, 0.99),
       
       ALPHA0.Well_3km = rnorm(1),
       ALPHA.Dev_lo.Well_3km = rnorm(1),
       ALPHA.Dev_bg.Well_3km = rnorm(1),
       r.Well_3km = runif(1),
       
       ALPHA0.Well_1km = rnorm(1),
       ALPHA.Dev_lo.Well_1km = rnorm(1),
       ALPHA.Dev_bg.Well_1km = rnorm(1),
       r.Well_1km = runif(1),
       
       ALPHA0.Road_1km = rnorm(1),
       ALPHA.Dev_lo.Road_1km = rnorm(1),
       ALPHA.Dev_bg.Road_1km = rnorm(1),
       shape.Road_1km = runif(1),
       
       alpha0.AHerb = rnorm(1),
       alpha.Well_1km.AHerb = rnorm(1),
       alpha.Road_1km.AHerb = rnorm(1),
       phi.AHerb = rgamma(1, 1),
       
       alpha0.Road_125m = rnorm(1),
       alpha.Dev_lo.Road_125m = rnorm(1),
       alpha.Dev_bg.Road_125m = rnorm(1),
       shape.Road_125m = runif(1))

# Assemble the initial values for JAGS.
z.init <- Y.sum #array(n.point, dim = c(n.grdyr, n.spp)) #array(1, dim = c(n.point, n.year, n.spp))
Z.init <- (Y.sum > 0) * 1 #array(1, dim = c(n.grid, n.year, n.spp))
w.init <- apply(Y.sum, 2, function(x) any(x > 0) * 1)

# More thoughtful initial values #
logit <- function(x) log(x/(1-x))
BETA0.init <- apply(Y.sum, 2, function(x) sum(x > 0) / length(x))
BETA0.init <- ifelse(BETA0.init == 0, logit(0.0001), logit(BETA0.init))
beta0.init <- apply(Y.sum, 2, function(x) x / n.point)
beta0.init[which(beta0.init == 0)] <- NA
beta0.init <- apply(beta0.init, 2, function(x) ifelse(any(!is.na(x)), mean(x, na.rm = T), 0))
beta0.init <- ifelse(beta0.init == 0, logit(0.0001), logit(beta0.init))

data <- list()
for(i in 1:length(data.nams)) data[[length(data) + 1]] <- eval(as.name(data.nams[i]))
names(data) <- data.nams

constants <- list()
for(i in 1:length(constant.nams)) constants[[length(constants) + 1]] <- eval(as.name(constant.nams[i]))
names(constants) <- constant.nams

# Fit model
source(model.file)
rm(.Random.seed, envir=.GlobalEnv)

# All at once using built-in wrapper... #
# st.time <- Sys.time()
# out <- nimbleMCMC(code = model,
#                   constants = constants,
#                   data=data,
#                   inits=inits, 
#                   nchains = nc,
#                   nburnin = nb,
#                   niter = ni,
#                   thin = nt,
#                   samplesAsCodaMCMC = T,
#                   summary = FALSE,
#                   WAIC = FALSE,
#                   monitors = parameters)
# end.time <- Sys.time()
# run.time <- end.time - st.time
# run.time
# rm(st.time,end.time)
# R.utils::saveObject(out, str_c(mod.nam, "_samples")) # If running chains in parallel.


#...or step through for customization #
# Step 1
st.time <- Sys.time()
mod <- nimbleModel(code = model,
                   constants = constants,
                   data = data,
                   inits = inits(),
                   returnDef = FALSE,
                   check = FALSE)
end.time <- Sys.time()
end.time - st.time

# Step 2
st.time <- Sys.time()
mcmc <- configureMCMC(mod,
                      monitors = parameters,
                      enableWAIC = F,
                      control = list(adaptInterval = 50))
#mcmc$removeSamplers(c("...")) # Switch out samplers here as desired.
#mcmc$addSampler(c("..."), type = "...")
end.time <- Sys.time()
end.time - st.time

# Step 3
st.time <- Sys.time()
mcmc <- buildMCMC(mcmc)
end.time <- Sys.time()
end.time - st.time

# Step 4
st.time <- Sys.time()
Cmodel <- compileNimble(mod)
end.time <- Sys.time()
end.time - st.time

# Step 5
st.time <- Sys.time()
Cmcmc <- compileNimble(mcmc, project = Cmodel)
end.time <- Sys.time()
end.time - st.time

# Step 6
st.time <- Sys.time()
out <- runMCMC(Cmcmc, nchains = nc, nburnin = nb, niter = ni, thin = nt)
#Cmcmc$run(ni) # Not understanding how this works.
#out <- Cmcmc$mvSamples[...]
end.time <- Sys.time()
end.time - st.time

# Step 7 (Continue if not converged; Not clear to me how this works.)
# Cmcmc$mvSamples # Does this have samples in it??
#Cmcmc$run(50000, reset = FALSE) # ??

# Save chain
R.utils::saveObject(out, str_c(mod.nam, "_chain", chain)) # If running chains in parallel.
rm(out)

# Recover independent chains #
if(nc == 1) {
  out_chain1 <- R.utils::loadObject(str_c(mod.nam, "_chain1"))
  out_chain2 <- R.utils::loadObject(str_c(mod.nam, "_chain2"))
  out.samples <- list(chain1 = out_chain1, chain2 = out_chain2)
  #out$samples$chain1 <- out$samples$chain1[-c(1:nb),][seq(1,(ni-nb),by=10),] # If burnin and thinning not implemented with model run.
  #out$samples$chain2 <- out$samples$chain2[-c(1:nb),][seq(1,(ni-nb),by=10),] # If burnin and thinning not implemented with model run.
  mod.raw <- coda::as.mcmc.list(lapply(out.samples, coda::mcmc))
}

library(mcmcOutput)
#if(nc > 1) mod.raw <- coda::as.mcmc.list(lapply(out$samples, coda::mcmc))
if(nc > 1) {
  mod <- mcmcOutput(out)
} else {
  mod <- mcmcOutput(mod.raw)
}
sumTab <- summary(mod, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = F)
sumTab <- sumTab %>%
  as_tibble() %>%
  mutate(Parameter = row.names(sumTab)) %>%
  select(Parameter, mean:f)

library(MCMCvis)
MCMCvis::MCMCtrace(mod.raw, params = "zeta1[100, 1]", pdf = F, ISB = F)
mod <- list(mcmcOutput = mod, summary = sumTab)

R.utils::saveObject(mod, mod.nam) # If running all in one.

#   R.utils::saveObject(out, mod.nam)

#################################################################
# Parallelization and updating code from NIMBLE user group post #
#################################################################

# library(parallel)
# library(coda)
# #set.seed(22)
# cl<-makeCluster(nc,timeout=5184000)
# inits <- function() {list(a = rnorm(1), b = rnorm(1))}
# nimbledata <- list(y = data)
# nimbleconstants <- list(length_y = 1000)
# params <- c("a", "b")
# clusterExport(cl, c("myCode", "inits", "nimbledata", "nimbleconstants", "params"))
# for (j in seq_along(cl)) {
#   set.seed(j)
#   init <- inits()
#   clusterExport(cl[j], "init")
# }
# out <- clusterEvalQ(cl, {
#   library(nimble)
#   library(coda)
#   model <- nimbleModel(code = myCode, name = "myCode",
#                        constants = nimbleconstants, data = nimbledata,
#                        inits = init)
#   Cmodel <- compileNimble(model)
#   modelConf <- configureMCMC(model)
#   modelConf$addMonitors(params)
#   modelMCMC <- buildMCMC(modelConf)
#   CmodelMCMC <- compileNimble(modelMCMC, project = myCode)
#   out1 <- runMCMC(CmodelMCMC, niter = 10000)
#   return(as.mcmc(out1))
# })
# 
# out.mcmc <- as.mcmc(out)
# traceplot(out.mcmc[, "b"])
#           
# ## If has not converged, continue sampling
# start <- Sys.time()
# out2 <- clusterEvalQ(cl, {
#   out1 <- runMCMC(CmodelMCMC, niter = 20000)
#   return(as.mcmc(out1))
# })
# 
# out.mcmc.update1 <- as.mcmc(out2)
# 
# out.mcmc.bind <- mcmc.list()
# for (i in seq_len(nc)) {
#   out.mcmc.bind[[i]] <- mcmc(rbind(out.mcmc[[i]], out.mcmc.update1[[i]]))
# }
# traceplot(out.mcmc.bind[, "b"])