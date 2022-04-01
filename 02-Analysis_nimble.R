library(nimble)
library(stringr)
library(tidyr)
library(dplyr)

setwd("C:/Users/Quresh.Latif/files/projects/BLM/ARIM")
load("Data_compiled.RData")

#_____ Script inputs _____#
scripts.loc <- "ARIM-analysis/"
model.file <- str_c(scripts.loc, "model_path_trend_flatten_to_grid.nimble")
mod.nam <- "mod_path_flatten_to_grid"
reduce.data.aug <- F # If TRUE, reduce data augmentation to 10 additional species.
development <- F # Set to TRUE for running test model with only develop.spp, and FALSE to run the full model.
develop.spp <- c("MODO", "MOPL", "HOLA", "WEME")

# MCMC values
nc <- 2#3 # number of chains
nb <- 10 #1000 # burn in
ni <- 20 #3000 #11000 # number of iterations
nt <- 1#10 # thinning
#_________________________#

# Compile data #
source(str_c(scripts.loc, "Data_processing_path_flatten_to_grid.R"))

# Data objects to send to NIMBLE
# Parameters to set up variance-covariance prior for PSI, psi, and p
#R <- matrix(c(5, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
#Sig_df <- 4

data.nams <- c("Y.sum", "Y.int", "X.PSI.raw", "X.psi.raw")

constant.nams <- c("n.grdyr", "yearID.grdyr", #"gridID", "yearID", "pointID", "Y.point", "Y.grid", "Y.spp", "TPeriod", 
                   "n.year", "n.point", "K", "n.spp", #"n.grid", "n.pntyr", "gridXyrID",  "R", "Sig_df",
                   "n.int", "grdyrID.int", "sppID.int",
                   
                   "X.PSI", "p.PSI", "ind.PSI.offset", "ind.PSI.no_offset",
                   
                   "X.psi", "X.psi.dyn", "p.psi.init", "p.psi.dyn", "X.psi.dyn",
                   "ind.psi.init.offset", "ind.psi.init.no_offset",
                   "ind.psi.dyn.offset", "ind.psi.dyn.no_offset",
                   
                   "X.zeta", "p.zeta",
                   
                   "guildMem", "n.guild")
if(mod.nam == "mod_path_flatten_to_grid") {
  constant.nams <- c(constant.nams,
                 
                 "ind.PSI.Dev_bg", "ind.PSI.Dev_lo", "ind.Well_3km",
                 "ind.Well_1km", "ind.Road_1km",
                 
                 "ind.Well_125m", "ind.Road_125m", "ind.AHerb") # "ind.D_Road_125m", "ind.PA_Road_125m",
}

# Stuff to save from NIMBLE
parameters <- c(# Bird community parameters
  "omega", "rho.zb", "rho.bB",# "rho.zB",
  
  "BETA0.mu", "sigma.BETA0", "sigma.B0", "BETA1.mu", "sigma.BETA1",
  "DELTA0.mu", "sigma.DELTA0", "DELTA1.mu", "sigma.DELTA1",
  "beta0.mu", "sigma.beta0", "sigma.b0", "beta1.mu", "sigma.beta1",
  "delta0.mu", "sigma.delta0", "delta1.mu", "sigma.delta1",
  "zeta0.mu", "sigma.zeta0", "sigma.z0", "zeta1.mu", "sigma.zeta1",
  
  "BETA0", "dev.BETA", "BETA1", "BETA1.offset",
  "DELTA0", "DELTA1", "DELTA1.offset",
  "beta0", "dev.beta", "beta1", "beta1.offset",
  "delta0", "delta1", "delta1.offset",
  "zeta0", "dev.zeta", "zeta1",

  # Interm path parameters #
  # Interm path parameters #
  "ALPHA0.Well_3km", "ALPHA.Dev_lo.Well_3km",
  "ALPHA.Dev_bg.Well_3km", "r.Well_3km",
  
  "ALPHA0.Well_1km", "ALPHA.Dev_lo.Well_1km",
  "ALPHA.Dev_bg.Well_1km", "r.Well_1km",
  
  "ALPHA0.Road_1km", "ALPHA.Dev_lo.Road_1km",
  "ALPHA.Dev_bg.Road_1km", "shape.Road_1km",
  
  "alpha0.Well_125m", "alpha.Dev_lo.Well_125m",
  "alpha.Dev_bg.Well_125m",
  
  "alpha0.Road_125m",
  "alpha.Dev_lo.Road_125m", "alpha.Dev_bg.Road_125m",
  "shape.Road_125m",
  
  "alpha0.AHerb", "alpha.WellA_125m.AHerb",
  "alpha.Road_125m.AHerb")

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
       tvar.delta1 = rnorm(p.psi.dyn))

# Assemble the initial values for JAGS.
z.init <- Y.sum #array(n.point, dim = c(n.grdyr, n.spp)) #array(1, dim = c(n.point, n.year, n.spp))
Z.init <- (Y.sum > 0) * 1 #array(1, dim = c(n.grid, n.year, n.spp))
w.init <- apply(Y.sum, 2, function(x) any(x > 0) * 1)

data <- list()
for(i in 1:length(data.nams)) data[[length(data) + 1]] <- eval(as.name(data.nams[i]))
names(data) <- data.nams

constants <- list()
for(i in 1:length(constant.nams)) constants[[length(constants) + 1]] <- eval(as.name(constant.nams[i]))
names(constants) <- constant.nams

# Fit model
source(model.file)

# All at once using built-in wrapper... #
st.time <- Sys.time()
out <- nimbleMCMC(code = model,
                  constants = constants,
                  data=data,
                  inits=inits, 
                  nchains = nc,
                  niter = ni,
                  summary = TRUE,
                  WAIC = FALSE,
                  monitors = parameters)
end.time <- Sys.time()
run.time <- end.time - st.time
run.time
rm(st.time,end.time)

# ...or step through for customization #
# Step 1
st.time <- Sys.time()
mod <- nimbleModel(code = model,
                   constants = constants,
                   data = data,
                   inits = inits(),
                   returnDef = FALSE,
                   check = FALSE)
end.time <- Sys.time()
time <- end.time - st.time

# Step 2
st.time <- Sys.time()
mcmc <- configureMCMC(mod,
                      monitors = parameters,
                      enableWAIC = F)
mcmc$removeSamplers(c("z"))
mcmc$addSampler(c("z"), type = "RW_block", try = 4)
end.time <- Sys.time()
time <- end.time - st.time

# Step 3
st.time <- Sys.time()
mcmc <- buildMCMC(mcmc)
end.time <- Sys.time()
time <- end.time - st.time

# Step 4
st.time <- Sys.time()
Cmodel <- compileNimble(mod)
end.time <- Sys.time()
time <- end.time - st.time

# Step 5
st.time <- Sys.time()
Cmcmc <- compileNimble(mcmc, project = Cmodel)
end.time <- Sys.time()
time <- end.time - st.time

# Step 6
st.time <- Sys.time()
out <- runMCMC(Cmcmc, nchains = nc, nburnin = nb, niter = ni, thin = nt)
#Cmcmc_default$run(1000)
end.time <- Sys.time()
time <- end.time - st.time

#   R.utils::saveObject(out, mod.nam)
