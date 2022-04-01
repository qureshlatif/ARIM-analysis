library(stringr)
library(tidyr)
library(dplyr)

setwd("C:/Users/Quresh.Latif/files/projects/BLM/ARIM")
load("Data_compiled.RData")

#_____ Script inputs _____#
scripts.loc <- "ARIM-analysis/"
model.file <- str_c(scripts.loc, "model_community_mech_only_flat.jags") # "model_path_trend_flattened_marg.jags"
saveJAGS.loc <- "saveJAGS/"
package <- "saveJAGS" # Set to jagsUI or saveJAGS
mod.nam <- "mod_community_mech_only_flat"
reduce.data.aug <- F # If TRUE, reduce data augmentation to 10 additional species.
development <- F # Set to TRUE for running test model with only develop.spp, and FALSE to run the full model.
develop.spp <- c("MODO", "MOPL", "HOLA", "WEME")

# MCMC values
nc <- 2 #3 # number of chains
nb <- 10 #1000 # burn in
ni <- 20 #11000 # number of iterations
nt <- 1 #10 # thinning
ns <- 220 # number of saveJAGS chunks
#_________________________#

if(package == "jagsUI") {
  library(jagsUI)
} else {
  library(saveJAGS)
}

# Compile data #
source(str_c(scripts.loc, "Data_processing_path_flatten_to_grid.R"))

# Data objects to send to JAGS
  # Parameters to set up variance-covariance prior for PSI, psi, and p
#R <- matrix(c(5, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
#Sig_df <- 4

data.nams <- c("Y.sum", "n.grdyr", "yearID.grdyr", #"gridID", "yearID", "pointID", "Y.point", "Y.grid", "Y.spp", "TPeriod", 
               "n.year", "n.point", "K", "n.spp", #"n.grid", "n.pntyr", "gridXyrID",  "R", "Sig_df",
               "Y.int", "n.int", "grdyrID.int", "sppID.int",
               
               "X.PSI", "X.PSI.raw", "p.PSI", "ind.PSI.offset", "ind.PSI.no_offset",
               
               "X.psi", "X.psi.raw", "ind.psi.dyn", "p.psi.init", "p.psi.dyn",
               "ind.psi.init.offset", "ind.psi.init.no_offset",
               "ind.psi.dyn.offset", "ind.psi.dyn.no_offset",
               
               "X.zeta", "p.zeta",
               
               "guildMem", "n.guild")
if(mod.nam == "mod_path_flatten_to_grid") {
  data.nams <- c(data.nams,
                 
                 "ind.PSI.Dev_bg", "ind.PSI.Dev_lo", "ind.Well_3km",
                 "ind.Well_1km", "ind.Road_1km",
                 
                 "ind.Well_125m", "ind.Road_125m", "ind.AHerb") # "ind.D_Road_125m", "ind.PA_Road_125m",
}
# "X.ETA", "p.ETA", "X.DELTA", "p.DELTA", "X.eta", "p.eta", "X.delta", "p.delta",

# Stuff to save from JAGS
parameters <- c(# Bird community parameters
                "omega", "rho.zb", "rho.bB", "rho.zB",
                
                "BETA0.mu", "sigma.BETA0", "sigma.B0", "BETA1.mu", "sigma.BETA1",
                "DELTA0.mu", "sigma.DELTA0", "DELTA1.mu", "sigma.DELTA1",
                "beta0.mu", "sigma.beta0", "sigma.b0", "beta1.mu", "sigma.beta1",
                "delta0.mu", "sigma.delta0", "delta1.mu", "sigma.delta1",
                "zeta0.mu", "sigma.zeta0", "sigma.z0", "zeta1.mu", "sigma.zeta1",
                
                "BETA0", "dev.BETA", "BETA1", "BETA1.offset",
                "DELTA0", "DELTA1", "DELTA1.offset",
                "beta0", "dev.beta", "beta1", "tpi.threshold", "beta1.offset",
                "delta0", "delta1", "delta1.offset",
                "zeta0", "dev.zeta", "zeta1",
                
                "occ.pnt.yr", "occ.grd.yr", "SR.pnt.yr",
                "SR.pnt.gld.yr", "SR.grd.yr", "SR.grd.gld.yr",
                
                # Interm path parameters #
                "ALPHA0.Well_3km", "ALPHA.Dev_lo.Well_3km",
                "ALPHA.Dev_bg.Well_3km", "r.Well_3km",

                "ALPHA0.Well_1km", "ALPHA.Dev_lo.Well_1km",
                "ALPHA.Dev_bg.Well_1km", "r.Well_1km",

                "ALPHA0.Road_1km", "ALPHA.Dev_lo.Road_1km",
                "ALPHA.Dev_bg.Road_1km", "shape.Road_1km",
                
                "alpha0.Well_125m", "alpha.Dev_lo.Well_125m",
                "alpha.Dev_bg.Well_125m",
                
                #"alpha0.PA_Road_125m", "alpha.Dev_lo.PA_Road_125m",
                #"alpha.Dev_bg.PA_Road_125m",
                #"alpha0.D_Road_125m",
                #"alpha.Dev_lo.D_Road_125m", "alpha.Dev_bg.D_Road_125m",
                #"shape.D_Road_125m",
                "alpha0.Road_125m",
                "alpha.Dev_lo.Road_125m", "alpha.Dev_bg.Road_125m",
                "shape.Road_125m",
                
                "alpha0.AHerb", "alpha.WellA_125m.AHerb",
                "alpha.Road_125m.AHerb",
                
                "test.WellD_3km", "test.WellA_1km", "test.Road_1km",
                "test.WellA_125m", "test.Road_125m", "test.AHerb")

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

# Fit model
st.time <- Sys.time()
if(package == "jagsUI") {
  out <- jagsUI(data, inits, parameters.to.save = parameters, model.file, n.thin=nt, n.chains=nc, # jagsUI abandoned after multiple crashes
                n.burnin=nb, n.iter=ni, parallel=TRUE)
} else {
  if(!file.exists(str_c(saveJAGS.loc, mod.nam))) dir.create(str_c(saveJAGS.loc, mod.nam))
  out <- saveJAGS(data = data, inits = inits, params = parameters, modelFile = model.file, thin = nt, chains = nc, # taking up saveJAGS instead to save as we go.
                  burnin = nb, sample2save = ((ni/nt)/ns), nSaves = ns, fileStub = str_c(saveJAGS.loc, mod.nam, "/modsave"))
}
end.time <- Sys.time()
run.time <- end.time - st.time
run.time
rm(st.time,end.time)

if(package == "jagsUI") {
  R.utils::saveObject(out, mod.nam)
}
#out <- resumeJAGS(fileStub = str_c(saveJAGS.loc, mod.nam, "/modsave"), nSaves = 40)

# Save summaries for simple model #
out$summary[which(!(str_sub(row.names(out$summary), 1, 6) %in%
                      c("BETA0[", "beta0[", "zeta0["))),] %>%
  round(digits = 5) %>%
  write.csv("Mod_sum_simple.csv")

# Gather, combine, and summarize JAGS saves from hard drive #
rsav <- recoverSaves(str_c(saveJAGS.loc, mod.nam, "/modsave"))
mod.raw <- combineSaves(rsav, burnFiles = 0, thin = 1)
mod <- mcmcOutput(mod.raw)
sumTab <- summary(mod, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = F)
sumTab <- sumTab %>%
  as_tibble() %>%
  mutate(Parameter = row.names(sumTab)) %>%
  select(Parameter, mean:f)
mod <- list(mcmcOutput = mod, sims.list = simsList(mod.raw), summary = sumTab)

# Check the basics
#max(mod$summary[-which(str_sub(dimnames(mod$summary)[[1]], 1, 2) == "z["), "Rhat"], na.rm = T)
max(mod$summary[, "Rhat"], na.rm = T)
min(mod$summary[, "n.eff"], na.rm = T)
#mod$summary %>% slice(order(mod$summary$Rhat, decreasing = T)[1:50]) %>% View
#traceplot(mod.raw[, "sigma.alpha1[1]"])

# traceplots (***not sure this will work if z's are saved and included.) #
pdf(file=str_c(mod.nam, '_traceplots.pdf'))
plot.params <- params.saved <- parameters
for(i in 1:length(plot.params)) {
  par.i <- plot.params[i]
  pars.lst <- params.saved[which(substring(params.saved,1,nchar(par.i))==par.i)]
  pars.cols <- c()
  if(any(substring(colnames(mod.raw$AA),1,nchar(par.i))==par.i)) {
    for(j in 1:length(pars.lst)) {
      pars.cols <- c(pars.cols,which(substring(colnames(mod.raw$AA),1,nchar(pars.lst[j]))==pars.lst[j]))
    }

    if(length(pars.cols)<=9) {
      par(mfrow=c(ceiling(sqrt(length(pars.cols))),ceiling(sqrt(length(pars.cols)))),xpd=NA,mar=c(5,4,1,1))
      for(j in pars.cols) {
        matplot(cbind(mod.raw$AA[,j],mod.raw$AB[,j],mod.raw$AC[,j]),type='l',lty=1,ylab=colnames(mod.raw$AA)[j])
      }
    }
    if(length(pars.cols)>9) {
      for(j in 1:length(pars.cols)) {
        if((j-1)%%9==0) {par(mfrow=c(3,3),xpd=NA,mar=c(5,4,1,1))}
        matplot(cbind(mod.raw$AA[,pars.cols[j]],mod.raw$AB[,pars.cols[j]],mod.raw$AC[,pars.cols[j]]),type='l',lty=1,ylab=colnames(mod.raw$AA)[pars.cols[j]])
      }
    }
  }
}
dev.off()

# Save output
library(R.utils)
saveObject(mod, mod.nam)
