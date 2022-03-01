library(stringr)
library(tidyr)
library(dplyr)

setwd("C:/Users/Quresh.Latif/files/projects/BLM/ARIM")
load("Data_compiled.RData")

#_____ Script inputs _____#
scripts.loc <- "ARIM-analysis/"
model.file <- str_c(scripts.loc, "model.jags")
saveJAGS.loc <- "saveJAGS/"
package <- "jagsUI" # Set to jagsUI or saveJAGS
mod.nam <- "mod_trt"

# MCMC values
nc <- 2 #3 # number of chains
nb <- 10#1000 # burn in
ni <- 20#11000 # number of iterations
nt <- 1 #10 # thinning
#_________________________#

if(package == "jagsUI") {
  library(jagsUI)
} else {
  library(saveJAGS)
}

# Compile data #
source(str_c(scripts.loc, "Data_processing.R"))

# Data objects to send to JAGS
data.nams <- list("Y", "TPeriod", "gridID", "yearID", "pointID", "n.grid", "n.year", "n.point",
                  "n.pntyr", "n.spp", "X.beta", "p.beta", "tpi", "X.delta", "p.delta",
                  "X.eta", "p.eta", "X.BETA", "p.BETA", "X.DELTA", "p.DELTA",
                  "X.ETA", "p.ETA", "X.zeta", "p.zeta")

# Stuff to save from JAGS
parameters <- c("omega", "rho.zb", "rho.bB",
                
                "BETA0.mu", "sigma.BETA0", "BETA1.mu", "sigma.BETA1",
                "DELTA0.mu", "sigma.DELTA0", "sigma.D0", "DELTA1.mu", "sigma.DELTA1",
                "ETA0.mu", "sigma.ETA0", "sigma.E0", "ETA1.mu", "sigma.ETA1",
                "beta0.mu", "sigma.beta0", "beta1.mu", "sigma.beta1",
                "delta0.mu", "sigma.delta0", "sigma.d0", "delta1.mu", "sigma.delta1",
                "eta0.mu", "sigma.eta0", "sigma.e0", "eta1.mu", "sigma.eta1",
                "zeta0.mu", "sigma.zeta0", "sigma.z0", "zeta1.mu", "sigma.zeta1",
                
                "BETA0", "BETA1",
                "DELTA0", "dev.DELTA", "DELTA1",
                "ETA0", "dev.ETA", "ETA1",
                "beta0", "beta1", "beta.gully", "tpi.threshold",
                "delta0", "dev.delta", "delta1",
                "eta0", "dev.eta", "eta1",
                "zeta0", "dev.zeta", "zeta1")

# Function for setting initial values in JAGS
inits <- function()
  list(z=z.init, Z=Z.init, w=w.init,
       tvar.sigma.zeta0 = rnorm(1),
       tvar.sigma.z0 = rnorm(1),
       tvar.zeta1 = rnorm(p.zeta),
       
       tvar.sigma.BETA0 = rnorm(1),
       tvar.BETA1 = rnorm(p.BETA),
       
       tvar.sigma.DELTA0 = rnorm(1),
       tvar.sigma.D0 = rnorm(1),
       tvar.DELTA1 = rnorm(p.DELTA),
       
       tvar.sigma.ETA0 = rnorm(1),
       tvar.sigma.E0 = rnorm(1),
       tvar.ETA1 = rnorm(p.ETA),
       
       tvar.sigma.beta0 = rnorm(1),
       tvar.beta1 = rnorm(p.beta),
       tvar.beta.gully = rnorm(1),
       
       tvar.sigma.delta0 = rnorm(1),
       tvar.sigma.d0 = rnorm(1),
       tvar.delta1 = rnorm(p.delta),
       
       tvar.sigma.eta0 = rnorm(1),
       tvar.sigma.e0 = rnorm(1),
       tvar.eta1 = rnorm(p.eta))

# Assemble the initial values for JAGS.
z.init <- array(1, dim = c(n.point, n.year, n.spp))
Z.init <- array(1, dim = c(n.grid, n.year, n.spp))
w.init <- rep(1, n.spp)

data <- list()
for(i in 1:length(data.nams)) data[[length(data) + 1]] <- eval(as.name(data.nams[[i]]))
names(data) <- data.nams

# Fit model
st.time <- Sys.time()
if(package == "jagsUI") {
  out <- jagsUI(data, inits, parameters.to.save = parameters, model.file, n.thin=nt, n.chains=nc, # jagsUI abandoned after multiple crashes
                n.burnin=nb, n.iter=ni, parallel=TRUE)
} else {
  if(!file.exists(str_c(saveJAGS.loc, mod.nam))) dir.create(str_c(saveJAGS.loc, mod.nam))
  out <- saveJAGS(data = data, inits = inits, params = parameters, modelFile = model.file, thin = nt, chains = nc, # taking up saveJAGS instead to save as we go.
                  burnin = nb, sample2save = ((ni/nt)/50), nSaves = 50, fileStub = str_c(saveJAGS.loc, mod.nam, "/modsave"))
}
end.time <- Sys.time()
run.time <- end.time - st.time
run.time
rm(st.time,end.time)

#out <- resumeJAGS(fileStub = str_c(saveJAGS.loc, mod.nam, "/modsave"), nSaves = 40)

# Gather, combine, and summarize JAGS saves from hard drive #
rsav <- recoverSaves(str_c(saveJAGS.loc, mod.nam, "/modsave"))
mod.raw <- combineSaves(rsav, burnFiles = 20, thin = 10)
mod <- mcmcOutput(mod.raw)
sumTab <- summary(mod, MCEpc = F, n.eff = T, f = T, overlap0 = T, verbose = F)
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
