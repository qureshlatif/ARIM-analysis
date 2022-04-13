# Primary data objects
if(mod.nam %in% c("mod_path", "mod_interm_paths")) {
  PSI.vars.cntrl <- c("PJ_area", "NDVI")
  PSI.vars.trt <- c("Dev_bg", "Dev_lo", "Well_3km", "Road_1km")
  
  psi.vars.cntrl <- c("TPI_min", "Sage", "Herb")
  psi.vars.trt <- c("Dev_bg", "Dev_lo", "Well_1km", "Road_125m", "AHerb")
}

# if(mod.nam == "mod_community_trend") {
#   PSI.vars.cntrl <- c("PJ_area", "NDVI")
#   PSI.vars.trt <- c("Dev_bg", "Dev_lo")
#   
#   beta.vars.cntrl <- c("TPI_min", "Sage", "Herb")
#   beta.vars.trt <- c("Dev_bg", "Dev_lo")
#   
#   psi_dyn.vars.cntrl <- c("Sage", "Herb")
#   psi_dyn.vars.trt <- c("Dev_bg", "Dev_lo")
# }
# 
# if(mod.nam == "mod_community_mech") {
#   PSI.vars.cntrl <- c("PJ_area", "NDVI")
#   PSI.vars.trt <- c("Well_3km", "Road_1km")
#   
#   beta.vars.cntrl <- c("TPI_min", "Sage", "Herb")
#   beta.vars.trt <- c("Well_1km", "Road_125m", "AHerb")
#   
#   psi_dyn.vars.cntrl <- c("Sage", "Herb")
#   psi_dyn.vars.trt <- c("Well_1km", "Road_125m", "AHerb")
# }

zeta.vars <- c("CanCov", "ShrubCov", "DOY", "Time_ssr")
zeta.quad <- c(F, F, T, T)

# Flatten and attach point-level covariate values #
Cov_point_flat <- matrix(NA, nrow = nrow(cov_pntyr), ncol = dim(Cov_point)[3],
                         dimnames = list(NULL, dimnames(Cov_point)[[3]]))
for(i in 1:nrow(cov_pntyr)) Cov_point_flat[i,] <-
  Cov_point[cov_pntyr$pointIndex[i], cov_pntyr$YearInd[i], ]

# Detection data #
TPeriod <- TR.mat
gridID.py <- cov_pntyr[, "gridIndex"]
yearID.py <- cov_pntyr[, "YearInd"]
grdyrID.py <- str_c(gridID.py, yearID.py, sep = "_") %>%
  as.factor() %>% as.integer()
observerID.py <- cov_pntyr[, "ObserverID"]
n.grdyr <- max(grdyrID.py)
n.spp <- dim(Y.mat)[2]
K <- max(TPeriod)
Y <- array(NA, dim = c(n.grdyr, n.spp, K))
for(j in 1:n.grdyr) for(k in 1:K)
  Y[j,,k] <-
  apply(Y.mat[which(grdyrID.py == j),] *
          (TPeriod[which(grdyrID.py == j),] == k), 2, sum)
Y.sum <- apply(Y, c(1, 2), sum)
gridID.grdyr <- tapply(gridID.py, grdyrID.py, unique)
yearID.grdyr <- tapply(yearID.py, grdyrID.py, unique)
observerID <- tapply(observerID.py, grdyrID.py, unique)
n.year <- max(yearID.py)
n.point <- tapply(grdyrID.py, grdyrID.py, length)
n.observer <- max(observerID)

  # Detections within survey intervals #
n.int <- sum(Y.sum > 0)
grdyrID.int.long <- rep(1:n.grdyr, n.spp)
sppID.int.long <- rep(1:n.spp, each = n.grdyr)
Y.long <- Y[,1,]
for(i in 2:n.spp) {
  Y.long <- rbind(Y.long, Y[,i,])
}
ind.det <- which(Y.sum > 0)
grdyrID.int <- grdyrID.int.long[ind.det]
sppID.int <- sppID.int.long[ind.det]
Y.int <- Y.long[ind.det,]
rm(ind.det, Y.long, grdyrID.int.long, sppID.int.long, i)

# # Guild matrix #
# guilds <- unique(spp.out$Guild)
# n.guild <- length(guilds)
# guildMem <- matrix(0, nrow = length(spp.list), ncol = n.guild,
#                    dimnames = list(spp.list, guilds))
# for(g in 1:n.guild) guildMem[which(spp.out$Guild == guilds[g]), g] <- 1
# guildMem[which(spp.out$Guild == "Sagebrush"), "Shrubland"] <- 1 # Sagebrush dependent species are also shrubland species.

# Subset for development runs #
if(development) {
  ind.spp <- which(spp.list %in% develop.spp)
  Y <- Y[,ind.spp,]
  n.spp <- dim(Y)[2]
  # guildMem <- guildMem[ind.spp,]

  Y.sum <- apply(Y, c(1, 2), sum)
  n.int <- sum(Y.sum > 0)
  grdyrID.int.long <- rep(1:n.grdyr, n.spp)
  sppID.int.long <- rep(1:n.spp, each = n.grdyr)
  Y.long <- Y[,1,]
  for(i in 2:n.spp) {
    Y.long <- rbind(Y.long, Y[,i,])
  }
  ind.det <- which(Y.sum > 0)
  grdyrID.int <- grdyrID.int.long[ind.det]
  sppID.int <- sppID.int.long[ind.det]
  Y.int <- Y.long[ind.det,]
  rm(ind.det, Y.long, grdyrID.int.long, sppID.int.long, i)
}

# Compile covariates #
X.scale.fn <- function(X) {
  mns <- apply(X, 3, mean, na.rm = T) %>%
    array(c(dim(X)[c(3, 1, 2)])) %>% aperm(perm = c(2, 3, 1))
  sds <- apply(X, 3, sd, na.rm = T) %>%
    array(c(dim(X)[c(3, 1, 2)])) %>% aperm(perm = c(2, 3, 1))
  X <- (X - mns) / sds
  return(X)
}

  # Grid cell
vars <- c(PSI.vars.cntrl, PSI.vars.trt)
p.PSI <- length(PSI.vars.cntrl) + length(PSI.vars.trt)
X.PSI.raw <- Cov_grid[,,vars]
X.PSI <- X.scale.fn(X.PSI.raw)
X.PSI.mns <- apply(X.PSI.raw, 3, mean, na.rm = T) # Save for unscaling later
X.PSI.sd <- apply(X.PSI.raw, 3, sd, na.rm = T) # Save for unscaling later
dimnames(X.PSI)[[3]] <- vars
#ind.PSI.offset <- which(vars %in% PSI.vars.trt)
#ind.PSI.no_offset <- which(vars %in% PSI.vars.cntrl)

    # Additional indices for intermediate path models #
if(mod.nam %in% c("mod_path", "mod_interm_paths")) {
  ind.PSI.Dev_bg <- which(vars == "Dev_bg")
  ind.PSI.Dev_lo <- which(vars == "Dev_lo")
  ind.PSI.Well_3km <- which(vars == "Well_3km")
  ind.PSI.Road_1km <- which(vars == "Road_1km")
  X.PSI.raw[,,"Road_1km"] <- X.PSI.raw[,,"Road_1km"] + 0.01

  X.PSI.raw <- abind::abind(X.PSI.raw, Cov_grid[,,"Well_1km"], along = 3)
  ind.PSI.Well_1km <- dim(X.PSI.raw)[3]
  dimnames(X.PSI.raw)[[3]][ind.PSI.Well_1km] <- "Well_1km"
  }

    # Flatten to grdXyr
X.PSI.flat <- matrix(NA, nrow = n.grdyr, ncol = dim(X.PSI)[3])
X.PSI.raw.flat <- matrix(NA, nrow = n.grdyr, ncol = dim(X.PSI.raw)[3])
for(j in 1:n.grdyr) {
  X.PSI.flat[j, ] <- X.PSI[gridID.grdyr[j], yearID.grdyr[j],]
  X.PSI.raw.flat[j, ] <- X.PSI.raw[gridID.grdyr[j], yearID.grdyr[j],]
}
dimnames(X.PSI.flat)[[2]] <- dimnames(X.PSI)[[3]]
dimnames(X.PSI.raw.flat)[[2]] <- dimnames(X.PSI.raw)[[3]]
X.PSI <- X.PSI.flat
X.PSI.raw <- X.PSI.raw.flat
rm(X.PSI.flat, X.PSI.raw.flat)

X.LAMBDA <- X.PSI[,PSI.vars.trt]
p.LAMBDA <- length(PSI.vars.trt)

  # Point
vars <- c(psi.vars.cntrl, psi.vars.trt)
X.psi.raw <- array(NA, dim = c(dim(Cov_point_flat)[1], length(vars)))
dimnames(X.psi.raw)[[2]] <- vars
#ind.psi.dyn <- which(vars %in% c(psi_dyn.vars.cntrl, psi_dyn.vars.trt))
ind.psi.point.vars <- which(vars %in% dimnames(Cov_point_flat)[[2]])
X.psi.raw[,ind.psi.point.vars] <- Cov_point_flat[,vars[ind.psi.point.vars]]
for(k in unique(gridID.py)) for(t in unique(yearID.py)) {
  vals <- Cov_grid[k,t,vars[-ind.psi.point.vars]] %>%
    array(dim = c(length(vars[-ind.psi.point.vars]), sum(gridID.py == k & yearID.py == t))) %>%
    aperm(perm = c(2, 1))
  X.psi.raw[which(gridID.py == k & yearID.py == t),-ind.psi.point.vars] <- vals
}
# if(mod.nam %in% c("mod_path", "mod_community_mech", "mod_interm_paths")) # Not needed after dropping Well_125m
#   X.psi.raw[,"Well_125m"] <- (X.psi.raw[,"Well_125m"] > 0) * 1 # Convert 125 m scale well pad count to presence/absence.
X.psi.raw <- apply(X.psi.raw, 2, function(x) tapply(x, grdyrID.py, mean))
# if(mod.nam %in% c("mod_path", "mod_community_mech", "mod_interm_paths")) # Not needed after dropping Well_125m
#   X.psi.raw[,"Well_125m"] <- X.psi.raw[,"Well_125m"] *
#     tapply(grdyrID.py, grdyrID.py, length) # Make this the sum rather than mean of point values
X.psi <- X.psi.raw %>% apply(2, function(x) (x - mean(x)) / sd(x))
p.psi <- length(vars)
X.lambda <- X.psi[,psi.vars.trt]
p.lambda <- length(psi.vars.trt)
# ind.psi.init.offset <- which(vars %in% beta.vars.trt)
# ind.psi.init.no_offset <- which(vars %in% beta.vars.cntrl)
# ind.psi.dyn.offset <- which(vars[ind.psi.dyn] %in% psi_dyn.vars.trt)
# ind.psi.dyn.no_offset <- which(vars[ind.psi.dyn] %in% psi_dyn.vars.cntrl)

#X.psi.dyn <- X.psi[, ind.psi.dyn]

if(mod.nam %in% c("mod_path", "mod_interm_paths")) {
  # Additional variables for intermediate path models #
  #ind.psi.Well_125m <- which(vars == "Well_125m") # Not needed after dropping Well_125m
  
  ind.psi.Well_1km <- which(vars == "Well_1km")
  
  X.psi.raw[,"Road_125m"] <- X.psi.raw[,"Road_125m"] + 0.01
  ind.psi.Road_125m <- which(vars == "Road_125m")
  
  X.psi.raw[,"AHerb"] <- ((X.psi.raw[,"AHerb"] + 0.001) / 100)
  ind.psi.AHerb <- which(vars == "AHerb")
}

  # Trend covariate
X.trend <- (years - mean(years)) / sd(years)

  # Detection
X.zeta <- cov_pntyr[, zeta.vars] %>%
  data.matrix %>%
  apply(2, function(x) tapply(x, grdyrID.py, mean)) %>%
  apply(2, (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)))
for(v in 1:length(zeta.vars)) {
  if(zeta.quad[v]) {
    zeta.vars <- c(zeta.vars, str_c(zeta.vars[v], "2"))
    X.zeta <- cbind(X.zeta, X.zeta[,v]^2)
  }
}
dimnames(X.zeta)[[2]] <- zeta.vars
X.zeta[which(is.na(X.zeta))] <- 0
#X.zeta <- cbind(X.zeta, X.trend[cov_pntyr$YearInd])
p.zeta <- ncol(X.zeta)
#dimnames(X.zeta)[[2]][p.zeta] <- "Trend"
rm(v)