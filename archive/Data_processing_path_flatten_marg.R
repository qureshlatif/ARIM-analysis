# Primary data objects
PSI.vars.cntrl <- c("PJ_area", "NDVI")
PSI.vars.trt <- c("Dev_bg", "Dev_lo", "Well_3km", "Road_1km")

beta.vars.cntrl <- c("TPI_min", "Sage", "Herb")
beta.vars.trt <- c("Dev_bg", "Dev_lo", "Well_1km", "Well_125m", "Road_125m", "AHerb")

psi_dyn.vars.cntrl <- c("Sage", "Herb")
psi_dyn.vars.trt <- c("Dev_bg", "Dev_lo", "Well_1km", "Well_125m", "Road_125m", "AHerb")

zeta.vars <- c("CanCov", "ShrubCov", "DOY", "Time_ssr")
zeta.quad <- c(F, F, T, T)

# Flatten and attach point-level covariate values #
Cov_point_flat <- matrix(NA, nrow = nrow(cov_pntyr), ncol = dim(Cov_point)[3],
                         dimnames = list(NULL, dimnames(Cov_point)[[3]]))
for(i in 1:nrow(cov_pntyr)) Cov_point_flat[i,] <-
  Cov_point[cov_pntyr$pointIndex[i], cov_pntyr$YearInd[i], ]

# Detection data #
Y <- Y.mat
Y.point <- Y
TPeriod <- TR.mat
gridID <- cov_pntyr[, "gridIndex"]
yearID <- cov_pntyr[, "YearInd"]
pointID <- cov_pntyr[, "pointIndex"]
n.grid <- max(gridID)
n.year <- max(yearID)
n.point <- max(pointID)
n.pntyr <- dim(Y)[1]
n.spp <- dim(Y)[2]

gridXyrID <- str_c(gridID, yearID, sep = "_") %>% as.factor() %>% as.integer()
Y.grid <- apply(Y, 2, function(x)
  tapply(x, gridXyrID, function(y) any(y == 1) * 1))
Y.spp <- apply(Y, 2, function(x) any(x == 1) * 1)

# Guild matrix #
guilds <- unique(spp.out$Guild)
n.guild <- length(guilds)
guildMem <- matrix(0, nrow = length(spp.list), ncol = n.guild,
                   dimnames = list(spp.list, guilds))
for(g in 1:n.guild) guildMem[which(spp.out$Guild == guilds[g]), g] <- 1
guildMem[which(spp.out$Guild == "Sagebrush"), "Shrubland"] <- 1 # Sagebrush dependent species are also shrubland species.

# Reduce data augmentation as selected #
if(reduce.data.aug) {
  ind.rm <- which(apply(Y, 2, function(x) !any(x == 1))) %>%
    (function(x) x[11:length(x)])
  ind.spp <- (1:n.spp)[-ind.rm]
  Y <- Y[,ind.spp]
  Y.grid <- Y.grid[,ind.spp]
  Y.spp <- Y.spp[ind.spp]
  #spp.list <- spp.list[ind.spp]
  n.spp <- dim(Y)[2]
  guildMem <- guildMem[ind.spp,]
  rm(ind.rm, ind.spp)
}

# Subset for development runs #
if(development) {
  ind.spp <- which(spp.list %in% develop.spp)
  Y <- Y[,ind.spp]
  Y.grid <- Y.grid[,ind.spp]
  Y.spp <- Y.spp[ind.spp]
  #spp.list <- spp.list[ind.spp]
  n.spp <- dim(Y)[2]
  guildMem <- guildMem[ind.spp,]
  rm(ind.spp)
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
dimnames(X.PSI)[[3]] <- vars
ind.PSI.offset <- which(vars %in% PSI.vars.trt)
ind.PSI.no_offset <- which(vars %in% PSI.vars.cntrl)

    # Additional indices for intermediate path models #
ind.PSI.Dev_bg <- which(vars == "Dev_bg")
ind.PSI.Dev_lo <- which(vars == "Dev_lo")
ind.Well_3km <- which(vars == "Well_3km")
ind.Road_1km <- which(vars == "Road_1km")
X.PSI.raw[,,"Road_1km"] <- X.PSI.raw[,,"Road_1km"] + 0.01

X.PSI.raw <- abind::abind(X.PSI.raw, Cov_grid[,,"Well_1km"], along = 3)
ind.Well_1km <- dim(X.PSI.raw)[3]

  # Point
vars <- c(beta.vars.cntrl, beta.vars.trt)
X.psi.raw <- array(NA, dim = c(dim(Cov_point_flat)[1], length(vars)))
dimnames(X.psi.raw)[[2]] <- vars
ind.psi.dyn <- which(vars %in% c(psi_dyn.vars.cntrl, psi_dyn.vars.trt))
ind.psi.point.vars <- which(vars %in% dimnames(Cov_point_flat)[[2]])
X.psi.raw[,ind.psi.point.vars] <- Cov_point_flat[,vars[ind.psi.point.vars]]
for(k in unique(gridID)) for(t in unique(yearID)) {
  vals <- Cov_grid[k,t,vars[-ind.psi.point.vars]] %>%
    array(dim = c(length(vars[-ind.psi.point.vars]), sum(gridID == k & yearID == t))) %>%
    aperm(perm = c(2, 1))
  X.psi.raw[which(gridID == k & yearID == t),-ind.psi.point.vars] <- vals
}
X.psi.raw[,"Well_125m"] <- (X.psi.raw[,"Well_125m"] > 0) * 1 # Convert 125 m scale well pad count to presence/absence.
X.psi <- X.psi.raw %>% apply(2, function(x) (x - mean(x)) / sd(x))
p.psi.init <- length(vars)
p.psi.dyn <- length(c(psi_dyn.vars.cntrl, psi_dyn.vars.trt))
ind.psi.init.offset <- which(vars %in% beta.vars.trt)
ind.psi.init.no_offset <- which(vars %in% beta.vars.cntrl)
ind.psi.dyn.offset <- which(vars[ind.psi.dyn] %in% psi_dyn.vars.trt)
ind.psi.dyn.no_offset <- which(vars[ind.psi.dyn] %in% psi_dyn.vars.cntrl)

X.psi.dyn <- X.psi[, ind.psi.dyn]

    # Additional variables for intermediate path models #
ind.Well_125m <- which(vars == "Well_125m")
ind.Road_125m <- which(vars == "Road_125m")

D_Road_125m <- replace(X.psi.raw[,"Road_125m"],
                       which(X.psi.raw[,"Road_125m"] == 0),
                       NA)
X.psi.raw <- cbind(X.psi.raw, D_Road_125m)
ind.D_Road_125m <- dim(X.psi.raw)[2]

PA_Road_125m <- (X.psi.raw[,"Road_125m"] > 0) * 1
X.psi.raw <- cbind(X.psi.raw, PA_Road_125m)
ind.PA_Road_125m <- dim(X.psi.raw)[2]

X.psi.raw[,"AHerb"] <- ((X.psi.raw[,"AHerb"] + 0.001) / 100)
ind.AHerb <- which(vars == "AHerb")

  # Detection
X.zeta <- cov_pntyr[, zeta.vars] %>%
  data.matrix %>%
  apply(2, (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)))
for(v in 1:length(zeta.vars)) {
  if(zeta.quad[v]) {
    zeta.vars <- c(zeta.vars, str_c(zeta.vars[v], "2"))
    X.zeta <- cbind(X.zeta, X.zeta[,v]^2)
  }
}
dimnames(X.zeta)[[2]] <- zeta.vars
X.zeta[which(is.na(X.zeta))] <- 0
p.zeta <- ncol(X.zeta)
rm(v)
