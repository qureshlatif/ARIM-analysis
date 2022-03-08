# Primary data objects
PSI.vars.cntrl <- c("PJ_area", "NDVI")
PSI.vars.trt <- c("Dev_bg", "Dev_lo", "WellA_3km",
                   "WellD_3km", "Road_1km")
PSI.additional.path <- c("WellA_1km", "WellD_1km")

beta.vars.cntrl <- c("vrm_125m", "TPI_min", "Sage", "Herb")
beta.vars.trt <- c("Dev_bg", "Dev_lo", "WellA_125m", "Road_125m",
                   "AHerb", "WellA_1km", "Road_1km")

psi_dyn.vars.cntrl <- c("Sage", "Herb")
psi_dyn.vars.trt <- c("Dev_bg", "Dev_lo", "WellA_125m", "Road_125m",
                      "AHerb", "WellA_1km", "Road_1km")

zeta.vars <- c("CanCov", "ShrubCov", "DOY", "Time_ssr")
zeta.quad <- c(F, F, T, T)

# Detection data #
Y <- Y.mat
TPeriod <- TR.mat
gridID <- cov_point[, "gridIndex"]
yearID <- cov_pntyr[, "YearInd"]
pointID <- cov_pntyr[, "pointIndex"]
n.grid <- max(gridID)
n.year <- max(yearID)
n.point <- max(pointID)
n.pntyr <- dim(Y)[1]
n.spp <- dim(Y)[2]

# Guild matrix #
guilds <- unique(spp.out$Guild)
n.guild <- length(guilds)
guildMem <- matrix(0, nrow = length(spp.list), ncol = n.guild,
                   dimnames = list(spp.list, guilds))
for(g in 1:n.guild) guildMem[which(spp.out$Guild == guilds[g]), g] <- 1
guildMem[which(spp.out$Guild == "Sagebrush"), "Shrubland"] <- 1 # Sagebrush dependent species are also shrubland species.

# Subset for development runs #
if(development) {
  ind.spp <- which(spp.list %in% develop.spp)
  Y <- Y[,ind.spp]
  #spp.list <- spp.list[ind.spp]
  n.spp <- dim(Y)[2]
  guildMem <- guildMem[ind.spp,]
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
vars <- c(PSI.vars.cntrl, PSI.vars.trt, PSI.additional.path)
p.PSI <- length(PSI.vars.cntrl) + length(PSI.vars.trt)
X.PSI.raw <- Cov_grid[,,vars]
X.PSI <- X.scale.fn(X.PSI.raw)
dimnames(X.PSI)[[3]] <- vars
ind.PSI.offset <- which(vars %in% PSI.vars.trt)
ind.PSI.no_offset <- which(vars %in% PSI.vars.cntrl)

    # Additional indices for intermediate path models #
ind.PSI.Dev_bg <- which(vars == "Dev_bg")
ind.PSI.Dev_lo <- which(vars == "Dev_lo")
ind.WellA_3km <- which(vars == "WellA_3km")
ind.WellD_3km <- which(vars == "WellD_3km")
ind.WellA_1km <- which(vars == "WellA_1km")
ind.Road_1km <- which(vars == "Road_1km")
X.PSI.raw[,,"Road_1km"] <- X.PSI.raw[,,"Road_1km"] + 0.01

  # Point
vars <- c(beta.vars.cntrl, beta.vars.trt)
X.psi.raw <- array(NA, dim = c(dim(Cov_point)[1:2], length(vars)))
dimnames(X.psi.raw)[[2]] <- dimnames(Cov_point)[[2]]
dimnames(X.psi.raw)[[3]] <- vars
ind.psi.tpi <- which(vars == "TPI_min")
ind.psi.not_tpi <- which(vars != "TPI_min")
ind.psi.dyn <- which(vars %in% c(psi_dyn.vars.cntrl, psi_dyn.vars.trt))
ind.psi.point.vars <- which(vars %in% dimnames(Cov_point)[[3]])
X.psi.raw[,,ind.psi.point.vars] <- Cov_point[,,vars[ind.psi.point.vars]]
X.psi.raw[,,-ind.psi.point.vars] <- Cov_grid[gridID,,vars[-ind.psi.point.vars]]
X.psi.raw[,,"WellA_125m"] <- (X.psi.raw[,,"WellA_125m"] > 0) * 1 # Convert 125 m scale active well pad count to presence/absence.
X.psi <- X.scale.fn(X.psi.raw)
p.psi.init <- length(vars)
p.psi.dyn <- length(c(psi_dyn.vars.cntrl, psi_dyn.vars.trt))
ind.psi.init.offset <- which(vars %in% beta.vars.trt)
ind.psi.init.no_offset <- which(vars %in% beta.vars.cntrl)
ind.psi.dyn.offset <- which(vars[ind.psi.dyn] %in% psi_dyn.vars.trt)
ind.psi.dyn.no_offset <- which(vars[ind.psi.dyn] %in% psi_dyn.vars.cntrl)

    # Additional variables for intermediate path models #
ind.WellA_125m <- which(vars == "WellA_125m")
ind.Road_125m <- which(vars == "Road_125m")

D_Road_125m <- replace(X.psi.raw[,,"Road_125m"],
                       which(X.psi.raw[,,"Road_125m"] == 0),
                       NA)
X.psi.raw <- abind::abind(X.psi.raw, D_Road_125m, along = 3)
ind.D_Road_125m <- dim(X.psi.raw)[3]

PA_Road_125m <- (X.psi.raw[,,"Road_125m"] > 0) * 1
X.psi.raw <- abind::abind(X.psi.raw, PA_Road_125m, along = 3)
ind.PA_Road_125m <- dim(X.psi.raw)[3]

X.psi.raw[,,"AHerb"] <- ((X.psi.raw[,,"AHerb"] + 0.001) / 100)
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
