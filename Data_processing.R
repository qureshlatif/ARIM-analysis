#*** Need to convert to presence / absence of active well pad (Only 17 points with decom well pad within 125 m).

# Primary data objects
BETA.vars.cntrl <- c("PJ_area", "NDVI")
BETA.vars.trt <- c("Develop_low", "Develop_high")
BETA.vars.mech <- c("WellA_3x3km", "WellD_3x3km", "Road_1km")

DELTA.vars.cntrl <- c("PJ_area", "NDVI")
DELTA.vars.trt <- c("Develop_low", "Develop_high")
DELTA.vars.mech <- c("WellA_3x3km", "WellD_3x3km", "Road_1km")

ETA.vars.cntrl <- c("PJ_area", "NDVI")
ETA.vars.trt <- c("Develop_low", "Develop_high")
ETA.vars.mech <- c("WellA_3x3km", "WellD_3x3km", "Road_1km")

beta.vars.cntrl <- c("vrm_125m", "Sage", "Herb")
beta.vars.trt <- c("Develop_low", "Develop_high")
beta.vars.mech <- c("WellA_125m", "Road_125m", "AHerb", "WellA_1km", "WellD_1km", "Road_1km")

delta.vars.cntrl <- c("Sage", "Herb")
delta.vars.trt <- c("Develop_low", "Develop_high")
delta.vars.mech <- c("WellA_125m", "Road_125m", "AHerb", "WellA_1km", "WellD_1km", "Road_1km")

eta.vars.cntrl <- c("Sage", "Herb")
eta.vars.trt <- c("Develop_low", "Develop_high")
eta.vars.mech <- c("WellA_125m", "Road_125m", "AHerb", "WellA_1km", "WellD_1km", "Road_1km")

zeta.vars <- c("CanCov", "ShrubCov", "DayOfYear", "Time_ssr")
zeta.quad <- c(F, F, T, T)

#p.beta <- length(beta.vars) # Number of covariates on grid-level occupancy
#p.zeta <- length(zeta.vars) + sum(zeta.quad) # Number of covariates on detection

# Detection data #
Y <- Y.mat
TPeriod <- TR.mat
gridID <- Cov_point[, "gridIndex"]
yearID <- Cov_point[, "YearInd"]
pointID <- Cov_point[, "pointInd"]
#yearID.grid <- landscape_data$YearInd
#gridID.grid <- landscape_data$gridIndex
n.grid <- max(gridID)
n.year <- max(yearID)
n.point <- max(pointID)
#n.gridXyear <- nrow(landscape_data)
n.pntyr <- dim(Y)[1]
n.spp <- dim(Y)[2]

X.scale.fn <- function(X) {
  mns <- apply(X, 3, mean, na.rm = T) %>%
    array(c(dim(X)[c(3, 1, 2)])) %>% aperm(perm = c(2, 3, 1))
  sds <- apply(X, 3, sd, na.rm = T) %>%
    array(c(dim(X)[c(3, 1, 2)])) %>% aperm(perm = c(2, 3, 1))
  X <- (X - mns) / sds
  return(X)
}

X.impute.fn <- function(X, groupID = NULL) {
  require(QSLpersonal)
  ind.all0s <- which(apply(X, 1, function(x) sum(!is.na(x))) == 0)
  if(length(ind.all0s) > 0) X[ind.all0s,] <- mean(X, na.rm = T)
  X.vec <- as.numeric(X)
  rowID <- rep(1:nrow(X), ncol(X))
  colID <- rep(1:ncol(X), each = nrow(X))
  X.row.mean <- tapply(X.vec, rowID, function(x) mean(x, na.rm = T))[rowID]
  X.col.mean <- tapply(X.vec, colID, function(x) mean(x, na.rm = T))[colID]
  if(!is.null(groupID)) {
    groupID <- rep(groupID, ncol(X))
    X.group.mean <- tapply(X.vec, groupID, function(x) mean(x, na.rm = T))[groupID]
    dat <- cbind(X.vec, X.row.mean, X.col.mean, X.group.mean)
    v.inf <- c("X.row.mean", "X.col.mean", "X.group.mean")
  } else {
    v.inf <- c("X.row.mean", "X.col.mean")
    dat <- cbind(X.vec, X.row.mean, X.col.mean)
  }
  dat <- Impute_missing_covs_rf(dat, v.fill = "X.vec",
                                v.inform = v.inf)
  X.vec.new <- dat[,"X.vec"]
  X <- matrix(X.vec.new, nrow = max(rowID), ncol = max(colID))
  return(X)
}

# Compile covariates #
  # Grid cell
if(mod.nam == "mod_trt") {
  vars <- c(BETA.vars.trt, BETA.vars.cntrl)
} else {
  vars <- c(BETA.vars.mech, BETA.vars.cntrl)
}
X.BETA <- X.scale.fn(Cov_grid[,,vars])
dimnames(X.BETA)[[3]] <- vars
p.BETA <- length(vars)
if(any(vars %in% c("Develop_low", "Develop_high")))
  for(v in vars[which(vars %in% c("Develop_low", "Develop_high"))]) {
    vals <- apply(X.BETA[,,v], 1, function(x) unique(na.omit(x)))
    X.BETA[,,v] <- vals
  }
for(v in 1:length(vars)) if(any(is.na(X.BETA[,,vars[v]]))) # Impute missing values where needed.
  X.BETA[,,vars[v]] <- X.impute.fn(X.BETA[,,vars[v]])
X.BETA <- X.BETA[,1,] # Just keep values for the first year

if(mod.nam == "mod_trt") {
  vars <- c(DELTA.vars.trt, DELTA.vars.cntrl)
} else {
  vars <- c(DELTA.vars.mech, DELTA.vars.cntrl)
}
X.DELTA <- X.scale.fn(Cov_grid[,,vars])
dimnames(X.DELTA)[[3]] <- vars
p.DELTA <- length(vars)
if(any(vars %in% c("Develop_low", "Develop_high")))
  for(v in vars[which(vars %in% c("Develop_low", "Develop_high"))]) {
    vals <- apply(X.DELTA[,,v], 1, function(x) unique(na.omit(x)))
    X.DELTA[,,v] <- vals
  }
for(v in 1:length(vars)) if(any(is.na(X.DELTA[,,vars[v]]))) # Impute missing values where needed.
  X.DELTA[,,vars[v]] <- X.impute.fn(X.DELTA[,,vars[v]])
X.DELTA <- X.DELTA[,-1,] # Just keep values for the first year

if(mod.nam == "mod_trt") {
  vars <- c(ETA.vars.trt, ETA.vars.cntrl)
} else {
  vars <- c(ETA.vars.mech, ETA.vars.cntrl)
}
X.ETA <- X.scale.fn(Cov_grid[,,vars])
dimnames(X.ETA)[[3]] <- vars
p.ETA <- length(vars)
if(any(vars %in% c("Develop_low", "Develop_high")))
  for(v in vars[which(vars %in% c("Develop_low", "Develop_high"))]) {
    vals <- apply(X.ETA[,,v], 1, function(x) unique(na.omit(x)))
    X.ETA[,,v] <- vals
  }
for(v in 1:length(vars)) if(any(is.na(X.ETA[,,vars[v]]))) # Impute missing values where needed.
  X.ETA[,,vars[v]] <- X.impute.fn(X.ETA[,,vars[v]])
X.ETA <- X.ETA[,-1,] # Just keep values for the first year

  # Point
if(mod.nam == "mod_trt") {
  vars <- c(beta.vars.trt, beta.vars.cntrl)
} else {
  vars <- c(beta.vars.mech, beta.vars.cntrl)
}
X.beta.long <- Cov_point[, vars] %>%
  apply(2, (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))) %>%
  (function(x) ifelse(is.na(x), 0, x))
X.beta <- array(NA, dim = c(n.point, n.year, length(vars)))
for(j in 1:n.pntyr) X.beta[pointID[j], yearID[j], ] <- X.beta.long[j,]
dimnames(X.beta)[[3]] <- vars
if(any(vars %in% c("Develop_low", "Develop_high"))) {
  for(v in vars[which(vars %in% c("Develop_low", "Develop_high"))]) {
    vals <- apply(X.beta[,,v], 1, function(x) unique(na.omit(x)))
    X.beta[,,v] <- vals
  }
}
grpID <- tapply(gridID, pointID, unique)
for(v in 1:length(vars)) if(any(is.na(X.beta[,,vars[v]]))) # Impute missing values where needed.
  X.beta[,,vars[v]] <- X.impute.fn(X.beta[,,vars[v]], grpID)
X.beta <- X.beta[,1,]
rm(X.beta.long, grpID)
p.beta <- dim(X.beta)[2]
tpi <- Cov_point[,"TPI_min"] %>%
  (function(x) (x - mean(x)) / sd(x)) %>%
  tapply(pointID, function(x) x[1])

if(mod.nam == "mod_trt") {
  vars <- c(delta.vars.trt, delta.vars.cntrl)
} else {
  vars <- c(delta.vars.mech, delta.vars.cntrl)
}
X.delta.long <- Cov_point[, vars] %>%
  apply(2, (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))) %>%
  (function(x) ifelse(is.na(x), 0, x))
X.delta <- array(NA, dim = c(n.point, n.year, length(vars)))
for(j in 1:n.pntyr) X.delta[pointID[j], yearID[j], ] <- X.delta.long[j,]
dimnames(X.delta)[[3]] <- vars
if(any(vars %in% c("Develop_low", "Develop_high"))) {
  for(v in vars[which(vars %in% c("Develop_low", "Develop_high"))]) {
    vals <- apply(X.delta[,,v], 1, function(x) unique(na.omit(x)))
    X.delta[,,v] <- vals
  }
}
grpID <- tapply(gridID, pointID, unique)
for(v in 1:length(vars)) if(any(is.na(X.delta[,,vars[v]]))) # Impute missing values where needed.
  X.delta[,,vars[v]] <- X.impute.fn(X.delta[,,vars[v]], grpID)
X.delta <- X.delta[,-1,]
rm(X.delta.long, grpID)
p.delta <- dim(X.delta)[3]

if(mod.nam == "mod_trt") {
  vars <- c(eta.vars.trt, eta.vars.cntrl)
} else {
  vars <- c(eta.vars.mech, eta.vars.cntrl)
}
X.eta.long <- Cov_point[, vars] %>%
  apply(2, (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))) %>%
  (function(x) ifelse(is.na(x), 0, x))
X.eta <- array(NA, dim = c(n.point, n.year, length(vars)))
for(j in 1:n.pntyr) X.eta[pointID[j], yearID[j], ] <- X.eta.long[j,]
dimnames(X.eta)[[3]] <- vars
if(any(vars %in% c("Develop_low", "Develop_high"))) {
  for(v in vars[which(vars %in% c("Develop_low", "Develop_high"))]) {
    vals <- apply(X.eta[,,v], 1, function(x) unique(na.omit(x)))
    X.eta[,,v] <- vals
  }
}
grpID <- tapply(gridID, pointID, unique)
for(v in 1:length(vars)) if(any(is.na(X.eta[,,vars[v]]))) # Impute missing values where needed.
  X.eta[,,vars[v]] <- X.impute.fn(X.eta[,,vars[v]], grpID)
X.eta <- X.eta[,-1,]
rm(X.eta.long, grpID)
p.eta <- dim(X.eta)[3]

  # Detection
X.zeta <- Cov_point[, zeta.vars] %>%
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
