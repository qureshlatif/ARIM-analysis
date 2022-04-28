# Primary data objects
PSI.vars <- c("PJ_area", "NDVI", "Dev_bg", "Dev_lo", "Well_3km", "Road_1km")
LAMBDA.vars <- c("Dev_bg", "Dev_lo", "Well_3km")
  
psi.vars <- c("TPI_min", "Sage", "Herb", "Dev_bg", "Dev_lo",
              "Well_1km", "Road_125m", "AHerb")
lambda.vars <- c("Dev_bg", "Dev_lo", "Well_1km", "Road_125m", "AHerb")

zeta.vars <- c("ShrubCov", "DOY", "Time_ssr")
zeta.quad <- c(F, T, F)

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
vars <- unique(c(PSI.vars, LAMBDA.vars))
p.PSI <- length(PSI.vars)
p.LAMBDA <- length(LAMBDA.vars)
X.raw <- Cov_grid[,,vars]
X <- X.scale.fn(X.raw)
X.mns <- apply(X.raw, 3, mean, na.rm = T) # Save for unscaling later
X.sd <- apply(X.raw, 3, sd, na.rm = T) # Save for unscaling later
dimnames(X)[[3]] <- vars

  # Additional indices for intermediate path models #
ind.PSI.Dev_bg <- which(vars == "Dev_bg")
ind.PSI.Dev_lo <- which(vars == "Dev_lo")
ind.PSI.Well_3km <- which(vars == "Well_3km")
ind.PSI.Road_1km <- which(vars == "Road_1km")
X.raw[,,"Road_1km"] <- X.raw[,,"Road_1km"] + 0.01

X.raw <- abind::abind(X.raw, Cov_grid[,,"Well_1km"], along = 3)
ind.PSI.Well_1km <- dim(X.raw)[3]
dimnames(X.raw)[[3]][ind.PSI.Well_1km] <- "Well_1km"

    # Flatten to grdXyr
X.flat <- matrix(NA, nrow = n.grdyr, ncol = dim(X)[3])
X.raw.flat <- matrix(NA, nrow = n.grdyr, ncol = dim(X.raw)[3])
for(j in 1:n.grdyr) {
  X.flat[j, ] <- X[gridID.grdyr[j], yearID.grdyr[j],]
  X.raw.flat[j, ] <- X.raw[gridID.grdyr[j], yearID.grdyr[j],]
}
dimnames(X.flat)[[2]] <- dimnames(X)[[3]]
dimnames(X.raw.flat)[[2]] <- dimnames(X.raw)[[3]]
X <- X.flat
X.raw <- X.raw.flat
X.PSI <- X[,PSI.vars]
X.LAMBDA <- X[,LAMBDA.vars]
X.PSI.raw <- X.raw
rm(X.flat, X.raw.flat, X, X.raw)

  # Point
vars <- unique(c(psi.vars, lambda.vars))
p.psi <- length(psi.vars)
p.lambda <- length(lambda.vars)
X.raw <- array(NA, dim = c(dim(Cov_point_flat)[1], length(vars)))
dimnames(X.raw)[[2]] <- vars
ind.point.vars <- which(vars %in% dimnames(Cov_point_flat)[[2]])
X.raw[,ind.point.vars] <- Cov_point_flat[,vars[ind.point.vars]]
for(k in unique(gridID.py)) for(t in unique(yearID.py)) {
  vals <- Cov_grid[k,t,vars[-ind.point.vars]] %>%
    array(dim = c(length(vars[-ind.point.vars]),
                  sum(gridID.py == k & yearID.py == t))) %>%
    aperm(perm = c(2, 1))
  X.raw[which(gridID.py == k & yearID.py == t),-ind.point.vars] <-
    vals
}
X.raw <- apply(X.raw, 2, function(x) tapply(x, grdyrID.py, mean))
X <- X.raw %>% apply(2, function(x) (x - mean(x)) / sd(x))
X.psi <- X[,psi.vars]
X.lambda <- X[,lambda.vars]
X.psi.raw <- X.raw
rm(X, X.raw)

ind.psi.Well_1km <- which(vars == "Well_1km")
  
X.psi.raw[,"Road_125m"] <- X.psi.raw[,"Road_125m"] + 0.01
ind.psi.Road_125m <- which(vars == "Road_125m")
  
X.psi.raw[,"AHerb"] <- ((X.psi.raw[,"AHerb"] + 0.001) / 100)
ind.psi.AHerb <- which(vars == "AHerb")

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