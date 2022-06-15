## Grid cell ##
X.PSI.pred.hi <- X.PSI.pred.lo <- X.PSI.pred.bg <-
  matrix(NA, nrow = nsims, ncol = p.PSI)
dimnames(X.PSI.pred.hi)[[2]] <- 
  dimnames(X.PSI.pred.lo)[[2]] <- 
  dimnames(X.PSI.pred.bg)[[2]] <- 
  dimnames(X.PSI)[[2]]
X.PSI.pred.hi[, c("PJ_area", "NDVI", "Road_1km")] <-
  X.PSI.pred.lo[, c("PJ_area", "NDVI", "Road_1km")] <-
  X.PSI.pred.bg[, c("PJ_area", "NDVI", "Road_1km")] <- 0

# High development occupancy #
X.Dev <- X.PSI[which(X.PSI.raw[,"Dev_bg"] == 0 &
                       X.PSI.raw[,"Dev_lo"] == 0),
               c("Dev_lo", "Dev_bg")][1,]
X.PSI.pred.hi[, "Dev_lo"] <- X.Dev[1]
X.PSI.pred.hi[, "Dev_bg"] <- X.Dev[2]
X.PSI.pred.hi[, "Well_3km"] <- exp(mod$mcmcOutput$ALPHA0.Well_3km +
  mod$mcmcOutput$ALPHA.Dev_lo.Well_3km * X.Dev[1] +
  mod$mcmcOutput$ALPHA.Dev_bg.Well_3km * X.Dev[2]) %>%
  (function(x) (x - X.mns["Well_3km"]) / X.sd["Well_3km"])

# Low development occupancy #
X.Dev <- X.PSI[which(X.PSI.raw[,"Dev_lo"] == 1), c("Dev_lo", "Dev_bg")][1,]
X.PSI.pred.lo[, "Dev_lo"] <- X.Dev[1]
X.PSI.pred.lo[, "Dev_bg"] <- X.Dev[2]
X.PSI.pred.lo[, "Well_3km"] <- exp(mod$mcmcOutput$ALPHA0.Well_3km +
                                     mod$mcmcOutput$ALPHA.Dev_lo.Well_3km * X.Dev[1] +
                                     mod$mcmcOutput$ALPHA.Dev_bg.Well_3km * X.Dev[2]) %>%
  (function(x) (x - X.mns["Well_3km"]) / X.sd["Well_3km"])

# Background occupancy #
X.Dev <- X.PSI[which(X.PSI.raw[,"Dev_bg"] == 1), c("Dev_lo", "Dev_bg")][1,]
X.PSI.pred.bg[, "Dev_lo"] <- X.Dev[1]
X.PSI.pred.bg[, "Dev_bg"] <- X.Dev[2]
X.PSI.pred.bg[, "Well_3km"] <- exp(mod$mcmcOutput$ALPHA0.Well_3km +
                                     mod$mcmcOutput$ALPHA.Dev_lo.Well_3km * X.Dev[1] +
                                     mod$mcmcOutput$ALPHA.Dev_bg.Well_3km * X.Dev[2]) %>%
  (function(x) (x - X.mns["Well_3km"]) / X.sd["Well_3km"])

# Trend #
X.LAMBDA.pred.hi <- X.PSI.pred.hi[, dimnames(X.LAMBDA)[[2]]]
X.LAMBDA.pred.lo <- X.PSI.pred.lo[, dimnames(X.LAMBDA)[[2]]]
X.LAMBDA.pred.bg <- X.PSI.pred.bg[, dimnames(X.LAMBDA)[[2]]]


## Point-scale ##
X.psi.pred.hi <- X.psi.pred.lo <- X.psi.pred.bg <-
  matrix(NA, nrow = nsims, ncol = p.psi)
dimnames(X.psi.pred.hi)[[2]] <- 
  dimnames(X.psi.pred.lo)[[2]] <- 
  dimnames(X.psi.pred.bg)[[2]] <- 
  dimnames(X.psi)[[2]]
X.psi.pred.hi[, c("TPI_min", "Sage", "Herb")] <-
  X.psi.pred.lo[, c("TPI_min", "Sage", "Herb")] <-
  X.psi.pred.bg[, c("TPI_min", "Sage", "Herb")] <- 0

# High development occupancy #
X.Dev <- X.psi[which(X.psi.raw[,"Dev_bg"] == 0 & X.psi.raw[,"Dev_lo"] == 0),
               c("Dev_lo", "Dev_bg")][1,]
X.DEV <- X.PSI[which(X.PSI.raw[,"Dev_bg"] == 0 & X.PSI.raw[,"Dev_lo"] == 0),
               c("Dev_lo", "Dev_bg")][1,]
X.psi.pred.hi[, "Dev_lo"] <- X.Dev[1]
X.psi.pred.hi[, "Dev_bg"] <- X.Dev[2]
X.psi.pred.hi[, "Well_1km"] <- exp(mod$mcmcOutput$ALPHA0.Well_1km +
                                     mod$mcmcOutput$ALPHA.Dev_lo.Well_1km * X.DEV[1] +
                                     mod$mcmcOutput$ALPHA.Dev_bg.Well_1km * X.DEV[2]) %>%
  (function(x) (x - mean(X.psi.raw[,"Well_1km"])) / sd(X.psi.raw[,"Well_1km"]))
X.psi.pred.hi[, "Road_125m"] <- exp(mod$mcmcOutput$alpha0.Road_125m +
                                     mod$mcmcOutput$alpha.Dev_lo.Road_125m * X.DEV[1] +
                                     mod$mcmcOutput$alpha.Dev_bg.Road_125m * X.DEV[2]) %>%
  (function(x) (x - mean(X.psi.raw[,"Road_125m"])) / sd(X.psi.raw[,"Road_125m"]))
X.psi.pred.hi[, "AHerb"] <- QSLpersonal::expit(mod$mcmcOutput$alpha0.AHerb +
                                                 mod$mcmcOutput$alpha.Dev_lo.AHerb * X.DEV[1] +
                                                 mod$mcmcOutput$alpha.Dev_bg.AHerb * X.DEV[2] +
                                                 mod$mcmcOutput$alpha.Well_1km.AHerb * X.psi.pred.hi[, "Well_1km"] +
                                                 mod$mcmcOutput$alpha.Road_125m.AHerb * X.psi.pred.hi[, "Road_125m"]) %>%
  (function(x) (x - mean(X.psi.raw[,"AHerb"])) / sd(X.psi.raw[,"AHerb"]))

# Low development occupancy #
X.Dev <- X.psi[which(X.psi.raw[,"Dev_lo"] == 1), c("Dev_lo", "Dev_bg")][1,]
X.DEV <- X.PSI[which(X.PSI.raw[,"Dev_lo"] == 1), c("Dev_lo", "Dev_bg")][1,]
X.psi.pred.lo[, "Dev_lo"] <- X.Dev[1]
X.psi.pred.lo[, "Dev_bg"] <- X.Dev[2]
X.psi.pred.lo[, "Well_1km"] <- exp(mod$mcmcOutput$ALPHA0.Well_1km +
                                     mod$mcmcOutput$ALPHA.Dev_lo.Well_1km * X.DEV[1] +
                                     mod$mcmcOutput$ALPHA.Dev_bg.Well_1km * X.DEV[2]) %>%
  (function(x) (x - mean(X.psi.raw[,"Well_1km"])) / sd(X.psi.raw[,"Well_1km"]))
X.psi.pred.lo[, "Road_125m"] <- exp(mod$mcmcOutput$alpha0.Road_125m +
                                      mod$mcmcOutput$alpha.Dev_lo.Road_125m * X.DEV[1] +
                                      mod$mcmcOutput$alpha.Dev_bg.Road_125m * X.DEV[2]) %>%
  (function(x) (x - mean(X.psi.raw[,"Road_125m"])) / sd(X.psi.raw[,"Road_125m"]))
X.psi.pred.lo[, "AHerb"] <- QSLpersonal::expit(mod$mcmcOutput$alpha0.AHerb +
                                                 mod$mcmcOutput$alpha.Dev_lo.AHerb * X.DEV[1] +
                                                 mod$mcmcOutput$alpha.Dev_bg.AHerb * X.DEV[2] +
                                                 mod$mcmcOutput$alpha.Well_1km.AHerb * X.psi.pred.lo[, "Well_1km"] +
                                                 mod$mcmcOutput$alpha.Road_125m.AHerb * X.psi.pred.lo[, "Road_125m"]) %>%
  (function(x) (x - mean(X.psi.raw[,"AHerb"])) / sd(X.psi.raw[,"AHerb"]))

# Background occupancy #
X.Dev <- X.psi[which(X.psi.raw[,"Dev_bg"] == 1), c("Dev_lo", "Dev_bg")][1,]
X.DEV <- X.PSI[which(X.PSI.raw[,"Dev_bg"] == 1), c("Dev_lo", "Dev_bg")][1,]
X.psi.pred.bg[, "Dev_lo"] <- X.Dev[1]
X.psi.pred.bg[, "Dev_bg"] <- X.Dev[2]
X.psi.pred.bg[, "Well_1km"] <- exp(mod$mcmcOutput$ALPHA0.Well_1km +
                                     mod$mcmcOutput$ALPHA.Dev_lo.Well_1km * X.DEV[1] +
                                     mod$mcmcOutput$ALPHA.Dev_bg.Well_1km * X.DEV[2]) %>%
  (function(x) (x - mean(X.psi.raw[,"Well_1km"])) / sd(X.psi.raw[,"Well_1km"]))
X.psi.pred.bg[, "Road_125m"] <- exp(mod$mcmcOutput$alpha0.Road_125m +
                                      mod$mcmcOutput$alpha.Dev_lo.Road_125m * X.DEV[1] +
                                      mod$mcmcOutput$alpha.Dev_bg.Road_125m * X.DEV[2]) %>%
  (function(x) (x - mean(X.psi.raw[,"Road_125m"])) / sd(X.psi.raw[,"Road_125m"]))
X.psi.pred.bg[, "AHerb"] <- QSLpersonal::expit(mod$mcmcOutput$alpha0.AHerb +
                                                 mod$mcmcOutput$alpha.Dev_lo.AHerb * X.DEV[1] +
                                                 mod$mcmcOutput$alpha.Dev_bg.AHerb * X.DEV[2] +
                                                 mod$mcmcOutput$alpha.Well_1km.AHerb * X.psi.pred.bg[, "Well_1km"] +
                                                 mod$mcmcOutput$alpha.Road_125m.AHerb * X.psi.pred.bg[, "Road_125m"]) %>%
  (function(x) (x - mean(X.psi.raw[,"AHerb"])) / sd(X.psi.raw[,"AHerb"]))

# Trend #
X.lambda.pred.hi <- X.psi.pred.hi[, dimnames(X.lambda)[[2]]]
X.lambda.pred.lo <- X.psi.pred.lo[, dimnames(X.lambda)[[2]]]
X.lambda.pred.bg <- X.psi.pred.bg[, dimnames(X.lambda)[[2]]]

rm(X.Dev, X.DEV)