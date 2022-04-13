model <<- nimbleCode({
  ########################
  # Bird community model #
  ########################

  for (i in 1:n.spp) {
    for(j in 1:n.grdyr) {
      # Observation process (total detected)
      Y.sum[j, i] ~ dbin(p_star[j, i], z[j, i])
      p_star[j, i] <- 1 - pow(1 - p[j, i], K)
      logit(p[j, i]) <- zeta0[i] + dev.zeta[i, yearID.grdyr[j]] +
        inprod(zeta1[i, 1:p.zeta], X.zeta[j, 1:p.zeta])
      
      # Point-level ecological process
      z[j, i] ~ dbin(psi[j, i] * Z[j, i], n.point[j])
      logit(psi[j, i]) <- beta0[i] + dev.beta[i, yearID.grdyr[j]] +
        inprod(beta1[i, 1:p.psi], X.psi[j, 1:p.psi]) +
        lambda[j, i] * X.trend[yearID.grdyr[j]]
      lambda[j, i] <- delta0[i] +
        inprod(delta1[i, 1:p.lambda], X.lambda[j, 1:p.lambda])
        
      # Grid-level ecological process
      Z[j, i] ~ dbern(PSI[j, i] * w[i])
      logit(PSI[j, i]) <- BETA0[i] + dev.BETA[i, yearID.grdyr[j]] +
        inprod(BETA1[i, 1:p.PSI], X.PSI[j, 1:p.PSI]) +
        LAMBDA[j, i] * X.trend[yearID.grdyr[j]]
      LAMBDA[j, i] <- DELTA0[i] +
        inprod(DELTA1[i, 1:p.LAMBDA], X.LAMBDA[j, 1:p.LAMBDA])
      }
      
    w[i] ~ dbern(omega) # Species inclusion in community

    # Species-specific intercepts
    BETA0[i] ~ dnorm(BETA0.mu, pow(sigma.BETA0, -2))
    DELTA0[i] ~ dnorm(DELTA0.mu, pow(sigma.DELTA0, -2))
    beta0[i] ~ dnorm(beta0.mu + (rho.bB*sigma.beta0 / sigma.BETA0) * (BETA0[i] - BETA0.mu),
      pow(sigma.beta0, -2) / (1 - pow(rho.bB, 2)))
    delta0[i] ~ dnorm(delta0.mu, pow(sigma.delta0, -2))
    zeta0[i] ~ dnorm(zeta0.mu + (rho.zb*sigma.zeta0 / sigma.beta0) * (beta0[i] - beta0.mu),
      pow(sigma.zeta0, -2)/(1 - pow(rho.zb, 2)))

    # Species-specific yearly deviations
    for(t in 1:n.year) {
      dev.zeta[i, t] ~ dnorm(0, pow(sigma.z0, -2)) # Year-specific deviations in detectability from mean for each species
      dev.beta[i, t] ~ dnorm(0, pow(sigma.b0, -2)) # Year-specific deviations in point occupancy from mean for each species
      dev.BETA[i, t] ~ dnorm(0, pow(sigma.B0, -2)) # Year-specific deviations in grid occupancy from mean for each species
    }
    
    # Species-specific covariate relationships
    for(b in 1:p.zeta) {
      zeta1[i, b] ~ dnorm(zeta1.mu[b], pow(sigma.zeta1[b], -2))
    }

    for(b in 1:p.psi) {
      beta1[i, b] ~ dnorm(beta1.mu[b], pow(sigma.beta1[b], -2))
    }

    for(b in 1:p.lambda) {
      delta1[i, b] ~ dnorm(delta1.mu[b], pow(sigma.delta1[b], -2))
    }

    for(b in 1:p.PSI) {
      BETA1[i, b] ~ dnorm(BETA1.mu[b], pow(sigma.BETA1[b], -2))
    }

    for(b in 1:p.LAMBDA) {
      DELTA1[i, b] ~ dnorm(DELTA1.mu[b], pow(sigma.DELTA1[b], -2))
    }
  }
  
  # Observation process (detections within minute intervals)
  for(i in 1:n.int) {
    Y.int[i, 1:K] ~ dmulti(pclass[i, 1:K], Y.sum[grdyrID.int[i], sppID.int[i]])
    for(k in 1:K) {
      pint[i, k] <- pow(1 - p[grdyrID.int[i], sppID.int[i]], (k - 1)) * p[grdyrID.int[i], sppID.int[i]]
      pclass[i, k] <- pint[i, k] / sum(pint[i, 1:K])
    }
  }

  ### Prior distributions ###
  
  # parameter correlations
  
  rho.zb ~ dunif(-1, 1) # Correlation between detectability and point occupancy
  rho.bB ~ dunif(-1, 1) # Correlation between point and grid occupancy

  zeta0.mu ~  dt(0, pow(t.sigma, -2), t.nu) # logit detection probability intercept
  tvar.sigma.zeta0 ~ dt(0,1,1)  # Cauchy distribution
  sigma.zeta0 <- abs(tvar.sigma.zeta0)  # half-Cauchy distribution
    # Year effect detectability hyper-parameters #
  tvar.sigma.z0 ~ dt(0,1,1)  # Cauchy distribution
  sigma.z0 <- abs(tvar.sigma.z0)  # half-Cauchy distribution

  beta0.mu ~  dt(0, pow(t.sigma, -2), t.nu) # logit point occupancy intercept
  sigma.beta0 <- abs(tvar.sigma.beta0)  # half-Cauchy distribution
  tvar.sigma.beta0 ~ dt(0,1,1)  # Cauchy distribution
    # Year effect point occupancy hyper-parameters #
  tvar.sigma.b0 ~ dt(0,1,1)  # Cauchy distribution
  sigma.b0 <- abs(tvar.sigma.b0)  # half-Cauchy distribution

  delta0.mu ~  dt(0, pow(t.sigma, -2), t.nu) # logit point occupancy intercept
  sigma.delta0 <- abs(tvar.sigma.delta0)  # half-Cauchy distribution
  tvar.sigma.delta0 ~ dt(0,1,1)  # Cauchy distribution

  BETA0.mu ~  dt(0, pow(t.sigma, -2), t.nu) # logit point occupancy intercept
  sigma.BETA0 <- abs(tvar.sigma.BETA0)  # half-Cauchy distribution
  tvar.sigma.BETA0 ~ dt(0,1,1)  # Cauchy distribution
    # Year effect point occupancy hyper-parameters #
  tvar.sigma.B0 ~ dt(0,1,1)  # Cauchy distribution
  sigma.B0 <- abs(tvar.sigma.B0)  # half-Cauchy distribution

  DELTA0.mu ~  dt(0, pow(t.sigma, -2), t.nu) # logit point occupancy intercept
  sigma.DELTA0 <- abs(tvar.sigma.DELTA0)  # half-Cauchy distribution
  tvar.sigma.DELTA0 ~ dt(0,1,1)  # Cauchy distribution

  # Probability of species inclusion in the metacommunity
  omega ~ dunif(0,1)
  
  # covariates for p (zeta)
  for(b in 1:p.zeta) {
    zeta1.mu[b] ~ dt(0, pow(t.sigma, -2), t.nu)
    sigma.zeta1[b] <- abs(tvar.zeta1[b]) # half-Cauchy distribution
    tvar.zeta1[b] ~ dt(0,1,1)  # Cauchy distribution
  }

  # covariates for psi (beta)
  for(b in 1:p.psi) {
    beta1.mu[b] ~ dt(0, pow(t.sigma, -2), t.nu)
    sigma.beta1[b] <- abs(tvar.beta1[b]) # half-Cauchy distribution
    tvar.beta1[b] ~ dt(0,1,1)  # Cauchy distribution
  }

  # covariates for lambda (delta)
  for(b in 1:p.lambda) {
    delta1.mu[b] ~ dt(0, pow(t.sigma, -2), t.nu)
    sigma.delta1[b] <- abs(tvar.delta1[b]) # half-Cauchy distribution
    tvar.delta1[b] ~ dt(0,1,1)  # Cauchy distribution
  }

  # covariates for PSI (BETA)
  for(b in 1:p.PSI) {
    BETA1.mu[b] ~ dt(0, pow(t.sigma, -2), t.nu)
    sigma.BETA1[b] <- abs(tvar.BETA1[b]) # half-Cauchy distribution
    tvar.BETA1[b] ~ dt(0,1,1)  # Cauchy distribution
  }

  # covariates for LAMBDA (DELTA)
  for(b in 1:p.LAMBDA) {
    DELTA1.mu[b] ~ dt(0, pow(t.sigma, -2), t.nu)
    sigma.DELTA1[b] <- abs(tvar.DELTA1[b]) # half-Cauchy distribution
    tvar.DELTA1[b] ~ dt(0,1,1)  # Cauchy distribution
  }

  t.nu <- 7.763179      # Uniform prior
  t.sigma <- 1.566267   # Uniform prior
  
  #####################################
  # Intermediate mechanistic pathways #
  #####################################
  
  for(j in 1:n.grdyr) {
    # Grid cell scale pathways #
    X.PSI.raw[j, ind.PSI.Well_3km] ~ dnegbin(pred.Well_3km[j], r.Well_3km)
    pred.Well_3km[j] <- r.Well_3km / (r.Well_3km + lambda.Well_3km[j])
    log(lambda.Well_3km[j]) <- ALPHA0.Well_3km +
      ALPHA.Dev_lo.Well_3km * X.PSI[j, ind.PSI.Dev_lo] +
      ALPHA.Dev_bg.Well_3km * X.PSI[j, ind.PSI.Dev_bg]

    X.PSI.raw[j, ind.PSI.Well_1km] ~ dnegbin(pred.Well_1km[j], r.Well_1km)
    pred.Well_1km[j] <- r.Well_1km / (r.Well_1km + lambda.Well_1km[j])
    log(lambda.Well_1km[j]) <- ALPHA0.Well_1km +
      ALPHA.Dev_lo.Well_1km * X.PSI[j, ind.PSI.Dev_lo] +
      ALPHA.Dev_bg.Well_1km * X.PSI[j, ind.PSI.Dev_bg]

    X.PSI.raw[j, ind.PSI.Road_1km] ~ dgamma(shape.Road_1km,
      shape.Road_1km / exp(pred.Road_1km[j]))
    pred.Road_1km[j] <- ALPHA0.Road_1km +
      ALPHA.Dev_lo.Road_1km * X.PSI[j, ind.PSI.Dev_lo] +
      ALPHA.Dev_bg.Road_1km * X.PSI[j, ind.PSI.Dev_bg]

    # Point scale pathways #
    X.psi.raw[j, ind.psi.AHerb] ~ dbeta(a.AHerb[j], b.AHerb[j])
    a.AHerb[j] <- pred.AHerb[j] * phi.AHerb
    b.AHerb[j] <- (1 - pred.AHerb[j]) * phi.AHerb
    logit(pred.AHerb[j]) <- alpha0.AHerb +
      alpha.Well_1km.AHerb * X.psi[j, ind.psi.Well_1km] +
      alpha.Road_1km.AHerb * X.psi[j, ind.psi.Road_125m]

    X.psi.raw[j, ind.psi.Road_125m] ~ dgamma(shape.Road_125m,
      shape.Road_125m / exp(pred.Road_125m[j]))
    pred.Road_125m[j] <- alpha0.Road_125m +
      alpha.Dev_lo.Road_125m * X.PSI[j, ind.PSI.Dev_lo] +
      alpha.Dev_bg.Road_125m * X.PSI[j, ind.PSI.Dev_bg]
  }

  ### Prior distributions ###
  ## Grid cell level models ##
  ALPHA0.Well_3km ~ dnorm(0, 0.66667)
  ALPHA.Dev_lo.Well_3km ~ dnorm(0, 0.66667)
  ALPHA.Dev_bg.Well_3km ~ dnorm(0, 0.66667)
  r.Well_3km ~ dunif(0, 50)

  ALPHA0.Well_1km ~ dnorm(0, 0.66667)
  ALPHA.Dev_lo.Well_1km ~ dnorm(0, 0.66667)
  ALPHA.Dev_bg.Well_1km ~ dnorm(0, 0.66667)
  r.Well_1km ~ dunif(0, 50)

  ALPHA0.Road_1km ~ dnorm(0, 0.5)
  ALPHA.Dev_lo.Road_1km ~ dnorm(0, 0.5)
  ALPHA.Dev_bg.Road_1km ~ dnorm(0, 0.5)
  shape.Road_1km ~ dunif(0, 100)

  ## Point level models ##
  alpha0.AHerb ~ dnorm(0, 0.66667)
  alpha.Well_1km.AHerb ~ dnorm(0, 0.66667)
  alpha.Road_1km.AHerb ~ dnorm(0, 0.66667)
  phi.AHerb ~ dgamma(.1,.1)
  
  alpha0.Road_125m ~ dnorm(0, 1)
  alpha.Dev_lo.Road_125m ~ dnorm(0, 1)
  alpha.Dev_bg.Road_125m ~ dnorm(0, 1)
  shape.Road_125m ~ dunif(0, 100)
})
