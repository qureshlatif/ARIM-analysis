model <<- nimbleCode({
  ########################
  # Bird community model #
  ########################

  for (i in 1:n.spp) {
    for(j in 1:n.pntyr) {
      # Observation process
      Y[j, i] ~ dbinom(prob.y[j, i], TPeriod[j, i])
      prob.y[j, i] <- (Y.point[j, i] * Y.grid[gridXyrID[j], i] * Y.spp[i] * prob.y1Y1YS1[j, i]) +
        ((1 - Y.point[j, i]) * Y.grid[gridXyrID[j], i] * Y.spp[i] * prob.y0Y1YS1[j, i]) +
        ((1 - Y.grid[gridXyrID[j], i]) * Y.spp[i] * prob.y0Y0YS1[j, i]) +
        ((1 - Y.spp[i]) * prob.y0Y0YS0[j, i])
      prob.y1Y1YS1[j, i] <- p[j, i] * psi[j, i] * PSI[gridID[j], yearID[j], i] * omega
      prob.y0Y1YS1[j, i] <- pow(1 - p[j, i], TPeriod[j, i]) * psi[j, i] *
        PSI[gridID[j], yearID[j], i] * omega +
        (1 - psi[j, i]) * PSI[gridID[j], yearID[j], i] * omega
      prob.y0Y0YS1[j, i] <- pow(1 - p[j, i], TPeriod[j, i]) * psi[j, i] *
        PSI[gridID[j], yearID[j], i] * omega +
        (1 - psi[j, i]) * PSI[gridID[j], yearID[j], i] * omega +
        (1 - PSI[gridID[j], yearID[j], i]) * omega
      prob.y0Y0YS0[j, i] <- pow(1 - p[j, i], TPeriod[j, i]) * psi[j, i] *
        PSI[gridID[j], yearID[j], i] * omega +
        (1 - psi[j, i]) * PSI[gridID[j], yearID[j], i] * omega +
        (1 - PSI[gridID[j], yearID[j], i]) * omega +
        (1 - omega)
        
      logit(p[j, i]) <- zeta0[i] + dev.zeta[i, yearID[j]] +
        inprod(zeta1[i, 1:p.zeta], X.zeta[j, 1:p.zeta])
      
      # Point-level ecological process
      logit(psi[j, i]) <- beta0[i] + dev.beta[i, yearID[j]] +
        inprod(beta1[i, 1:p.psi.init], X.psi[j, 1:p.psi.init]) +
        lambda[j, i] * (yearID[j] - 1)
      lambda[j, i] <- delta0[i] +
        inprod(delta1[i, 1:p.psi.dyn], X.psi.dyn[j, 1:p.psi.dyn])
    }

    # Grid-level ecological process
    for(k in 1:n.grid) {
      for(t in 1:n.year) {
        logit(PSI[k, t, i]) <- BETA0[i] + dev.BETA[i, t] +
          inprod(BETA1[i, 1:p.PSI], X.PSI[k, 1, 1:p.PSI]) +
          LAMBDA[k, t, i] * (t - 1)
        LAMBDA[k, t, i] <- DELTA0[i] +
          inprod(DELTA1[i, 1:p.PSI], X.PSI[k, t, 1:p.PSI])
      }
    }

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

    for(b in 1:p.psi.init) {
      beta1[i, b] <- beta1.base[i, b] +
        inprod(beta1.offset[1:(n.guild-1), b], guildMem[i, 2:n.guild])
      beta1.base[i, b] ~ dnorm(beta1.mu[b], pow(sigma.beta1[b], -2))
    }

    for(b in 1:p.psi.dyn) {
      delta1[i, b] <- delta1.base[i, b] +
        inprod(delta1.offset[1:(n.guild-1), b], guildMem[i, 2:n.guild])
      delta1.base[i, b] ~ dnorm(delta1.mu[b], pow(sigma.delta1[b], -2))
    }

    for(b in 1:p.PSI) {
      BETA1[i, b] <- BETA1.base[i, b] +
        inprod(BETA1.offset[1:(n.guild-1), b], guildMem[i, 2:n.guild])
      BETA1.base[i, b] ~ dnorm(BETA1.mu[b], pow(sigma.BETA1[b], -2))
    }

    for(b in 1:p.PSI) {
      DELTA1[i, b] <- DELTA1.base[i, b] +
        inprod(DELTA1.offset[1:(n.guild-1), b], guildMem[i, 2:n.guild])
      DELTA1.base[i, b] ~ dnorm(DELTA1.mu[b], pow(sigma.DELTA1[b], -2))
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
  for(b in 1:p.psi.init) {
    beta1.mu[b] ~ dt(0, pow(t.sigma, -2), t.nu)
    sigma.beta1[b] <- abs(tvar.beta1[b]) # half-Cauchy distribution
    tvar.beta1[b] ~ dt(0,1,1)  # Cauchy distribution
  }
    # Guild-level offsets
  for(ib in 1:length(ind.psi.init.offset)) {
    for(g in 1:(n.guild-1)) {
      beta1.offset[g, ind.psi.init.offset[ib]] ~ dt(0, pow(t.sigma, -2), t.nu)
    }
  }
  for(ib in 1:length(ind.psi.init.no_offset)) {
    for(g in 1:(n.guild-1)) {
      beta1.offset[g, ind.psi.init.no_offset[ib]] <- 0
    }
  }

  # covariates for lambda (delta)
  for(b in 1:p.psi.dyn) {
    delta1.mu[b] ~ dt(0, pow(t.sigma, -2), t.nu)
    sigma.delta1[b] <- abs(tvar.delta1[b]) # half-Cauchy distribution
    tvar.delta1[b] ~ dt(0,1,1)  # Cauchy distribution
  }
    # Guild-level offsets
  for(ib in 1:length(ind.psi.dyn.offset)) {
    for(g in 1:(n.guild-1)) {
      delta1.offset[g, ind.psi.dyn.offset[ib]] ~ dt(0, pow(t.sigma, -2), t.nu)
    }
  }
  for(ib in 1:length(ind.psi.dyn.no_offset)) {
    for(g in 1:(n.guild-1)) {
      delta1.offset[g, ind.psi.dyn.no_offset[ib]] <- 0
    }
  }

  # covariates for PSI (BETA)
  for(b in 1:p.PSI) {
    BETA1.mu[b] ~ dt(0, pow(t.sigma, -2), t.nu)
    sigma.BETA1[b] <- abs(tvar.BETA1[b]) # half-Cauchy distribution
    tvar.BETA1[b] ~ dt(0,1,1)  # Cauchy distribution
  }
    # Guild-level offsets
  for(ib in 1:length(ind.PSI.offset)) {
    for(g in 1:(n.guild-1)) {
      BETA1.offset[g, ind.PSI.offset[ib]] ~ dt(0, pow(t.sigma, -2), t.nu)
    }
  }
  for(ib in 1:length(ind.PSI.no_offset)) {
    for(g in 1:(n.guild-1)) {
      BETA1.offset[g, ind.PSI.no_offset[ib]] <- 0
    }
  }

  # covariates for LAMBDA (DELTA)
  for(b in 1:p.PSI) {
    DELTA1.mu[b] ~ dt(0, pow(t.sigma, -2), t.nu)
    sigma.DELTA1[b] <- abs(tvar.DELTA1[b]) # half-Cauchy distribution
    tvar.DELTA1[b] ~ dt(0,1,1)  # Cauchy distribution
  }
    # Guild-level offsets
  for(ib in 1:length(ind.PSI.offset)) {
    for(g in 1:(n.guild-1)) {
      DELTA1.offset[g, ind.PSI.offset[ib]] ~ dt(0, pow(t.sigma, -2), t.nu)
    }
  }
  for(ib in 1:length(ind.PSI.no_offset)) {
    for(g in 1:(n.guild-1)) {
      DELTA1.offset[g, ind.PSI.no_offset[ib]] <- 0
    }
  }

  t.nu <- 7.763179      # Uniform prior
  t.sigma <- 1.566267   # Uniform prior
  
  #####################################
  # Intermediate mechanistic pathways #
  #####################################
  
  # Grid cell scale pathways #
  for(k in 1:n.grid) {
    for(t in 1:n.year) {
      X.PSI.raw[k, t, ind.Well_3km] ~ dnegbin(pred.Well_3km[k, t], r.Well_3km)
      pred.Well_3km[k, t] <- r.Well_3km / (r.Well_3km + lambda.Well_3km[k, t])
      log(lambda.Well_3km[k, t]) <- ALPHA0.Well_3km +
        ALPHA.Dev_lo.Well_3km * X.PSI[k, t, ind.PSI.Dev_lo] +
        ALPHA.Dev_bg.Well_3km * X.PSI[k, t, ind.PSI.Dev_bg]

      X.PSI.raw[k, t, ind.Well_1km] ~ dnegbin(pred.Well_1km[k, t], r.Well_1km)
      pred.Well_1km[k, t] <- r.Well_1km / (r.Well_1km + lambda.Well_1km[k, t])
      log(lambda.Well_1km[k, t]) <- ALPHA0.Well_1km +
        ALPHA.Dev_lo.Well_1km * X.PSI[k, t, ind.PSI.Dev_lo] +
        ALPHA.Dev_bg.Well_1km * X.PSI[k, t, ind.PSI.Dev_bg]
    }
    
    for(t in 1:2) { # Years 2-n.year are redundant.
      X.PSI.raw[k, t, ind.Road_1km] ~ dgamma(shape.Road_1km,
        shape.Road_1km / exp(pred.Road_1km[k, t]))
      pred.Road_1km[k, t] <- ALPHA0.Road_1km +
        ALPHA.Dev_lo.Road_1km * X.PSI[k, t, ind.PSI.Dev_lo] +
        ALPHA.Dev_bg.Road_1km * X.PSI[k, t, ind.PSI.Dev_bg]
    }
  }
    
  # Point scale pathways #
  for(j in 1:n.pntyr) {
    X.psi.raw[j, ind.Well_125m] ~ dbern(pred.Well_125m[j])
    logit(pred.Well_125m[j]) <- alpha0.Well_125m +
      alpha.Dev_lo.Well_125m * X.PSI[gridID[j], yearID[j], ind.PSI.Dev_lo] +
      alpha.Dev_bg.Well_125m * X.PSI[gridID[j], yearID[j], ind.PSI.Dev_bg]

    X.psi.raw[j, ind.AHerb] ~ dbeta(a.AHerb[j], b.AHerb[j])
    a.AHerb[j] <- pred.AHerb[j] * phi.AHerb
    b.AHerb[j] <- (1 - pred.AHerb[j]) * phi.AHerb
    logit(pred.AHerb[j]) <- alpha0.AHerb +
      alpha.WellA_125m.AHerb * X.psi[j, ind.Well_125m] +
      alpha.Road_125m.AHerb * X.psi[j, ind.Road_125m]

    X.psi.raw[j, ind.D_Road_125m] ~ dgamma(shape.D_Road_125m,
      shape.D_Road_125m / exp(pred.D_Road_125m[j])) # Need to compile in R. Should be NA where roads are absent.
    X.psi.raw[j, ind.PA_Road_125m] ~ dbern(pred.PA_Road_125m[j]) # Need to compile this in R.
    logit(pred.PA_Road_125m[j]) <- alpha0.PA_Road_125m +
      alpha.Dev_lo.PA_Road_125m * X.PSI[gridID[j], yearID[j], ind.PSI.Dev_lo] +
      alpha.Dev_bg.PA_Road_125m * X.PSI[gridID[j], yearID[j], ind.PSI.Dev_bg]
    pred.D_Road_125m[j] <- alpha0.D_Road_125m +
      alpha.Dev_lo.D_Road_125m * X.PSI[gridID[j], yearID[j], ind.PSI.Dev_lo] +
      alpha.Dev_bg.D_Road_125m * X.PSI[gridID[j], yearID[j], ind.PSI.Dev_bg]
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
  alpha0.Well_125m ~ dnorm(0, 0.66667)
  alpha.Dev_lo.Well_125m ~ dnorm(0, 0.66667)
  alpha.Dev_bg.Well_125m ~ dnorm(0, 0.66667)

  alpha0.PA_Road_125m ~ dnorm(0, 0.66667)
  alpha.Dev_lo.PA_Road_125m ~ dnorm(0, 0.66667)
  alpha.Dev_bg.PA_Road_125m ~ dnorm(0, 0.66667)
  alpha0.D_Road_125m ~ dnorm(0, 1)
  alpha.Dev_lo.D_Road_125m ~ dnorm(0, 1)
  alpha.Dev_bg.D_Road_125m ~ dnorm(0, 1)
  shape.D_Road_125m ~ dunif(0, 100)

  alpha0.AHerb ~ dnorm(0, 0.66667)
  alpha.WellA_125m.AHerb ~ dnorm(0, 0.66667)
  alpha.Road_125m.AHerb ~ dnorm(0, 0.66667)
  phi.AHerb ~ dgamma(.1,.1)
})
