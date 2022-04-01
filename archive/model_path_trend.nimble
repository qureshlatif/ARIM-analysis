model <<- nimbleCode({
  ########################
  # Bird community model #
  ########################

  for (i in 1:n.spp) {
    # Observation process
    for(j in 1:n.pntyr) {
      Y[j, i] ~ dbinom(p[j, i]*z[pointID[j], yearID[j], i], TPeriod[j, i])
      logit(p[j, i]) <- zeta0[i] + dev.zeta[i, yearID[j]] +
        inprod(zeta1[i, 1:p.zeta], X.zeta[j, 1:p.zeta])
    }

    # Ecological state processes
    for(j in 1:n.point) {
      for(t in 1:n.year) {
        z[j, t, i] ~ dbern(psi[j, t, i] * Z[gridID[j], t, i])
        logit(psi[j, t, i]) <- beta0[i] + dev.beta[i, t] +
          inprod(beta1[i, 1:p.psi.init], X.psi[j, 1, 1:p.psi.init]) +
          lambda[j, t, i] * (t - 1)
        lambda[j, t, i] <- delta0[i] +
          inprod(delta1[i, 1:p.psi.dyn], X.psi[j, t, ind.psi.dyn])
          
        #___ For deriving yearly occupancy rates ___#
        z_bkg[j, t, i] <-
          z[j, t, i] * X.PSI[gridID[j], 1, ind.PSI.Dev_bg]
        z_dlo[j, t, i] <-
          z[j, t, i] * X.PSI[gridID[j], 1, ind.PSI.Dev_lo]
        z_dhi[j, t, i] <-
          z[j, t, i] * (1 - (X.PSI[gridID[j], 1, ind.PSI.Dev_bg] +
          X.PSI[gridID[j], 1, ind.PSI.Dev_lo]))
        #___________________________________________#
      }
    }
    
    for(k in 1:n.grid) {
      for(t in 1:n.year) {
        Z[k, t, i] ~ dbern(PSI[k, t, i] * w[i])
        logit(PSI[k, t, i]) <- BETA0[i] + dev.BETA[i, t] +
          inprod(BETA1[i, 1:p.PSI], X.PSI[k, 1, 1:p.PSI]) +
          LAMBDA[k, t, i] * (t - 1)
        LAMBDA[k, t, i] <- DELTA0[i] +
          inprod(DELTA1[i, 1:p.PSI], X.PSI[k, t, 1:p.PSI])
        
        #___ For deriving yearly occupancy rates ___#
        Z_bkg[k, t, i] <-
          Z[k, t, i] * X.PSI[k, 1, ind.PSI.Dev_bg]
        Z_dlo[k, t, i] <-
          Z[k, t, i] * X.PSI[k, 1, ind.PSI.Dev_lo]
        Z_dhi[k, t, i] <-
          Z[k, t, i] * (1 - (X.PSI[k, 1, ind.PSI.Dev_bg] +
            X.PSI[k, 1, ind.PSI.Dev_lo]))
        #___________________________________________#
      }
    }
    
    #___ Derived yearly finite-sample occupancy by development level___#
    for(t in 1:n.year) {
      occ.pnt.yr[t, i, 1] <- sum(z_bkg[1:n.point, t, i]) / n.pnt.bkg
      occ.pnt.yr[t, i, 2] <- sum(z_dlo[1:n.point, t, i]) / n.pnt.dlo
      occ.pnt.yr[t, i, 3] <- sum(z_dhi[1:n.point, t, i]) / n.pnt.dhi
      occ.grd.yr[t, i, 1] <- sum(Z_bkg[1:n.grid, t, i]) / n.grd.bkg
      occ.grd.yr[t, i, 2] <- sum(Z_dlo[1:n.grid, t, i]) / n.grd.dlo
      occ.grd.yr[t, i, 3] <- sum(Z_dhi[1:n.grid, t, i]) / n.grd.dhi
    }
    #__________________________________________________________________#

    # Species-specific intercepts
    B[i, 1:3] ~ dmnorm(B.mu, invSigma.B) # Multivariate Wishart approach for initial occupancy and detectability #
    BETA0[i] <- B[i, 1]
    DELTA0[i] ~ dnorm(DELTA0.mu, pow(sigma.DELTA0, -2))
    beta0[i] <- B[i, 2]
    delta0[i] ~ dnorm(delta0.mu, pow(sigma.delta0, -2))
    zeta0[i] <- B[i, 3]

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

    w[i] ~ dbern(omega)
  }

  ### Prior distributions ###
  
  # parameter correlations
  
  rho.zb ~ dunif(-1,1) # Correlation between detectability and point occupancy
  rho.bB ~ dunif(-1,1) # Correlation between point and grid occupancy
  rho.zB ~ dunif(-1,1) # Correlation between detectability and grid occupancy
  
  # mean and precison for the parameter intercepts
  
  #____ Multivariate Wishart prior for variance-covariance of occupancy and detectability ____#
  R[1,1] <- pow(sigma.BETA0, 2)
  R[2,2] <- pow(sigma.beta0, 2)
  R[3,3] <- pow(sigma.zeta0, 2)
  R[1,2] <- rho.bB*sigma.BETA0*sigma.beta0
  R[2,1] <- rho.bB*sigma.BETA0*sigma.beta0
  R[2,3] <- rho.zb*sigma.beta0*sigma.zeta0
  R[3,2] <- rho.zb*sigma.beta0*sigma.zeta0
  R[1,3] <- rho.zB*sigma.BETA0*sigma.zeta0
  R[3,1] <- rho.zB*sigma.BETA0*sigma.zeta0
  invSigma.B ~ dwish(R, 3)
  Sigma.B <- inverse(invSigma.B) # Not sure this is needed except possibly to be saved.

  B.mu[1] <- BETA0.mu
  B.mu[2] <- beta0.mu
  B.mu[3] <- zeta0.mu
  #___________________________________________________________________________________________#

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
  for(b in ind.psi.init.offset) {
    for(g in 1:(n.guild-1)) {
      beta1.offset[g, b] ~ dt(0, pow(t.sigma, -2), t.nu)
    }
  }
  for(b in ind.psi.init.no_offset) {
    for(g in 1:(n.guild-1)) {
      beta1.offset[g, b] <- 0
    }
  }

  # covariates for gamma (delta)
  for(b in 1:p.psi.dyn) {
    delta1.mu[b] ~ dt(0, pow(t.sigma, -2), t.nu)
    sigma.delta1[b] <- abs(tvar.delta1[b]) # half-Cauchy distribution
    tvar.delta1[b] ~ dt(0,1,1)  # Cauchy distribution
  }
    # Guild-level offsets
  for(b in ind.psi.dyn.offset) {
    for(g in 1:(n.guild-1)) {
      delta1.offset[g, b] ~ dt(0, pow(t.sigma, -2), t.nu)
    }
  }
  for(b in ind.psi.dyn.no_offset) {
    for(g in 1:(n.guild-1)) {
      delta1.offset[g, b] <- 0
    }
  }

  # covariates for PSI (BETA)
  for(b in 1:p.PSI) {
    BETA1.mu[b] ~ dt(0, pow(t.sigma, -2), t.nu)
    sigma.BETA1[b] <- abs(tvar.BETA1[b]) # half-Cauchy distribution
    tvar.BETA1[b] ~ dt(0,1,1)  # Cauchy distribution
  }
    # Guild-level offsets
  for(b in ind.PSI.offset) {
    for(g in 1:(n.guild-1)) {
      BETA1.offset[g, b] ~ dt(0, pow(t.sigma, -2), t.nu)
    }
  }
  for(b in ind.PSI.no_offset) {
    for(g in 1:(n.guild-1)) {
      BETA1.offset[g, b] <- 0
    }
  }

  # covariates for GAMMA (DELTA)
  for(b in 1:p.PSI) {
    DELTA1.mu[b] ~ dt(0, pow(t.sigma, -2), t.nu)
    sigma.DELTA1[b] <- abs(tvar.DELTA1[b]) # half-Cauchy distribution
    tvar.DELTA1[b] ~ dt(0,1,1)  # Cauchy distribution
  }
    # Guild-level offsets
  for(b in ind.PSI.offset) {
    for(g in 1:(n.guild-1)) {
      DELTA1.offset[g, b] ~ dt(0, pow(t.sigma, -2), t.nu)
    }
  }
  for(b in ind.PSI.no_offset) {
    for(g in 1:(n.guild-1)) {
      DELTA1.offset[g, b] <- 0
    }
  }

  t.nu <- 7.763179      # Uniform prior
  t.sigma <- 1.566267   # Uniform prior
  
  ################################################
  # Derived yearly richness by development level #
  ################################################
  
  for(t in 1:n.year) {
    for(j in 1:n.point) {
      SR.point[j, t, 1] <- sum(z_bkg[j, t, 1:n.spp])
      SR.point[j, t, 2] <- sum(z_dlo[j, t, 1:n.spp])
      SR.point[j, t, 3] <- sum(z_dhi[j, t, 1:n.spp])
      for(g in 1:n.guild) {
        SR.pnt.guild[j, t, g, 1] <- sum(z_bkg[j, t, 1:n.spp] * guildMem[1:n.spp, g])
        SR.pnt.guild[j, t, g, 2] <- sum(z_dlo[j, t, 1:n.spp] * guildMem[1:n.spp, g])
        SR.pnt.guild[j, t, g, 3] <- sum(z_dhi[j, t, 1:n.spp] * guildMem[1:n.spp, g])
      }
    }
    SR.pnt.yr[t, 1] <- sum(SR.point[1:n.point, t, 1]) / n.pnt.bkg
    SR.pnt.yr[t, 2] <- sum(SR.point[1:n.point, t, 2]) / n.pnt.dlo
    SR.pnt.yr[t, 3] <- sum(SR.point[1:n.point, t, 3]) / n.pnt.dhi
    for(g in 1:n.guild) {
      SR.pnt.gld.yr[t, g, 1] <- sum(SR.pnt.guild[1:n.point, t, g, 1]) / n.pnt.bkg
      SR.pnt.gld.yr[t, g, 2] <- sum(SR.pnt.guild[1:n.point, t, g, 2]) / n.pnt.dlo
      SR.pnt.gld.yr[t, g, 3] <- sum(SR.pnt.guild[1:n.point, t, g, 3]) / n.pnt.dhi
    }
    
    for(k in 1:n.grid) {
      SR.grid[k, t, 1] <- sum(Z_bkg[k, t, 1:n.spp])
      SR.grid[k, t, 2] <- sum(Z_dlo[k, t, 1:n.spp])
      SR.grid[k, t, 3] <- sum(Z_dhi[k, t, 1:n.spp])
      for(g in 1:n.guild) {
        SR.grd.guild[k, t, g, 1] <- sum(z_bkg[k, t, 1:n.spp] * guildMem[1:n.spp, g])
        SR.grd.guild[k, t, g, 2] <- sum(z_dlo[k, t, 1:n.spp] * guildMem[1:n.spp, g])
        SR.grd.guild[k, t, g, 3] <- sum(z_dhi[k, t, 1:n.spp] * guildMem[1:n.spp, g])
      }
    }
    SR.grd.yr[t, 1] <- sum(SR.grid[1:n.grid, t, 1]) / n.grd.bkg
    SR.grd.yr[t, 2] <- sum(SR.grid[1:n.grid, t, 2]) / n.grd.dlo
    SR.grd.yr[t, 3] <- sum(SR.grid[1:n.grid, t, 3]) / n.grd.dhi
    for(g in 1:n.guild) {
      SR.grd.gld.yr[t, g, 1] <- sum(SR.grd.guild[1:n.grid, t, g, 1]) / n.grd.bkg
      SR.grd.gld.yr[t, g, 2] <- sum(SR.grd.guild[1:n.grid, t, g, 2]) / n.grd.dlo
      SR.grd.gld.yr[t, g, 3] <- sum(SR.grd.guild[1:n.grid, t, g, 3]) / n.grd.dhi
    }
  }
  
  #____ Sample sizes for derived occupancy estimates ____#
  n.pnt.bkg <- sum(X.PSI[gridID, 1, ind.PSI.Dev_bg])
  n.pnt.dlo <- sum(X.PSI[gridID, 1, ind.PSI.Dev_lo])
  n.pnt.dhi <- (n.point - sum(X.PSI[gridID, 1, ind.PSI.Dev_bg]) -
    sum(sum(X.PSI[gridID, 1, ind.PSI.Dev_lo])))
  n.grd.bkg <- sum(X.PSI[1:n.grid, 1, ind.PSI.Dev_bg])
  n.grd.dlo <- sum(X.PSI[1:n.grid, 1, ind.PSI.Dev_lo])
  n.grd.dhi <- (n.grid - sum(X.PSI[1:n.grid, 1, ind.PSI.Dev_bg]) -
    sum(sum(X.PSI[1:n.grid, 1, ind.PSI.Dev_lo])))
  #______________________________________________________#
  
  #####################################
  # Intermediate mechanistic pathways #
  #####################################
  
  # Grid cell scale pathways #
  for(k in 1:n.grid) {
    for(t in 1:n.year) {
      X.PSI.raw[k, t, ind.WellD_3km] ~ dnegbin(pred.WellD_3km[k, t], r.WellD_3km)
      pred.WellD_3km[k, t] <- r.WellD_3km / (r.WellD_3km + lambda.WellD_3km[k, t])
      log(lambda.WellD_3km[k, t]) <- ALPHA0.WellD_3km +
        ALPHA.Dev_lo.WellD_3km * X.PSI[k, t, ind.PSI.Dev_lo] +
        ALPHA.Dev_bg.WellD_3km * X.PSI[k, t, ind.PSI.Dev_bg]
      #____ GOF _____#
      LLobs.WellD_3km[k, t] <- logdensity.negbin(X.PSI.raw[k, t, ind.WellD_3km],
        pred.WellD_3km[k, t], r.WellD_3km)
      X.sim.WellD_3km[k, t] ~ dnegbin(pred.WellD_3km[k, t], r.WellD_3km)
      LLsim.WellD_3km[k, t] <- logdensity.negbin(X.sim.WellD_3km[k, t],
        pred.WellD_3km[k, t], r.WellD_3km)
      #______________#
        
      X.PSI.raw[k, t, ind.WellA_1km] ~ dnegbin(pred.WellA_1km[k, t], r.WellA_1km)
      pred.WellA_1km[k, t] <- r.WellA_1km / (r.WellA_1km + lambda.WellA_1km[k, t])
      log(lambda.WellA_1km[k, t]) <- ALPHA0.WellA_1km +
        ALPHA.Dev_lo.WellA_1km * X.PSI[k, t, ind.PSI.Dev_lo] +
        ALPHA.Dev_bg.WellA_1km * X.PSI[k, t, ind.PSI.Dev_bg]
      #____ GOF _____#
      LLobs.WellA_1km[k, t] <- logdensity.negbin(X.PSI.raw[k, t, ind.WellA_1km],
        pred.WellA_1km[k, t], r.WellA_1km)
      X.sim.WellA_1km[k, t] ~ dnegbin(pred.WellA_1km[k, t], r.WellA_1km)
      LLsim.WellA_1km[k, t] <- logdensity.negbin(X.sim.WellA_1km[k, t],
        pred.WellA_1km[k, t], r.WellA_1km)
      #______________#
    }
    
    for(t in 1:2) { # Years 2-n.year are redundant.
      X.PSI.raw[k, t, ind.Road_1km] ~ dgamma(shape.Road_1km,
        shape.Road_1km / exp(pred.Road_1km[k, t]))
      pred.Road_1km[k, t] <- ALPHA0.Road_1km +
        ALPHA.Dev_lo.Road_1km * X.PSI[k, t, ind.PSI.Dev_lo] +
        ALPHA.Dev_bg.Road_1km * X.PSI[k, t, ind.PSI.Dev_bg]
      #____ GOF _____#
      LLobs.Road_1km[k, t] <- logdensity.gamma(X.PSI.raw[k, t, ind.Road_1km],
        shape.Road_1km, shape.Road_1km / exp(pred.Road_1km[k, t]))
      X.sim.Road_1km[k, t] ~ dgamma(shape.Road_1km,
        shape.Road_1km / exp(pred.Road_1km[k, t]))
      LLsim.Road_1km[k, t] <- logdensity.gamma(X.sim.Road_1km[k, t],
        shape.Road_1km, shape.Road_1km / exp(pred.Road_1km[k, t]))
      #______________#
    }
  }
    
  # Point scale pathways #
  for(j in 1:n.point) {
    for(t in 1:n.year) {
      X.psi.raw[j, t, ind.WellA_125m] ~ dbern(pred.WellA_125m[j, t])
      logit(pred.WellA_125m[j, t]) <- alpha0.WellA_125m +
        alpha.Dev_lo.WellA_125m * X.PSI[gridID[j], t, ind.PSI.Dev_lo] +
        alpha.Dev_bg.WellA_125m * X.PSI[gridID[j], t, ind.PSI.Dev_bg]
      #____ GOF _____#
      LLobs.WellA_125m[j, t] <- logdensity.bern(X.psi.raw[j, t, ind.WellA_125m],
        pred.WellA_125m[j, t])
      X.sim.WellA_125m[j, t] ~ dbern(pred.WellA_125m[j, t])
      LLsim.WellA_125m[j, t] <- logdensity.bern(X.sim.WellA_125m[j, t],
        pred.WellA_125m[j, t])
      #______________#
      
      X.psi.raw[j, t, ind.AHerb] ~ dbeta(a.AHerb[j, t], b.AHerb[j, t])
      a.AHerb[j, t] <- pred.AHerb[j, t] * phi.AHerb
      b.AHerb[j, t] <- (1 - pred.AHerb[j, t]) * phi.AHerb
      logit(pred.AHerb[j, t]) <- alpha0.AHerb +
        alpha.WellA_125m.AHerb * X.psi.raw[j, t, ind.WellA_125m] +
        alpha.Road_125m.AHerb * X.psi.raw[j, t, ind.Road_125m]
      #____ GOF _____#
      LLobs.AHerb[j, t] <- logdensity.beta(X.psi.raw[j, t, ind.AHerb],
        a.AHerb[j, t], b.AHerb[j, t])
      X.sim.AHerb[j, t] ~ dbeta(a.AHerb[j, t], b.AHerb[j, t])
      LLsim.AHerb[j, t] <- logdensity.beta(X.sim.AHerb[j, t],
        a.AHerb[j, t], b.AHerb[j, t])
      #______________#
    }
    
    for(t in 1:2) { # Years 2-n.year are redundant.
      X.psi.raw[j, t, ind.D_Road_125m] ~ dgamma(shape.D_Road_125m,
        shape.D_Road_125m / exp(pred.D_Road_125m[j, t])) # Need to compile in R. Should be NA where roads are absent.
      X.psi.raw[j, t, ind.PA_Road_125m] ~ dbern(pred.PA_Road_125m[j, t]) # Need to compile this in R.
      logit(pred.PA_Road_125m[j, t]) <- alpha0.PA_Road_125m +
        alpha.Dev_lo.PA_Road_125m * X.PSI[gridID[j], t, ind.PSI.Dev_lo] +
        alpha.Dev_bg.PA_Road_125m * X.PSI[gridID[j], t, ind.PSI.Dev_bg]
      pred.D_Road_125m[j, t] <- alpha0.D_Road_125m +
        alpha.Dev_lo.D_Road_125m * X.PSI[gridID[j], t, ind.PSI.Dev_lo] +
        alpha.Dev_bg.D_Road_125m * X.PSI[gridID[j], t, ind.PSI.Dev_bg]
      #____ GOF _____#
      LLobs.Road_125m[j, t] <- logdensity.bern(X.psi.raw[j, t, ind.PA_Road_125m],
        pred.PA_Road_125m[j, t]) +
        ifelse(X.psi.raw[j, t, ind.PA_Road_125m] == 1,
          logdensity.gamma(X.psi.raw[j, t, ind.D_Road_125m],
          shape.D_Road_125m, shape.D_Road_125m / exp(pred.D_Road_125m[j, t])),
          0)
      X.sim.PA_Road_125m[j, t] ~ dbern(pred.PA_Road_125m[j, t])
      X.sim.D_Road_125m[j, t] ~ dgamma(shape.D_Road_125m,
        shape.D_Road_125m / exp(pred.D_Road_125m[j, t]))
      LLsim.Road_125m[j, t] <- logdensity.bern(X.sim.PA_Road_125m[j, t],
        pred.PA_Road_125m[j, t]) +
        ifelse(X.sim.PA_Road_125m[j, t] == 1,
          logdensity.gamma(X.sim.D_Road_125m[j, t],
          shape.D_Road_125m, shape.D_Road_125m / exp(pred.D_Road_125m[j, t])),
          0)
      #______________#
    }
  }
  
  #________ GOF __________#
  dev.obs.WellD_3km <- -2 * sum(LLobs.WellD_3km[1:n.grid, 1:n.year])
  dev.sim.WellD_3km <- -2 * sum(LLsim.WellD_3km[1:n.grid, 1:n.year])
  test.WellD_3km <- step(dev.sim.WellD_3km - dev.obs.WellD_3km)

  dev.obs.WellA_1km <- -2 * sum(LLobs.WellA_1km[1:n.grid, 1:n.year])
  dev.sim.WellA_1km <- -2 * sum(LLsim.WellA_1km[1:n.grid, 1:n.year])
  test.WellA_1km <- step(dev.sim.WellA_1km - dev.obs.WellA_1km)
  
  dev.obs.Road_1km <- -2 * sum(LLobs.Road_1km[1:n.grid, 1:2])
  dev.sim.Road_1km <- -2 * sum(LLsim.Road_1km[1:n.grid, 1:2])
  test.Road_1km <- step(dev.sim.Road_1km - dev.obs.Road_1km)

  dev.obs.WellA_125m <- -2 * sum(LLobs.WellA_125m[1:n.point, 1:n.year])
  dev.sim.WellA_125m <- -2 * sum(LLsim.WellA_125m[1:n.point, 1:n.year])
  test.WellA_125m <- step(dev.sim.WellA_125m - dev.obs.WellA_125m)

  dev.obs.Road_125m <- -2 * sum(LLobs.Road_125m[1:n.point, 1:2])
  dev.sim.Road_125m <- -2 * sum(LLsim.Road_125m[1:n.point, 1:2])
  test.Road_125m <- step(dev.sim.Road_125m - dev.obs.Road_125m)

  dev.obs.AHerb <- -2 * sum(LLobs.AHerb[1:n.point, 1:n.year])
  dev.sim.AHerb <- -2 * sum(LLsim.AHerb[1:n.point, 1:n.year])
  test.AHerb <- step(dev.sim.AHerb - dev.obs.AHerb)
  #_______________________#
  
  ### Prior distributions ###
  ## Grid cell level models ##
  ALPHA0.WellD_3km ~ dnorm(0, 0.66667)
  ALPHA.Dev_lo.WellD_3km ~ dnorm(0, 0.66667)
  ALPHA.Dev_bg.WellD_3km ~ dnorm(0, 0.66667)
  r.WellD_3km ~ dunif(0, 50)

  ALPHA0.WellA_1km ~ dnorm(0, 0.66667)
  ALPHA.Dev_lo.WellA_1km ~ dnorm(0, 0.66667)
  ALPHA.Dev_bg.WellA_1km ~ dnorm(0, 0.66667)
  r.WellA_1km ~ dunif(0, 50)

  ALPHA0.Road_1km ~ dnorm(0, 0.5)
  ALPHA.Dev_lo.Road_1km ~ dnorm(0, 0.5)
  ALPHA.Dev_bg.Road_1km ~ dnorm(0, 0.5)
  shape.Road_1km ~ dunif(0, 100)

  ## Point level models ##
  alpha0.WellA_125m ~ dnorm(0, 0.66667)
  alpha.Dev_lo.WellA_125m ~ dnorm(0, 0.66667)
  alpha.Dev_bg.WellA_125m ~ dnorm(0, 0.66667)

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
