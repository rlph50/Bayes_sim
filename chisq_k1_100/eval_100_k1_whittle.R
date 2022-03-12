## eval_100_k1_whittle.R
##
## distribution: ChiSq
## k=1 GAR(1) process
##
## Length: 100
##
## Whittle estimation


library(pracma)
library(Rsolnp)
library(nloptr)

setwd("~/chisq_k1_100")

# generic optimize function
sim_optim<-function(par, objective, lower, upper, ...) {
  # check pars
  for (i in 1:length(par)) if (par[i]<=lower[i]|par[i]>=upper[i]) {
    par[i]<- (lower[i]+upper[i])/2   # if below lower bound then set to middle value
  }
  
  res <- solnp(par, objective, LB=lower, UB=upper, control=list(tol=1e-12,trace=0,outer.iter=3), ...)
  if (res$convergence!=0)
    res<-cobyla(par, objective, lower=lower, upper=upper, control=list(maxeval=500), ...)
  else res$par <- res$pars

  return(res)
}

# Calculate Gegenbauer coeficients
ggbr.coef<-function(n,d,eta) {
  cf<-c(1,2*d*eta,2*d*(d+1)*eta^2-d)
  for (j in 3:(n-1)) cf[j+1]<-2*eta*((d-1)/j+1)*cf[j]-(2*(d-1)/j+1)*cf[j-1]
  return(cf)
}


# find 2 largest spikes in spectrum
big_spikes<-function(ss) {
  # find freq 0.1
  start_idx <- 1
  end_idx <- which.min(abs(ss$freq-0.4))
  spec <- ss$spec[start_idx:end_idx]
  freq <- ss$freq[start_idx:end_idx]
  
  spike1 <- which.max(spec)
  #for (i in (spike1-5):(spike1+5)) if (i>0&i<length(spec)) spec[i] <- 0
  #spike2 <- which.max(spec)
  return(list(f1=freq[spike1]))
}

## Whittle Estimate
# Objective Function for Whittle method - GAR(1)
whittle.obj<-function(theta,ss) {
  fd1 <- theta[1]
  f1  <- theta[2]
  phi <- theta[3]
  u1  <- cos(2*pi*f1)
  cos_2_pi_f <- cos(2.0*pi*ss$freq)
  excl_idx1 <- which.min(abs(u1-cos_2_pi_f))
  keep_idx <- setdiff(1:length(ss$spec),excl_idx1)
  cos_2_pi_f <- cos_2_pi_f[keep_idx]
  spec <- ss$spec[keep_idx]
  
  mod_phi <- (1+phi^2-2*phi*cos_2_pi_f)
  spec_den_inv <- 2.0*pi * (4.0*((cos_2_pi_f-u1)^2))^fd1 * mod_phi   # Inverse of spectral density
  
  spec_den_inv[is.infinite(spec_den_inv)] <- NA
  spec_den_inv[spec_den_inv<=0] <-NA
  I_f <- spec*spec_den_inv
  res <- sum(I_f-log(spec_den_inv),na.rm=TRUE)
  return(res)
}

# Find Whittle Estimates for parameters GAR(1)
whittle.est<-function(y) {
  ss<-spectrum(y,plot=FALSE,detrend=FALSE,demean=FALSE,method='pgram',taper=0,fast=FALSE)
  
  spikes <- big_spikes(ss)
  fit<-sim_optim(c(0.25,spikes[['f1']],0.0),
                 whittle.obj,
                 lower=c(0.01,0.01,-0.999),
                 upper=c(0.49,0.40, 0.999),
                 ss=ss)
  return(fit$par)
}

synProcess100 <- readRDS('syn_k1chisq100.rds')

res <- data.frame(whittle_d1=numeric(1000),whittle_f1=numeric(1000),
                  whittle_phi=numeric(1000), whittle_t=numeric(1000))

#options(warn=-1)
for (i in 1:1000) {
  y<-synProcess100[,paste0("y",i)]
  tt<-system.time(whittle<-whittle.est(y))
  res$whittle_d1[i]<-whittle[1]
  res$whittle_f1[i]<-whittle[2]
  res$whittle_phi[i]<-whittle[3]
  res$whittle_t[i]<-tt[['elapsed']]
  if ((i%%25)==0) cat(paste("running ",i,".\n"))
}

saveRDS(res,"res_chisq_100_k1_whittle.RDS")

