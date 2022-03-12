## eval_500_k2_wll.R
##
## distribution: ChiSq
## k=2 GAR(1) process
##
## Length: 500
##
## WLL estimation

library(pracma)
library(Rsolnp)
library(nloptr)

setwd("~/chisq_k2_500")


# generic optimize function
sim_optim<-function(par, objective, lower, upper, ...) {
  # check pars
  for (i in 1:length(par))
    if (par[i]<=lower[i]|par[i]>=upper[i]) par[i]<- (lower[i]+upper[i])/2   # if below lower bound then set to middle value
  
  res <- solnp(par, objective, LB=lower, UB=upper, control=list(tol=1e-12,trace=0), ...)
  if (res$convergence!=0)
    lbfgs(par, objective, lower=lower, upper=upper, control=list(maxeval=500), ...)
  else res$par <- res$pars
  return(res)
}

# Calculate Gegenbauer coeficients
ggbr.coef<-function(n,d,eta) {
  cf<-c(1,2*d*eta,2*d*(d+1)*eta^2-d)
  for (j in 3:(n-1)) cf[j+1]<-2*eta*((d-1)/j+1)*cf[j]-(2*(d-1)/j+1)*cf[j-1]
  return(cf)
}

wll.ggbr.obj<-function(par,ss) {
  fd1   <- par[1]
  f1    <- par[2]
  fd2   <- par[3]
  f2    <- par[4]
  phi   <- par[5]
  
  cos_2_pi_f <- cos(2.0*pi*ss$freq)
  u1 <- cos(2*pi*f1)
  u2 <- cos(2*pi*f2)
  
  excl_idx1 <- which.min(abs(u1-cos_2_pi_f))
  excl_idx2 <- which.min(abs(u2-cos_2_pi_f))
  keep_idx <- setdiff(1:length(ss$spec),c(excl_idx1,excl_idx2))
  cos_2_pi_f <- cos_2_pi_f[keep_idx]
  spec <- ss$spec[keep_idx]
  
  mod_phi <- (1+phi^2-2*phi*cos_2_pi_f)
  spec_den_inv <- 2.0*pi * (4.0*((cos_2_pi_f-u1)^2))^fd1 * (4.0*((cos_2_pi_f-u2)^2))^fd2 * mod_phi   # Inverse of spectral density
  
  spec_den_inv[is.infinite(spec_den_inv)] <- NA
  spec_den_inv[spec_den_inv<=0] <- NA
  I_f <- spec*spec_den_inv
  res <- sum((log(I_f))^2,na.rm=TRUE)
  return(res)
}


# find 2 largest spikes in spectrum
big_spikes<-function(ss) {
  # find freq 0.1
  start_idx <- 1
  end_idx <- which.min(abs(ss$freq-0.4))
  spec <- ss$spec[start_idx:end_idx]
  freq <- ss$freq[start_idx:end_idx]
  
  spike1 <- which.max(spec)
  for (i in (spike1-5):(spike1+5)) if (i>0&i<length(spec)) spec[i] <- 0
  spike2 <- which.max(spec)
  return(list(f1=freq[spike1],f2=freq[spike2]))
}

wll.est<-function(y) {
  ss<-spectrum(y,plot=FALSE,detrend=FALSE,demean=FALSE,method='pgram',taper=0,fast=FALSE)
  spikes <- big_spikes(ss)
  
  lb<-c(0.01,0.01,0.01,0.01,-0.999)
  ub<-c(0.49,0.4,0.49,0.4, 0.999)
  pars<-c(0.25,spikes[['f1']],0.25,spikes[['f2']],0.0) #initial.est(y)
  fit <- sim_optim(pars, wll.ggbr.obj, lower=lb, upper=ub, ss=ss)
  return(fit$par)  
}

synProcess500 <- readRDS('syn_k2chisq500.rds')

res <- data.frame(wll_d1=numeric(1000),wll_f1=numeric(1000),
                  wll_d2=numeric(1000),wll_f2=numeric(1000),
                  wll_phi=numeric(1000), wll_t=numeric(1000))

for (i in 1:1000) {
  y<-synProcess500[,paste0("y",i)]
  tt<-system.time(wll<-wll.est(y))
  # sort peaks into correct order
  if (wll[2]>wll[4]) {
    res$wll_d1[i]<-wll[3]
    res$wll_f1[i]<-wll[4]
    res$wll_d2[i]<-wll[1]
    res$wll_f2[i]<-wll[2]
  } else {
    res$wll_d1[i]<-wll[1]
    res$wll_f1[i]<-wll[2]
    res$wll_d2[i]<-wll[3]
    res$wll_f2[i]<-wll[4]
  }

  res$wll_phi[i]<-wll[5]
  res$wll_t[i]<-tt[['elapsed']]
  saveRDS(res,"res_chisq_500_k2_wll.RDS")
  if ((i%%25)==0) cat(paste("running ",i,".\n"))
}

saveRDS(res,"res_chisq_500_k2_wll.RDS")


