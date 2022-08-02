## eval_1000_k1_bw.R
##
## distribution: Gauss
## k=2 GARMA(1,1) process
##
## Length: 1000
##
## Bayesian (Stan) estimation using Whittle Likelihood

library(rstan)
library(rstantools)
library(bayesplot)
setwd("~/gauss_pq_k1_1000")
options(mc.cores=parallel::detectCores())

synProcess100 <- readRDS('syn_k1gauss1000.rds')


bayes_whittle<-function(y,runs=2000,...) {
  ss<-stats::spectrum(as.numeric(y),plot=FALSE,detrend=TRUE,demean=TRUE,method='pgram',taper=0,fast=FALSE)
  spec <- ss$spec/(2*pi)
  freq <- ss$freq
  
  stan_params <- list(N=length(spec),spec=spec,freq=freq)
  fit <- stan(
    file = "whittle_est_u.stan", # Stan program
    data = stan_params,        # named list of data
    chains = 4,                # number of Markov chains
    warmup = 1000,              # number of warmup iterations per chain
    iter = runs,               # total number of iterations per chain
    cores = 4,                 # number of cores (could use one per chain)
    seed=1795523,
    refresh = 0,                # no progress shown
    ...
  )
  return(fit)
}

# compile and test
y<-synProcess100$y30
system.time(fit <- bayes_whittle(y,runs=1500,control=list(adapt_delta=0.4)))
fit

res <-data.frame(bw_d1=numeric(1000),bw_f1=numeric(1000),bw_phi=numeric(1000),bw_theta=numeric(1000), bw_t=numeric(1000),
                 bw_d1_rhat=numeric(1000),bw_f1_rhat=numeric(1000),bw_phi_rhat=numeric(1000),bw_theta_rhat=numeric(1000),
                 bw_lp_rhat=numeric(1000),
                 bw_d1_neff=numeric(1000),bw_f1_neff=numeric(1000),bw_phi_neff=numeric(1000),bw_theta_neff=numeric(1000),
                 bw_lp_neff=numeric(1000),
                 bw_divergent=numeric(1000), bw_max_treedepth=numeric(1000),
                 bw_bfmi1=numeric(1000),bw_bfmi2=numeric(1000),bw_bfmi3=numeric(1000),bw_bfmi4=numeric(1000),
                 bw_d1_bulk_ess=numeric(1000),bw_f1_bulk_ess=numeric(1000),bw_phi_bulk_ess=numeric(1000),bw_theta_bulk_ess=numeric(1000),
                 bw_lp_bulk_ess=numeric(1000),
                 bw_d1_valid=numeric(1000),bw_f1_valid=numeric(1000),bw_phi_valid=numeric(1000),bw_theta_valid=numeric(1000),
                 bw_lp_valid=numeric(1000),
                 bw_attempt=numeric(1000),bw_adapt_delta=numeric(1000))


for (i in 1:1000) {
  y<-synProcess100[,paste0("y",i)]
  tt<-system.time(fit<-bayes_whittle(y,control=list(adapt_delta=0.3)))
  res$bw_attempt[i] <- 1
  res$bw_adapt_delta[i] <- 0.3
  res1 <- monitor(extract(fit, permuted = FALSE, inc_warmup = TRUE), print=FALSE) # matrix
  if (get_num_divergent(fit)>0) {
    tt<-system.time(fit<-bayes_whittle(y,control=list(adapt_delta=0.4)))
    res$bw_attempt[i] <- 2
    res$bw_adapt_delta[i] <- 0.4
    res1 <- monitor(extract(fit, permuted = FALSE, inc_warmup = TRUE), print=FALSE) # matrix
    if (get_num_divergent(fit)>0) {
      tt<-system.time(fit<-bayes_whittle(y,control=list(adapt_delta=0.5)))
      res$bw_attempt[i] <- 3
      res$bw_adapt_delta[i] <- 0.5
      res1 <- monitor(extract(fit, permuted = FALSE, inc_warmup = TRUE), print=FALSE) # matrix
      if (get_num_divergent(fit)>0) {
        tt<-system.time(fit<-bayes_whittle(y,control=list(adapt_delta=0.6)))
        res$bw_attempt[i] <- 4
        res$bw_adapt_delta[i] <- 0.6
        res1 <- monitor(extract(fit, permuted = FALSE, inc_warmup = TRUE), print=FALSE) # matrix
      }
    }
  }
  res$bw_t[i]<-tt[['elapsed']]
  res$bw_d1[i] <- res1["d1","mean"]
  res$bw_f1[i] <- res1["f1","mean"]
  res$bw_phi[i] <- res1["phi","mean"]
  res$bw_theta[i] <- res1["theta","mean"]
  res$bw_d1_rhat[i] <- res1["d1","Rhat"]
  res$bw_f1_rhat[i] <- res1["f1","Rhat"]
  res$bw_phi_rhat[i] <- res1["phi","Rhat"]
  res$bw_theta_rhat[i] <- res1["theta","Rhat"]
  res$bw_lp_rhat[i] <- res1["lp__","Rhat"]
  res$bw_d1_neff[i] <- res1["d1","n_eff"]
  res$bw_f1_neff[i] <- res1["f1","n_eff"]
  res$bw_phi_neff[i] <- res1["phi","n_eff"]
  res$bw_theta_neff[i] <- res1["theta","n_eff"]
  res$bw_lp_neff[i] <- res1["lp__","n_eff"]
  res$bw_divergent[i] <- get_num_divergent(fit)
  res$bw_max_treedepth[i] <- get_num_max_treedepth(fit)
  res$bw_d1_bulk_ess[i] <- res1["d1","Bulk_ESS"]
  res$bw_f1_bulk_ess[i] <- res1["f1","Bulk_ESS"]
  res$bw_phi_bulk_ess[i] <- res1["phi","Bulk_ESS"]
  res$bw_theta_bulk_ess[i] <- res1["theta","Bulk_ESS"]
  res$bw_lp_bulk_ess[i] <- res1["lp__","Bulk_ESS"]
  res$bw_d1_valid[i] <- res1["d1","valid"]
  res$bw_f1_valid[i] <- res1["f1","valid"]
  res$bw_phi_valid[i] <- res1["phi","valid"]
  res$bw_theta_valid[i] <- res1["theta","valid"]
  res$bw_lp_valid[i] <- res1["lp__","valid"]
  bfmi <- get_bfmi(fit)
  res$bw_bfmi1[i] <- bfmi[1]
  res$bw_bfmi2[i] <- bfmi[2]
  res$bw_bfmi3[i] <- bfmi[3]
  res$bw_bfmi4[i] <- bfmi[4]
  
  saveRDS(res,"res_gauss_1000_k1_ARMA_bw.RDS")
  cat(paste("Finished ",i,".\n"))
}

res <- readRDS("res_gauss_1000_k1_ARMA_bw.RDS")
for (i in 1:1000) if (res$bw_divergent[i]>0) {
  y<-synProcess100[,paste0("y",i)]
  tt<-system.time(fit<-bayes_whittle(y,control=list(adapt_delta=0.8)))
  res$bw_attempt[i] <- 5
  res$bw_adapt_delta[i] <- 0.8
  res1 <- monitor(extract(fit, permuted = FALSE, inc_warmup = TRUE), print=FALSE) # matrix
  res$bw_t[i]<-tt[['elapsed']]
  res$bw_d1[i] <- res1["d1","mean"]
  res$bw_f1[i] <- res1["f1","mean"]
  res$bw_phi[i] <- res1["phi","mean"]
  res$bw_theta[i] <- res1["theta","mean"]
  res$bw_d1_rhat[i] <- res1["d1","Rhat"]
  res$bw_f1_rhat[i] <- res1["f1","Rhat"]
  res$bw_phi_rhat[i] <- res1["phi","Rhat"]
  res$bw_theta_rhat[i] <- res1["theta","Rhat"]
  res$bw_lp_rhat[i] <- res1["lp__","Rhat"]
  res$bw_d1_neff[i] <- res1["d1","n_eff"]
  res$bw_f1_neff[i] <- res1["f1","n_eff"]
  res$bw_phi_neff[i] <- res1["phi","n_eff"]
  res$bw_theta_neff[i] <- res1["theta","n_eff"]
  res$bw_lp_neff[i] <- res1["lp__","n_eff"]
  res$bw_divergent[i] <- get_num_divergent(fit)
  res$bw_max_treedepth[i] <- get_num_max_treedepth(fit)
  res$bw_d1_bulk_ess[i] <- res1["d1","Bulk_ESS"]
  res$bw_f1_bulk_ess[i] <- res1["f1","Bulk_ESS"]
  res$bw_phi_bulk_ess[i] <- res1["phi","Bulk_ESS"]
  res$bw_theta_bulk_ess[i] <- res1["theta","Bulk_ESS"]
  res$bw_lp_bulk_ess[i] <- res1["lp__","Bulk_ESS"]
  res$bw_d1_valid[i] <- res1["d1","valid"]
  res$bw_f1_valid[i] <- res1["f1","valid"]
  res$bw_phi_valid[i] <- res1["phi","valid"]
  res$bw_theta_valid[i] <- res1["theta","valid"]
  res$bw_lp_valid[i] <- res1["lp__","valid"]
  bfmi <- get_bfmi(fit)
  res$bw_bfmi1[i] <- bfmi[1]
  res$bw_bfmi2[i] <- bfmi[2]
  res$bw_bfmi3[i] <- bfmi[3]
  res$bw_bfmi4[i] <- bfmi[4]
  
  saveRDS(res,"res_gauss_1000_k1_ARMA_bw.RDS")
  cat(paste("Finished ",i,".\n"))
}

