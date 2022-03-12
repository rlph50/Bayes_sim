library(rstan)
library(rstantools)
library(bayesplot)
setwd("~/gauss_k2_1000")
options(mc.cores=parallel::detectCores())

synProcess <- readRDS('syn_k2gauss1000.rds')


bayes_ml<-function(y,runs=2000,...) {
  stan_params <- list(N=length(y),y=y)
  fit <- stan(
    file = "exact_f.stan", # Stan program
    data = stan_params,        # named list of data
    chains = 4,                # number of Markov chains
    warmup = 1000,             # number of warmup iterations per chain
    iter = runs,               # total number of iterations per chain
    cores = 4,                 # number of cores (could use one per chain)
    seed=1795523,
    refresh = 0,                # no progress shown
    ...
  )
  return(fit)
}

# compile and test
y<-synProcess$y1
system.time(fit <- bayes_ml(y,control=list(adapt_delta=0.6)))
fit

res <-data.frame(bml_d1=numeric(1000),bml_f1=numeric(1000),bml_d2=numeric(1000),bml_f2=numeric(1000),bml_phi=numeric(1000), bml_t=numeric(1000),
                 bml_d1_rhat=numeric(1000),bml_f1_rhat=numeric(1000),bml_d2_rhat=numeric(1000),bml_f2_rhat=numeric(1000),bml_phi_rhat=numeric(1000),bml_lp_rhat=numeric(1000),
                 bml_d1_neff=numeric(1000),bml_f1_neff=numeric(1000),bml_d2_neff=numeric(1000),bml_f2_neff=numeric(1000),bml_phi_neff=numeric(1000),bml_lp_neff=numeric(1000),
                 bml_divergent=numeric(1000), bml_max_treedepth=numeric(1000),
                 bml_bfmi1=numeric(1000),bml_bfmi2=numeric(1000),bml_bfmi3=numeric(1000),bml_bfmi4=numeric(1000),
                 bml_d1_bulk_ess=numeric(1000),bml_f1_bulk_ess=numeric(1000),bml_d2_bulk_ess=numeric(1000),bml_f2_bulk_ess=numeric(1000),bml_phi_bulk_ess=numeric(1000),
                 bml_lp_bulk_ess=numeric(1000),bml_d1_valid=numeric(1000),bml_f1_valid=numeric(1000),bml_d2_valid=numeric(1000),bml_f2_valid=numeric(1000),
                 bml_phi_valid=numeric(1000),bml_lp_valid=numeric(1000),bml_attempt=numeric(1000),bml_adapt_delta=numeric(1000))

#res<-readRDS("res_gauss_1000_k2_bml.RDS")
for (i in 1:1000) {
  cat(sprintf("processing y%d.\n",i))
  y<-synProcess[,paste0("y",i)]
  tt<-system.time(fit<-bayes_ml(y,control=list(adapt_delta=0.6)))
  res1 <- monitor(extract(fit, permuted = FALSE, inc_warmup = TRUE), print=FALSE) # matrix
  res$bml_adapt_delta[i] <- 0.6
  res$bml_t[i]<-tt[['elapsed']]
  res$bml_d1[i] <- res1["d1","mean"]
  res$bml_f1[i] <- res1["f1","mean"]
  res$bml_d2[i] <- res1["d2","mean"]
  res$bml_f2[i] <- res1["f2","mean"]
  res$bml_phi[i] <- res1["phi","mean"]
  res$bml_d1_rhat[i] <- res1["d1","Rhat"]
  res$bml_f1_rhat[i] <- res1["f1","Rhat"]
  res$bml_d2_rhat[i] <- res1["d2","Rhat"]
  res$bml_f2_rhat[i] <- res1["f2","Rhat"]
  res$bml_phi_rhat[i] <- res1["phi","Rhat"]
  res$bml_lp_rhat[i] <- res1["lp__","Rhat"]
  res$bml_d1_neff[i] <- res1["d1","n_eff"]
  res$bml_f1_neff[i] <- res1["f1","n_eff"]
  res$bml_d2_neff[i] <- res1["d2","n_eff"]
  res$bml_f2_neff[i] <- res1["f2","n_eff"]
  res$bml_phi_neff[i] <- res1["phi","n_eff"]
  res$bml_lp_neff[i] <- res1["lp__","n_eff"]
  res$bml_divergent[i] <- get_num_divergent(fit)
  res$bml_max_treedepth[i] <- get_num_max_treedepth(fit)
  res$bml_d1_bulk_ess[i] <- res1["d1","Bulk_ESS"]
  res$bml_f1_bulk_ess[i] <- res1["f1","Bulk_ESS"]
  res$bml_d2_bulk_ess[i] <- res1["d2","Bulk_ESS"]
  res$bml_f2_bulk_ess[i] <- res1["f2","Bulk_ESS"]
  res$bml_phi_bulk_ess[i] <- res1["phi","Bulk_ESS"]
  res$bml_lp_bulk_ess[i] <- res1["lp__","Bulk_ESS"]
  res$bml_d1_valid[i] <- res1["d1","valid"]
  res$bml_f1_valid[i] <- res1["f1","valid"]
  res$bml_d2_valid[i] <- res1["d2","valid"]
  res$bml_f2_valid[i] <- res1["f2","valid"]
  res$bml_phi_valid[i] <- res1["phi","valid"]
  res$bml_lp_valid[i] <- res1["lp__","valid"]
  bfmi <- get_bfmi(fit)
  res$bml_bfmi1[i] <- bfmi[1]
  res$bml_bfmi2[i] <- bfmi[2]
  res$bml_bfmi3[i] <- bfmi[3]
  res$bml_bfmi4[i] <- bfmi[4]
  
  saveRDS(res,"res_gauss_1000_k2_bml.RDS")
}

