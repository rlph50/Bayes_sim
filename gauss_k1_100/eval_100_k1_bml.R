## eval_100_k1_bml.R
##
## distribution: Gauss
## k=1 GAR(1) process
##
## Length: 100
##
## Bayesian Maximum Likelihood estimation

library(rstan)
library(rstantools)
library(bayesplot)
setwd("~/gauss_k1_100")
options(mc.cores=parallel::detectCores())

synProcess <- readRDS('syn_k1gauss100.rds')


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
system.time(fit <- bayes_ml(y,control=list(adapt_delta=0.60)))
fit

res <-data.frame(bml_d=numeric(1000),bml_f=numeric(1000),bml_phi=numeric(1000), bml_t=numeric(1000),
                 bml_d_rhat=numeric(1000),bml_f_rhat=numeric(1000),bml_phi_rhat=numeric(1000),bml_lp_rhat=numeric(1000),
                 bml_d_neff=numeric(1000),bml_f_neff=numeric(1000),bml_phi_neff=numeric(1000),bml_lp_neff=numeric(1000),
                 bml_divergent=numeric(1000), bml_max_treedepth=numeric(1000),
                 bml_bfmi1=numeric(1000),bml_bfmi2=numeric(1000),bml_bfmi3=numeric(1000),bml_bfmi4=numeric(1000),
                 bml_d_bulk_ess=numeric(1000),bml_f_bulk_ess=numeric(1000),bml_phi_bulk_ess=numeric(1000),bml_lp_bulk_ess=numeric(1000),
                 bml_d_neff=numeric(1000),bml_f_neff=numeric(1000),bml_phi_neff=numeric(1000),bml_lp_neff=numeric(1000),
                 bml_d_valid=numeric(1000),bml_f_valid=numeric(1000),bml_phi_valid=numeric(1000),bml_lp_valid=numeric(1000))

for (i in 1:1000) {
  cat(sprintf("processing y%d.\n",i))
  y<-synProcess[,paste0("y",i)]
  tt<-system.time(fit<-bayes_ml(y,control=list(adapt_delta=0.60)))
  res1 <- monitor(extract(fit, permuted = FALSE, inc_warmup = TRUE), print=FALSE) # matrix
  res$bml_t[i]<-tt[['elapsed']]
  res$bml_d[i] <- res1["d","mean"]
  res$bml_f[i] <- res1["f","mean"]
  res$bml_phi[i] <- res1["phi","mean"]
  res$bml_d_rhat[i] <- res1["d","Rhat"]
  res$bml_f_rhat[i] <- res1["f","Rhat"]
  res$bml_phi_rhat[i] <- res1["phi","Rhat"]
  res$bml_lp_rhat[i] <- res1["lp__","Rhat"]
  res$bml_d_neff[i] <- res1["d","n_eff"]
  res$bml_f_neff[i] <- res1["f","n_eff"]
  res$bml_phi_neff[i] <- res1["phi","n_eff"]
  res$bml_lp_neff[i] <- res1["lp__","n_eff"]
  res$bml_divergent[i] <- get_num_divergent(fit)
  res$bml_max_treedepth[i] <- get_num_max_treedepth(fit)
  res$bml_d_bulk_ess[i] <- res1["d","Bulk_ESS"]
  res$bml_f_bulk_ess[i] <- res1["f","Bulk_ESS"]
  res$bml_phi_bulk_ess[i] <- res1["phi","Bulk_ESS"]
  res$bml_lp_bulk_ess[i] <- res1["lp__","Bulk_ESS"]
  res$bml_d_valid[i] <- res1["d","valid"]
  res$bml_f_valid[i] <- res1["f","valid"]
  res$bml_phi_valid[i] <- res1["phi","valid"]
  res$bml_lp_valid[i] <- res1["lp__","valid"]
  bfmi <- get_bfmi(fit)
  res$bml_bfmi1[i] <- bfmi[1]
  res$bml_bfmi2[i] <- bfmi[2]
  res$bml_bfmi3[i] <- bfmi[3]
  res$bml_bfmi4[i] <- bfmi[4]
  
  saveRDS(res,"res_gauss_100_k1_bml.RDS")
}

