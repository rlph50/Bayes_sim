library(ltsa)
library(pracma)
setwd("~/chisq_k2_1000")
options(mc.cores=parallel::detectCores())

# Calculate Gegenbauer coeficients
ggbr.coef<-function(n,d,eta) {
  cf<-c(1,2*d*eta,2*d*(d+1)*eta^2-d)
  for (j in 3:(n-1)) cf[j+1]<-2*eta*((d-1)/j+1)*cf[j]-(2*(d-1)/j+1)*cf[j-1]
  return(cf)
}


# Generate Realizations
# Theoretical Spectral Density
GgbrAR1SpecDen<-function(x,ll) {return (cos(ll$k*x)*(4*(cos(x)-cos(2*pi*ll$f1))^2)^(-ll$d1)*(4*(cos(x)-cos(2*pi*ll$f2))^2)^(-ll$d2)/(1-2*ll$phi*cos(x)+ll$phi^2));}

# trapezoidal quadrature
quad_trap<-function(fcn,lb,ub,...) {
  delta_x <- (ub-lb)/31415
  grid <- seq(lb,ub,by=delta_x)
  y <- fcn(grid,...)

  inf_idx <- which(is.infinite(y))
  y[inf_idx] <- (y[inf_idx-1]+y[inf_idx+1])/2
  y_avg <- (y[2:length(y)]+y[1:(length(y)-1)])/2
  area <- sum(y_avg*delta_x)
  return(area)
}


# Get acf of theoretical spectral density using numerical integration
ggbr.sim<-function(n,ll2) {
  # get acf by numerical integration
  g<-rep(0,n)
  for (k in 1:n) {
    ll2$k <- k-1
    g[k]<-quad_trap(GgbrAR1SpecDen, 0, 2*pi,ll=ll2)
  }
  g<-unlist(g)
  g<-g/g[1] #max(g)
  return(g)
}

ll1=list(d1=0.30,f1=1/28,d2=0.20,f2=1/7,phi=(-0.8),k=0,sigma2=1)

n<-1000
system.time(g <- ggbr.sim(n,ll2=ll1))

rchisq1 <- function(n) return(rchisq(n,1)-1)

# # data frame to store the simulated processes
syn_chisq1000<-data.frame(obs=1:n)
#periodograms for each of the above
sd1000<-data.frame(obs=1:as.integer(n/2))
i<-1
while (i<=1000) {
  # from "waveslim" package, use method of Hosking (1984) to generate a realization
  y2<-DHSimulate(2*n,g[1:n],rand.gen = rchisq1)
  syn_chisq1000[,paste0("y",i)]<-y2[(n+1):(2*n)]
  ss<-spectrum(y2[(n+1):(2*n)],plot=F,detrend=F)
  # check the peak in the spectral density is high enough; if not throw this one away.
  idx1 <- which.min(abs((1/7)-ss$freq))
  idx2 <- which.min(abs((1/28)-ss$freq))
  if (ss$spec[idx1]>25&ss$spec[idx2]>25) {
    if (i==1) sd1000$freq<-ss$freq
    sd1000[,paste0("spec",i)]<-ss$spec
    cat(sprintf("run %d.\n",i))
    i<-i+1
  }
}
# test
with(sd1000,plot(freq,spec1,type='l'))


saveRDS(syn_chisq1000,'syn_k2chisq1000.rds')

