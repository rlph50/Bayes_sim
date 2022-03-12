
library(ltsa)
library(pracma)
setwd("~/gauss_k2_1000")
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

#VecGgbrAR1SpecDen<-Vectorize(GgbrAR1SpecDen,vectorize.args ="x")

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



quadgk_custom <- function(f, a, b, tol = .Machine$double.eps^0.5, ...) {
  stopifnot(is.numeric(a), length(a) == 1,
            is.numeric(b), length(b) == 1)
  eps <- .Machine$double.eps
  
  fun <- match.fun(f)
  f <- function(x) fun(x, ...)
  
  if (a == b)     return(0)
  else if (a > b) return(-1 * quadgk(f, b, a, tol = tol))
  
  # Nodes and weights for Gauss-Kronrod (7, 15)
  n15 <- c(-0.9914553711208126, -0.9491079123427585, -0.8648644233597691,
           -0.7415311855993944, -0.5860872354676911, -0.4058451513773972,
           -0.2077849550078985,  0.0,                 0.2077849550078985,
           0.4058451513773972,  0.5860872354676911,  0.7415311855993944,
           0.8648644233597691,  0.9491079123427585,  0.9914553711208126)
  n7  <- c(-0.9491079123427585, -0.7415311855993944, -0.4058451513773972,
           0.0,
           0.4058451513773972, 0.7415311855993944,  0.9491079123427585)
  
  w15 <- c(0.02293532201052922, 0.06309209262997855,  0.1047900103222502,
           0.1406532597155259,  0.1690047266392679,   0.1903505780647854,
           0.2044329400752989,  0.2094821410847278,   0.2044329400752989,
           0.1903505780647854,  0.1690047266392679,   0.1406532597155259,
           0.1047900103222502,  0.06309209262997855,  0.02293532201052922)
  w7  <- c(0.1294849661688697,  0.2797053914892767,   0.3818300505051189,
           0.4179591836734694,
           0.3818300505051189,  0.2797053914892767,   0.1294849661688697)
  
  .gkadpt <- function(f, a, b, tol = tol) {
    # use nodes and weights from the environment
    x15 <- 0.5 * ((b - a) * n15 + b + a)
    x7  <- 0.5 * ((b - a) * n7  + b + a)
    
    f7<-f(x7)
    if (length(which(is.infinite(f7)))) {
      n<-which(is.infinite(f7))
      for (nn in n) if (nn==7) f7[7]<-f7[6] else f7[nn]<-f7[nn+1]
    }
    f15<-f(x15)
    if (length(which(is.infinite(f15)))) {
      n<-which(is.infinite(f15))
      for (nn in n) if (nn==15) f15[15]<-f15[14] else f15[nn]<-f15[nn+1]
    }
    
    Q7  <- sum(w7  * f7)  * (b-a)/2
    Q15 <- sum(w15 * f15) * (b-a)/2
    
    if (!is.finite(Q7) || !is.finite(Q15)) {
      #warning("Infinite or NA function value encountered.")
      return(Q15)
    } else if (abs(Q15 - Q7) < tol) {
      return(Q15)
    } else if (abs(b-a) < 16*eps) {
      #warning("Minimum step size reached; singularity possible.")
      return(Q15)
    } # else
    
    Q2 <- .gkadpt(f, (a+b)/2, b, tol = tol)
    Q1 <- .gkadpt(f, a, (a+b)/2, tol = tol)
    
    return(Q1 + Q2)
  }
  
  # start the recursive procedure
  .gkadpt(f, a, b, tol = tol)
}


# Get acf of theoretical spectral density using numerical integration
ggbr.sim<-function(n,ll2) {
  # get acf by numerical integration
  g<-rep(0,n)
  for (k in 1:n) {
    ll2$k <- k-1
    # g[k]<-quadgk(GgbrAR1SpecDen, 0, 2*pi, ll=ll2)
    # g[k]<-quadgk_custom(GgbrAR1SpecDen, 0, 2*pi, ll=ll2)
    # g[k]<-quadgr(GgbrAR1SpecDen, 0, 2*pi,ll=ll2)$value
    g[k]<-quad_trap(GgbrAR1SpecDen, 0, 2*pi,ll=ll2)
  }
  #for (k in 1:n) g[k]<-quadgr(VecGgbrAR1SpecDen, 0, 2*pi)$value
  g<-unlist(g)
  g<-g/g[1] #max(g)
  return(g)
}

ll1=list(d1=0.30,f1=1/28,d2=0.20,f2=1/7,phi=(-0.8),k=0,sigma2=1)

n<-1000
system.time(g <- ggbr.sim(n,ll2=ll1))
#plot(1:500,fft(g)[1:500],type='l')

# # data frame to store the simulated processes
syn_gauss1000<-data.frame(obs=1:n)
#periodograms for each of the above
sd1000<-data.frame(obs=1:as.integer(n/2))
i<-1
while (i<=1000) {
  # from "waveslim" package, use method of Hosking (1984) to generate a realization
  y2<-DHSimulate(2*n,g[1:n],rand.gen = rnorm)
  syn_gauss1000[,paste0("y",i)]<-y2[(n+1):(2*n)]
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


saveRDS(syn_gauss1000,'syn_k2gauss1000.rds')


