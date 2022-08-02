data {
  int<lower=1> N;            // num observations
  real y[N];                 // observed outputs
}

transformed data {
  real two_pi;

  two_pi = 2.0*pi();
}

parameters {
  real<lower=-0.999,upper=0.999> phi;    // autoregression coeff
  real<lower=-0.999,upper=0.999> theta;  // moving average coeff
  real<lower=0.01,upper=0.49> d1;
  real<lower=0.01,upper=0.40> f1;
}

model {
  vector[N] nu;
  vector[N] err;
  real g[N+1];
  real u;
  
  u=cos(2.0*pi()*f1);
  
  g[1] = 1.0;                // this and next few lines to derive the Ggbr Coefficients.
  g[2] = 2.0*u*d1;
  for (j in 3:(N+1)) 
  {
    real rj1;
    real temp;
    
    rj1 = j-1;
    temp = (d1-1.0)/rj1;
    g[j]=(2.0*u*(rj1+d1-1.0)*g[j-1]-(rj1+2.0*d1-2.0)*g[j-2])/rj1;
  }

  nu[1] = 0.0;
  err[1] = y[1] - nu[1];
  for (t in 2:N) {
    real sum_g;
    sum_g = 0.0;
    for (j in 1:(t-1)) sum_g += g[j+1]*err[t-j];

    nu[t] = phi * y[t-1] - theta * err[t-1]  + sum_g;
    err[t] = y[t] - nu[t];
  }
  phi ~ normal(0.0, 0.5);
  theta ~ normal(0.0, 0.5);
  d1  ~ uniform(0.01,0.49);
  f1  ~ normal(0.25, 0.10);
  err ~ normal(0.0, 1.0);
}

