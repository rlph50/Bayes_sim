data {
  int<lower=1> N;            // num observations
  real y[N];                 // observed outputs
}

transformed data {
  real two_pi;

  two_pi = 2.0*pi();
}

parameters {
  real<lower=-0.999,upper=0.999> phi;  // autoregression coeff
  real<lower=0.01,upper=0.49> d1;
  real<lower=0.01,upper=0.49> d2;
  real<lower=0.01,upper=0.2> f1;
  real<lower=f1+0.05,upper=0.4> f2;
}

model {
  vector[N] nu;
  vector[N] err;
  real g1[N+1];
  real g2[N+1];
  real g1g2[N+1];
  real u1;
  real u2;
  
  u1=cos(2.0*pi()*f1);
  u2=cos(2.0*pi()*f2);
  
  g1[1] = 1.0;                // this and next few lines to derive the Ggbr Coefficients.
  g1[2] = 2.0*u1*d1;
  g2[1] = 1.0;
  g2[2] = 2.0*u2*d2;
  for (j in 3:(N+1)) 
  {
    real rja;
    real temp;
    
    rja = j-1;
    temp = (d1-1.0)/rja;
    g1[j]=(2.0*u1*(rja+d1-1.0)*g1[j-1]-(rja+2.0*d1-2.0)*g1[j-2])/rja;
    temp = (d2-1.0)/rja;
    g2[j]=(2.0*u2*(rja+d2-1.0)*g2[j-1]-(rja+2.0*d2-2.0)*g2[j-2])/rja;
  }
  // Polynomial multiplication - first N terms only
  // g1g2 = g1*g2
  for (j in 1:(N+1)) {
    real sum1=0.0;
    for (i in 1:j) sum1 += g1[i]*g2[j-i+1];
    g1g2[j] = sum1;
  }

  nu[1] = 0.0;
  err[1] = y[1] - nu[1];
  for (t in 2:N) {
    real sum_g;
    sum_g = 0.0;
    for (j in 1:(t-1)) sum_g += g1g2[j+1]*err[t-j];

    nu[t] = phi * y[t-1] + sum_g;
    err[t] = y[t] - nu[t];
  }
  phi  ~ normal(0.0, 0.5);
  d1   ~ uniform(0.01,0.49);
  d2   ~ uniform(0.01,0.49);
  f1   ~ normal(0.25, 0.10);
  f2   ~ normal(0.25, 0.10);
  err  ~ normal(0.0, 1.0);
}


