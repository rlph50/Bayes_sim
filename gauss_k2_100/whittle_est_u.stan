// Whittle estimation

data {
  int<lower=1> N;            // num observations
  real spec[N];              // periodogram - but only from 0-pi not to 2pi.
  real freq[N];
}

transformed data {
  real two_pi;
  real cos_2_pi_f[N];

  two_pi = 2.0*pi();
  for (i in 1:N) cos_2_pi_f[i] = cos(two_pi*freq[i]);
}

parameters {
  real<lower=-0.999,upper=0.999> phi;  // autoregression coeff
  real<lower=0.01,upper=0.49> d1;
  real<lower=0.01,upper=0.49> d2;
  real<lower=0,upper=0.999> u1;  // limits correspond to f=0.1
  real<lower=u1+0.01,upper=0.9999> u2;  // limits correspond to f=0.1 u1<>u2
}

model {
  real log_likelihood;
  real f1;
  real f2;

  // Priors
  target += normal_lpdf(phi | 0.0, 0.50);
  target += uniform_lpdf(d1 | 0.01,0.49);
  target += uniform_lpdf(d2 | 0.01,0.49);
  target += normal_lpdf(u1 | 0.0,0.50);
  target += normal_lpdf(u2 | 0.0,0.50);

  // (approx.) log-likelihood.
  // Whittle Likelihood
  log_likelihood = 0.0;
  
  // next logic allows us to ensure we do not include the frequency cloests to  the current "u"
  f1 = acos(u1)/two_pi;
  if (f1<0.0) f1 += 0.5;
  if (f1>0.5) f1 -= 0.5;
  f2 = acos(u2)/two_pi;
  if (f2<0.0) f2 += 0.5;
  if (f2>0.5) f2 -= 0.5;

  for (i in 1:N) if ((fabs(freq[i]-f1)>0.001)&&(fabs(freq[i]-f2)>0.001)) {
    real spec_den_inv;
    real cos_diff1;
    real cos_diff2;
    
    cos_diff1 = cos_2_pi_f[i]-u1;
    cos_diff2 = cos_2_pi_f[i]-u2;
    spec_den_inv = two_pi * (1.0 + (phi^2) - 2.0*phi*cos_2_pi_f[i] );
    spec_den_inv *= (4.0*cos_diff1*cos_diff1)^d1;
    spec_den_inv *= (4.0*cos_diff2*cos_diff2)^d2;
    log_likelihood += ( (spec[i]*spec_den_inv)  - log(spec_den_inv)) ;
  }

  target += -(log_likelihood + log(two_pi));
}

generated quantities {
  real f1;
  real f2;

  f1 = acos(u1)/two_pi;
  if (f1<0.0) f1 += 0.5;
  if (f1>0.5) f1 -= 0.5;
  f2 = acos(u2)/two_pi;
  if (f2<0.0) f2 += 0.5;
  if (f2>0.5) f2 -= 0.5;
}
