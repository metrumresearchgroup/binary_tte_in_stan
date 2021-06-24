
data {
  int<lower=1> N;  // total number of observations
  int<lower=0, upper=1> Y[N];  // binary response variable
  vector[N] x;  // predictor
  int prior_only;  // should the likelihood be ignored?
}


parameters {
  real alpha;  // intercept
  real beta ;  // slope
}

model {
  // log-likelihood contributions
  if (prior_only==0) {
    Y ~ bernoulli_logit(alpha + beta*x);
  }
  // prior didstributions
  beta ~ normal(0,3);
  alpha ~ normal(0,3);
}

generated quantities {
  vector[N] Ysim;  // response variable

  for (n in 1:N) {
    Ysim[n] = bernoulli_logit_rng(alpha + beta*x[n]);
  }
  
}

