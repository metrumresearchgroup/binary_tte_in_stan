
data {
  int<lower=1> N;  // total number of observations
  int<lower=0, upper=1> Y[N];  // binary response variable
  vector[N] x;  // predictor
  int prior_only;  // should the likelihood be ignored?
}


parameters {
  real alpha;  // intercept
  real<lower=0> beta ;  // slope
  real<lower=0> gamma; // power parameter
}

model {
  vector[N] eta;
  
  // log-likelihood contributions
  if (prior_only==0) {
    for (n in 1:N) {
      eta[n] = alpha + beta*pow(x[n],gamma);
    }
    Y ~ bernoulli_logit(eta);
  }
  // prior didstributions
  beta ~ normal(0,1);
  alpha ~ normal(0,1);
  gamma ~ normal(1,1);
}

generated quantities {
  vector[N] Ysim;  // response variable
  vector[N] log_lik; // log-likelihood, for calculating LOO
  
  for (n in 1:N) {
    Ysim[n] = bernoulli_logit_rng(alpha + beta*pow(x[n],gamma));
    log_lik[n] = bernoulli_logit_lpmf(Y[n] | alpha + beta*pow(x[n],gamma));
  }
  
}

