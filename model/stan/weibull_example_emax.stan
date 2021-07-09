
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=-1,upper=2> cens[N];  // indicates censoring
  int<lower=1> K;  // number of population-level effects
  vector[N] RTS; // RTS
  matrix[N, K] X;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}


parameters {
  // haz(t) = lambda0 * shape * t^(shape-1) * exp(b*ECOG + emax*RTS / (e50 + RTS))
  vector[K] b;  // population-level effects
  real emax;  // maximum effect of RTS
  real<lower=0> e50; // E50
  real<lower=0> lambda0;  // 
  real<lower=0> shape;  // shape parameter
}

transformed parameters {
    vector[N] covariate_effects =  X * b;
    vector[N] hazard ;
    vector[N] cumulative_hazard ;
    
    for (n in 1:N) {
      hazard[n] = lambda0 * shape * pow(Y[n], shape-1) * exp(emax*RTS[n] / (e50+RTS[n]) + covariate_effects[n]);
      cumulative_hazard[n] = lambda0 * pow(Y[n], shape) * exp(emax*RTS[n] / (e50+RTS[n]) + covariate_effects[n]);
    }      
}

model {
  
  // log-likelihood contributions
    for (n in 1:N) {

      // target = log posterior.  We're adding the log-likelihood contributions.
      if (cens[n] == 0) { // Event
        target += log(hazard[n]) - cumulative_hazard[n];
      } else if (cens[n] == 1) { // Right censoring
        target += -cumulative_hazard[n];
      } 

    }

  // priors including all constants
  b ~ normal(0,3); 
  emax ~ normal(0,2);
  e50 ~ lognormal(0,.5);
  lambda0 ~ lognormal(0,5);
  shape ~ lognormal(0, 3);
}

generated quantities {
  vector[N] Ysim;  // response variable
  real u; 
  vector[N] log_lik;

  for (n in 1:N) {
    u = uniform_rng(0,1);
    Ysim[n] = pow(-log(u)/(lambda0 * exp(emax*RTS[n] / (e50+RTS[n]) + covariate_effects[n])), 1/shape);
   
   if (cens[n] == 0) { // Event
        log_lik[n] = log(hazard[n]) - cumulative_hazard[n];
      } else if (cens[n] == 1) { // Right censoring
        log_lik[n] = -cumulative_hazard[n];
      } 
    
  }
  
}

