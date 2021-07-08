//
// Gompertz model:
//   hazard function: h(t) = alpha * exp(-gamma * t)
//   Cumulative hazard function: H(t) = alpha / gamma * (1 - exp(-gamma * t))
//   CH between two observations: H(t2) - H(t1) = alpha / gamma * (exp(-gamma * t1) - exp(-gamma * t2)) (for t1 < t2 )
//

data {
  int<lower=0> Nsubj;
  int<lower=0> Nobs;
  int delta[Nobs];
  real<lower=0> time[Nobs];
  real<lower=0> prev_time[Nobs];
  int ID[Nobs];
  vector[Nobs] CAVGSS;
}

parameters {
  real log_alpha_pop;
  real log_gamma;
}

transformed parameters {
  vector[Nobs] cumulative_hazard;
  real<lower=0> gamma = exp(log_gamma);
  
  for (n in 1:Nobs) {
    cumulative_hazard[n] = exp(log_alpha_pop) / (gamma/90) * (exp(-gamma*(prev_time[n]/90)) - exp(-gamma*(time[n]/90)));
  }  
}


model {
  log_alpha_pop ~ normal(0, 2);
  log_gamma ~ normal(0,2);

  for (n in 1:Nobs) {
    target += delta[n] * (log_alpha_pop - gamma*(time[n]/90)) - cumulative_hazard[n];
  }
}

generated quantities {
  vector[Nsubj] log_like = rep_vector(0.0, Nsubj);
  
  for (n in 1:Nobs) {
    log_like[ID[n]] += delta[n] * (log_alpha_pop - gamma*(time[n]/90)) - cumulative_hazard[n];
  }
}

