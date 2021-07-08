
data {
  int<lower=0> Nsubj;
  int<lower=0> Nobs;
  int delta[Nobs];
  real<lower=0> time[Nobs];
  real<lower=0> prev_time[Nobs];
  int ID[Nobs];
  vector[Nsubj] CAVGSS;
}

parameters {
  real log_alpha_pop;
  real<lower=0> omega;
  real log_gamma;
  real log_alpha[Nsubj];
  real beta;
}

transformed parameters {
  
  vector[Nobs] cumulative_hazard;
  real gamma = exp(log_gamma);
  
  for (n in 1:Nobs) {
    cumulative_hazard[n] = exp(log_alpha[ID[n]]) / gamma * (exp(-gamma*prev_time[n]) - exp(-gamma*time[n]));
  }
}


model {

  log_alpha_pop ~ normal(0, 2);
  log_gamma ~ normal(0,2);
  beta ~ normal(0,2);
  omega ~ normal(0,2);
  
  // Taking advantage of vectorization to more efficently compute a Nsubj dimensional normal
  log_alpha ~ normal(log_alpha_pop + beta*CAVGSS, omega);
  

  for (n in 1:Nobs) {
        target += delta[n] * (log_alpha[ID[n]] - gamma*time[n]/90) - cumulative_hazard[n];
  }
}

generated quantities{
  vector[Nobs] log_like;
  real log_alpha_sim[Nsubj];
  
  for (n in 1:Nobs) {
        log_like[n] = delta[n] * (log_alpha[ID[n]] - gamma*time[n]/90) - cumulative_hazard[n];
  }
  
  for (i in 1:Nsubj){
    log_alpha_sim[i] = normal_rng(log_alpha_pop + beta*CAVGSS[i], omega);
    
  }
}

