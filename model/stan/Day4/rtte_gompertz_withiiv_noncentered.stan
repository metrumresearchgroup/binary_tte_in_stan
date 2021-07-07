
data {
  int<lower=0> Nsubj;
  int<lower=0> Nobs;
  int delta[Nobs];
  real<lower=0> time[Nobs];
  real<lower=0> prev_time[Nobs];
  int ID[Nobs];
  vector[Nobs] CAVGSS;
}

transformed data {
  real<lower=0> day[Nobs];
  real<lower=0> prev_day[Nobs];

  for (n in 1:Nobs){
    day[n] = time[n]/90;
    prev_day[n] = prev_time[n]/90;
  }
}

parameters {
  real log_alpha_pop;
  real<lower=0> omega;
  real log_gamma;
  real eta[Nsubj];
}

transformed parameters {

  real gamma = exp(log_gamma);
  vector[Nobs] cumulative_hazard;
  
  for (n in 1:Nobs) {
    cumulative_hazard[n] = exp(log_alpha_pop + omega*eta[ID[n]]) / gamma * (exp(-gamma*prev_time[n]) - exp(-gamma*time[n]));
  }
}


model {
  //vector[Nobs] cumulative_hazard;

  log_alpha_pop ~ normal(0, 2);
  eta ~ normal(0,1);
  log_gamma ~ normal(0,2);
  omega ~ normal(0,2);
  
   for (n in 1:Nobs) {
  //    cumulative_hazard[n] = exp(log_alpha_pop + omega*eta[ID[n]]) / (gamma/90) * (exp(-gamma*prev_day[n]) - exp(-gamma*day[n]));
      target += delta[n] * (log_alpha_pop + omega*eta[ID[n]] - gamma*day[n]) - cumulative_hazard[n];
   }
}

generated quantities {
  vector[Nsubj] log_like = rep_vector(0.0, Nsubj);
  
  for (n in 1:Nobs) {
    log_like[ID[n]] += delta[n] * (log_alpha_pop + omega*eta[ID[n]] - gamma*day[n]) - cumulative_hazard[n];
  }
}

