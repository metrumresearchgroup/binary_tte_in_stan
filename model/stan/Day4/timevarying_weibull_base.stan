
functions {
  
  real weibull_hazard(real x,          // Function argument
                    real xc,         // Complement of function argument on the domain (defined later)
                    real[] theta,    // parameters
                    real[] x_r,      // data (real)
                    int[] x_i) {
                      
                      real log_lambda = theta[1];
                      real log_shape = theta[2];
                      
                      real log_hazard;

                      log_hazard = log_shape - log_lambda + (exp(log_shape)-1)*(log(x)-log_lambda);

                      return( exp(log_hazard) ) ; 
                      
                    }
}

data {
  int<lower=1> N;  // total number of observations
  int<lower=1> Nsubj; // number of subjects
  
  // Data for calculating time-varying hazard
  int<lower=1> ID[N]; // subject ID
  vector[N] TIME;  // response variable
  vector[N] PREVTIME;  // for integrating
  vector[N] AUC0;  // AUC drug 0
  vector[N] AUC1;  // AUC drug 1
  vector[N] KG;  // tumor growth rate constant
  vector[N] KD0;  // tumor killing effect of drug 0
  vector[N] KD1;  // tumor killing effect of drug 1
  vector[N] TSBASE; // baseline tumor size
  vector[N] NLTIME; // new lesion time
  int NLIND[N]; // indicator for new lesion
  
  // Data for calculating likelihood
  int<lower=1> IDsubj[Nsubj]; // subject ID
  int<lower=0,upper=1> event[Nsubj];  // indicates censoring
  vector[Nsubj] event_time; // death or censoring time
  int<lower=1> row_first[Nsubj];
  int<lower=1> row_last[Nsubj];
  vector[Nsubj] last_obs_time;
}


parameters {
  // haz(t) = shape/lambda0 * (t/lambda0)^(shape-1)
  real log_lambda0;  // baseline scale parameter
  real log_shape;  // shape parameter
}


model {

    vector[N] hazard;
    vector[N] cumulative_hazard;
    vector[Nsubj] hazard_subject;
    vector[Nsubj] cumulative_hazard_subject;
    
    for (n in 1:N) {
      // Evaluate hazard function
      hazard[n] = weibull_hazard(TIME[n], 0.0, 
      {log_lambda0,log_shape}, 
      {0.0},{0}  );
      
      // Integrate hazard function
      cumulative_hazard[n] = integrate_1d(weibull_hazard, PREVTIME[n], TIME[n], 
      {log_lambda0,log_shape}, 
      {0.0},{0}, 1E-08 ) ;
    }


    for (i in 1:Nsubj) {
      hazard_subject[i] = hazard[row_last[i]];
      cumulative_hazard_subject[i] = sum(cumulative_hazard[row_first[i]:row_last[i]]);
    }

  // log-likelihood including all constants
  for (i in 1:Nsubj) {
        // Log-likelihood contributions in target +=
        target += event[i] * log(hazard_subject[i]) - cumulative_hazard_subject[i];
  }
  
  // priors including all constants
  log_lambda0 ~ normal(0,2);
  log_shape ~ normal(0, 2);
}


generated quantities {
  vector[Nsubj] log_like;
  vector[N] hazard;
  vector[N] cumulative_hazard;
  vector[Nsubj] hazard_subject;
  vector[Nsubj] cumulative_hazard_subject;
  
  vector[Nsubj] Ysim;
  real u;
  
  for (n in 1:N) {
    // Evaluate hazard function
    hazard[n] = weibull_hazard(TIME[n], 0.0, {log_lambda0,log_shape}, {0.0},{0}  );
    
    // Integrate hazard function
    cumulative_hazard[n] = integrate_1d(weibull_hazard, PREVTIME[n], TIME[n], 
    {log_lambda0,log_shape}, {0.0},{0}, 1E-08 ) ;
  }
  
  
  for (i in 1:Nsubj) {
    hazard_subject[i] = hazard[row_last[i]];
    cumulative_hazard_subject[i] = sum(cumulative_hazard[row_first[i]:row_last[i]]);
    
    // Log-likelihood contributions in target +=
    log_like[i] = event[i] * log(hazard_subject[i]) - cumulative_hazard_subject[i];
    
    u = uniform_rng(0,1);
    Ysim[i] = exp(log_lambda0) * pow(-log(u), 1/exp(log_shape));
    
  }
  
}

