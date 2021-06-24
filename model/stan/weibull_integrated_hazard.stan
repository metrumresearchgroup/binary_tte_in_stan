// generated with brms 2.14.4
functions {
  
  real weibull_hazard(real x,          // Function argument
                    real xc,         // Complement of function argument
                                     //  on the domain (defined later)
                    real[] theta,    // parameters
                    real[] x_r,      // data (real)
                    int[] x_i) {
                      
                      real llambda = theta[1];
                      real lshape = theta[2];
                      
                      return( lshape * llambda * pow(x, lshape-1) );
                      
                    }
}

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
  // haz(t) = lambda0 * shape * t^(shape-1) * exp(eta)
  vector[K] b;  // population-level effects
  real slope;  // maximum effect of RTS
  real<lower=0> lambda0;  // 
  real<lower=0> shape;  // shape parameter
}


model {

    vector[N] mu ;
    vector[N] lambda ;
    vector[N] hazard;
    vector[N] cumulative_hazard;
    
    for (n in 1:N) {
      mu[n]  =  (slope * RTS[n])  + (X * b)[n];
      lambda[n] = lambda0 * exp(mu[n]);

      // Evaluate hazard function
      hazard[n] = weibull_hazard(Y[n], 0.0, {lambda[n],shape}, {0.0}, {0}  );
      
      // Integrate hazard function
      cumulative_hazard[n] = integrate_1d(weibull_hazard, 1E-06, Y[n], {lambda[n],shape}, {0.0}, {0}, 1E-08 ) ;
    
    }

  // likelihood including all constants
  if (!prior_only) {
    for (n in 1:N) {
    // special treatment of censored data
      if (cens[n] == 0) {
        // Log-likelihood contributions in target +=
        target += log(hazard[n]) - cumulative_hazard[n];
      } else if (cens[n] == 1) {
        target += -cumulative_hazard[n];
      } 
    }
  }
  // priors including all constants
  b ~ normal(0,3); 
  slope ~ normal(0,2);
  lambda0 ~ lognormal(0,5);
  shape ~ lognormal(0, 3);
}

generated quantities {
  vector[N] Ysim;  // response variable
  real u; 
  vector[N] log_lik;
  vector[N] mu ;
  vector[N] lambda ;
  vector[N] hazard;
  vector[N] cumulative_hazard;
  
  for (n in 1:N) {
    mu[n]  =  (slope * RTS[n])  + (X * b)[n];
    lambda[n] = lambda0 * exp(mu[n]);
    
    // Evaluate hazard function
    hazard[n] = weibull_hazard(Y[n], 0.0, {lambda[n],shape}, {0.0}, {0}  );
    
    // Integrate hazard function
    cumulative_hazard[n] = integrate_1d(weibull_hazard, 1E-06, Y[n], {lambda[n],shape}, {0.0}, {0}, 1E-08 ) ;
    
  }
  
  for (n in 1:N) {
    u = uniform_rng(0,1);
    Ysim[n] = pow(-log(u)/lambda[n], 1/shape);
    
    if (cens[n] == 0) { // Event
    log_lik[n] = log(hazard[n]) - cumulative_hazard[n];
    } else if (cens[n] == 1) { // Right censoring
    log_lik[n] = -cumulative_hazard[n];
    } 
    
  }
  
}

