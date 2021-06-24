// generated with brms 2.14.4
functions {
  
  real[] ode(real t,        // time
           real[] y,      // state
           real[] theta,  // parameters
           real[] x_r,    // data (real)
           int[] x_i) {   // data (integer)
           
           real sigma = theta[1];  // scale for baseline hazard
           real alpha = theta[2];  // shape for baseline hazard
           real beta = theta[3];  // PH effect of rts on hazard
           
           real KG = x_r[1];
           real KD0 = x_r[2] ;
           real KD1 = x_r[3];
           real IBASE = x_r[4];
           real E0 = x_r[5] ;
           real E1 = x_r[6];
           
           real rts;
           
           real dydt[2];
           
           real baseline_hazard = (alpha/sigma) * pow(t/sigma, alpha-1);
           
           dydt[1] = ( KG/1000  - (KD0/1000 * E0 + KD1/100 * E1) ) * y[1] ;
           
           rts = y[1] / (IBASE * 1000); 

           dydt[2] = baseline_hazard * exp(beta * rts);
           
           return dydt;
           
           
                      
                    }
}

data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=-1,upper=2> cens[N];  // indicates censoring
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  print(Xc);

}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> shape;  // shape parameter
}
transformed parameters {
}
model {
  
  real hazard[N];
  real cumulative_hazard[N];
  
  // likelihood including all constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b;
    for (n in 1:N) {
      // apply the inverse link function
      mu[n] = exp(mu[n]) / tgamma(1 + 1 / shape);
    }
    for (n in 1:N) {
      // Evaluate hazard function
      //hazard[n] = weibull_hazard(Y[n], 0.0, {mu[n],shape}, {0.0}, {0}  );
      
      // Integrate hazard function
      cumulative_hazard[n] = integrate_ode_rk45(ode, {}, t0, ts, theta, x_r, x_i);
      
    // special treatment of censored data
      if (cens[n] == 0) {
        // Log-likelihood contributions in target +=
        target += log(hazard[n]) - cumulative_hazard[n];
      } else if (cens[n] == 1) {
        target += cumulative_hazard[n];
      } 
    }
  }
  // priors including all constants
  b ~ normal(0, 3);
  Intercept ~ student_t(3, 5.9, 2.5);
  shape ~ lognormal(0, 3);
}

generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}

