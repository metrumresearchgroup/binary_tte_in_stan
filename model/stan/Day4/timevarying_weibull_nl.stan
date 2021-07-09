
functions {
  
  real weibull_hazard(real x,          // Function argument
                    real xc,         // Complement of function argument on the domain (defined later)
                    real[] theta,    // parameters
                    real[] x_r,      // data (real)
                    int[] x_i) {
                      
                      real log_lambda = theta[1];
                      real log_shape = theta[2];
                      real rts_effect = theta[3];
                      real beta = theta[4];

                      real ts_rate = x_r[1];
                      real ts_base = x_r[2];
                      real ts_start = x_r[3];
                      real newl_time = x_r[4];
                      real time_start = x_r[5];
                      real bts_mean = x_r[6];

                      real tum_size = ts_start * exp( ts_rate* (x-time_start) );
                      real rel_tum_size = tum_size / ts_base;
                      real log_hazard;
                      
                      log_hazard = log_shape - log_lambda + (exp(log_shape)-1)*(log(x)-log_lambda);
                      log_hazard += rts_effect * log(rel_tum_size);
                      log_hazard += beta * log(ts_base/bts_mean);
                      
                      return( exp(log_hazard) ) ; 
                      
                    }
                    
    real sim_time_rng(real[] time_grid, real[] H_grid  ) {
      
      // Simulate random exponential (equivalent to -log(uniform(0,1)))
      real u = exponential_rng(1);
      real sim_time;
      real delta_time;
      real delta_H;
      
      // Length of time_grid and H_grid assumed to match
      // Assume first element of time_grid is 0
      
      int p = num_elements(time_grid);
      
      int index_lower;
      
      index_lower = 1;
      
      // Find index for interval in which u falls
      while(index_lower <= p) {
        int stop = 0;
        if ((index_lower < p) && (u > H_grid[index_lower]) && (u < H_grid[index_lower+1])) stop=1;
        if (index_lower == p)  stop=1;
        if (stop==1) break;
        index_lower += 1;
        }
        
      // Inverse function in that interval
      if (index_lower < p) {
        delta_time = time_grid[index_lower+1] - time_grid[index_lower];
        delta_H = H_grid[index_lower+1] - H_grid[index_lower];
        sim_time = time_grid[index_lower] + (u-H_grid[index_lower])*delta_time/delta_H;
      } else {
        sim_time = H_grid[p];
      }
      
      return(sim_time);
      
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

transformed data {
  vector[N] ts_rate;
  vector[N] ts_start;
  
  real bts_mean; 
  
  for (n in 1:N) {
    ts_rate[n] = ( KG[n]/1000  - (KD0[n]/1000 * AUC0[n] + KD1[n]/100 * AUC1[n]) );
  }
  
  for (i in 1:Nsubj){
    ts_start[row_first[i]] = TSBASE[row_first[i]];
    ts_start[row_last[i]] = TSBASE[row_first[i]] * exp(ts_rate[row_first[i]]*last_obs_time[i]);
  }

  bts_mean = sum(TSBASE)/N;
 
}

parameters {
  // haz(t) = lambda0 * shape * t^(shape-1) * exp(slope*log(RTS(t)) + beta * baselineTS)
  real beta1;  // effect of RTS
  real beta2;  // effect of baseline tumor size
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
      {log_lambda0,log_shape,beta1,beta2}, 
      {ts_rate[n],TSBASE[n],ts_start[n],NLTIME[n],PREVTIME[n], bts_mean},{NLIND[n]}  );
      
      // Integrate hazard function
      cumulative_hazard[n] = integrate_1d(weibull_hazard, PREVTIME[n], TIME[n], 
      {log_lambda0,log_shape,beta1,beta2}, 
      {ts_rate[n],TSBASE[n],ts_start[n],NLTIME[n],PREVTIME[n], bts_mean},{NLIND[n]}, 1E-08 ) ;
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
  beta1 ~ normal(0,2);
  beta2 ~ normal(0,2);
  log_lambda0 ~ normal(0,2);
  log_shape ~ normal(0, 2);
}



generated quantities {
  vector[Nsubj] log_like;
  vector[N] hazard;
  vector[N] cumulative_hazard;
  vector[Nsubj] hazard_subject;
  vector[Nsubj] cumulative_hazard_subject;
  
  for (n in 1:N) {
    hazard[n] = weibull_hazard(TIME[n], 0.0, 
    {log_lambda0,log_shape,beta1,beta2}, 
    {ts_rate[n],TSBASE[n],ts_start[n],NLTIME[n],PREVTIME[n], bts_mean},{NLIND[n]}  );
    
    // Integrate hazard function
    cumulative_hazard[n] = integrate_1d(weibull_hazard, PREVTIME[n], TIME[n], 
    {log_lambda0,log_shape,beta1,beta2}, 
    {ts_rate[n],TSBASE[n],ts_start[n],NLTIME[n],PREVTIME[n], bts_mean},{NLIND[n]}, 1E-08 ) ;
  }
  
  
  // Calulate likelihood contributions
  for (i in 1:Nsubj) {
    hazard_subject[i] = hazard[row_last[i]];
    cumulative_hazard_subject[i] = sum(cumulative_hazard[row_first[i]:row_last[i]]);
    
    // log-likelihood including all constants
    // Log-likelihood contributions in target +=
    log_like[i] = event[i] * log(hazard_subject[i]) - cumulative_hazard_subject[i];
  }
  
  {
    real maxtime = max(event_time);
    int grid_size = 100;
    real time_grid[grid_size+1];
    real H_grid[grid_size+1];
    
    time_grid[1] = 0.0;
    for (j in 1:grid_size) time_grid[j+1] = maxtime * j/grid_size;
    
  // Simulate event times from model
  for (i in 1:Nsubj) {
    H_grid[1] = 0.0;
    for (j in 1:grid_size) {}
    
  }
  }
}
