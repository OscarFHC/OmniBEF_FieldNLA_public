// This is Hierarchical model stan script
// In addition, I used the "generated quantitie" chunk to calculate the likelihood for the final model, so that wAIC or loo can be calculated using the "loo" package. 

data {
  int<lower=1> N; // the number of observations
  int<lower=1> K_Gzp; // number of column of the model matrix (6)
  matrix[N, K_Gzp] X_Gzp; // N (# of obs) by K1 (# of predictors; 6) model matrix 
  vector[N] Gzp; // response variable, Y
  
  vector[N] intn; // 
  int<lower=1> K_bintn; // number of column of the model matrix determining the interaction effects between zpSR and omni
  matrix[N, K_bintn] X_bintn; // the model matrix determining the effects of the interaction effects between zpSR and omni
}
parameters {
  vector[K_Gzp] gamma_Gzp;
  vector<lower=0>[K_Gzp] tau_Gzp;
  vector[K_Gzp] beta_raw_Gzp[N];
  
  vector[N] beta_intn;
  vector[K_bintn] gamma_beta_intn;
  vector<lower=0>[K_bintn] tau_beta_intn;
  vector[K_bintn] b_raw[N];
  
  real<lower=0> sigma_Gzp; //standard deviation of the individual observations
  real<lower=0> sigma_beta_intn; //standard deviation of the individual observations
}
transformed parameters{
  vector[K_Gzp] beta_Gzp[N]; //applying Matt's trick
  vector[K_bintn] b_bintn[N];  //applying Matt's trick
  
  vector[N] mu_Gzp;
  vector[N] mu_beta_intn;
  
  for(i in 1:N){
    beta_Gzp[i] = gamma_Gzp + tau_Gzp .* beta_raw_Gzp[i];
    b_bintn[i] = gamma_beta_intn + tau_beta_intn .* b_raw[i];
  }
    
  for (i in 1:N){
    mu_beta_intn[i] = X_bintn[i] * b_bintn[i];
    mu_Gzp[i] = X_Gzp[i] * beta_Gzp[i] + intn[i] * beta_intn[i];
  } 
}
model {
  //priors
  gamma_Gzp ~ normal(0, 10); //weakly informative priors on the regression coefficients
  tau_Gzp ~ normal(0, 10); //weakly informative priors
  gamma_beta_intn ~ normal(0, 10); //weakly informative priors on the regression coefficients
  tau_beta_intn ~ normal(0, 10); //weakly informative priors
  
  sigma_beta_intn ~ gamma(2, 0.5); //weakly informative priors
  sigma_Gzp ~ gamma(2, 0.5); //weakly informative priors
  
  //likelihood
  for (i in 1:N){
    beta_raw_Gzp[i] ~ normal(0, 1); 
    b_raw[i] ~ normal(0, 1);
  }
  
  beta_intn ~ normal(mu_beta_intn, sigma_beta_intn);  
  Gzp ~ normal(mu_Gzp, sigma_Gzp);
} 
generated quantities {
  vector[N] log_lik_Gzp;
  
  vector[N] Gzp_pred;
  
  for (i in 1:N){
    log_lik_Gzp[i] = normal_lpdf(Gzp[i] | mu_Gzp[i], sigma_Gzp);
    
    Gzp_pred[i] = X_Gzp[i] * beta_Gzp[i] + intn[i] * beta_intn[i];
  }
}
