// In this model, I let the coefficient of zpSR to depends on Only Omnivory but the influences of omnivory do not depend on the Omernik level 3 ecoregion. Other coefficients are region-dependent
// In addition, I used the "generated quantitie" chunk to calculate the likelihood for the final model, so that wAIC or loo can be calculated using the "loo" package. 

data {
  int<lower=1> N; // the number of observations
  int<lower=1> K1; // number of column of the model matrix (6)
  matrix[N, K1] X; // N (# of obs) by K (# of predictors) model matrix 
  int<lower=1> J; // the number of Omernik level 3 ecoregion
  int<lower=1,upper=J> region[N];
  vector[N] phyDen; // response variable, Y
  
  vector[N] zpSR; // 
  int<lower=1> K2; // number of column of the model matrix determining the effects of zpSR
  matrix[N, K2] X2; // the model matrix determining the effects of zpSR
  
  
  int<lower=1> newN; // the number of observations
  int<lower=1> newK1; // number of column of the model matrix (6)
  matrix[newN, newK1] newX; // N (# of obs) by K (# of predictors) model matrix 
  int<lower=1> newJ; // the number of Omernik level 3 ecoregion
  int<lower=1,upper=newJ> newregion[newN];
  
  vector[newN] newzpSR; // 
  int<lower=1> newK2; // number of column of the model matrix determining the effects of zpSR
  matrix[newN, newK2] newX2; // the model matrix determining the effects of zpSR
  
  
  int<lower=1> simN; // the number of observations
  int<lower=1> simK1; // number of column of the model matrix (6)
  matrix[simN, simK1] simX; // N (# of obs) by K (# of predictors) model matrix 
  int<lower=1> simJ; // the number of Omernik level 3 ecoregion
  int<lower=1,upper=simJ> simregion[simN];
  
  vector[simN] simzpSR; // 
  int<lower=1> simK2; // number of column of the model matrix determining the effects of zpSR
  matrix[simN, simK2] simX2; // the model matrix determining the effects of zpSR
}
parameters {
  vector[K1] gamma;
  vector<lower=0>[K1] tau;
  vector[K1] beta_raw[J]; //Applying Matt's trick
  
  vector[N] beta_zpSR;
  vector[K2] gamma_beta_zpSR;
  vector<lower=0>[K2] tau_beta_zpSR;
  vector[K2] b_raw[N]; //Applying Matt's trick

  real<lower=0> sigma_phyDen; //standard deviation of the individual observations
  real<lower=0> sigma_beta_zpSR; //standard deviation of the individual observations
}
transformed parameters{
  vector[K1] beta[J];
  vector[K2] b[N];
  vector[N] mu_phyD;
  vector[N] mu_beta_zpSR;
  
  for(j in 1:J){ // Region-dependent coefficients for coefficients other than ZpSR
    beta[j] = gamma + tau .* beta_raw[j]; 
  }

  for(i in 1:N){ // Region-independent coefficients for coefficient of ZpSR
    b[i] = gamma_beta_zpSR + tau_beta_zpSR .* b_raw[i];
  }
  
  for (i in 1:N){
    mu_beta_zpSR[i] = X2[i] * b[i];
    mu_phyD[i] = X[i] * beta[region[i]] + zpSR[i] * beta_zpSR[i];
  }
}

model {
  //priors
  gamma ~ normal(0, 10); //weakly informative priors on the regression coefficients
  tau ~ cauchy(0, 5); //weakly informative priors
  gamma_beta_zpSR ~ normal(0, 10); //weakly informative priors on the regression coefficients
  tau_beta_zpSR ~ cauchy(0, 5); //weakly informative priors
  sigma_beta_zpSR ~ gamma(2, 0.5); //weakly informative priors
  sigma_phyDen ~ gamma(2, 0.5); //weakly informative priors
  
  //likelihood
  for (i in 1:N){
    b_raw[i] ~ normal(0, 1);  
  }
  for (j in 1:J){
    beta_raw[j] ~ normal(0, 1);
  }
  
  beta_zpSR ~ normal(mu_beta_zpSR, sigma_beta_zpSR);
  phyDen ~ normal(mu_phyD, sigma_phyDen);
} 

generated quantities {
  vector[N] log_lik;
  vector[N] y_pred;
  vector[newN] y_12to07;
  vector[simN] y_sim;
  
  for (i in 1:N){
    log_lik[i] = normal_lpdf(phyDen[i] | mu_phyD[i], sigma_phyDen);
    y_pred[i] = X[i] * beta[region[i]] + zpSR[i] * beta_zpSR[i];
  }
  
  for (i in 1:newN){
    y_12to07[i] = newX[i] * beta[newregion[i]] + newzpSR[i] * beta_zpSR[i];
  }
   
  for (i in 1:simN){
    y_sim[i] = simX[i] * beta[simregion[i]] + simzpSR[i] * beta_zpSR[i];
  }
}

