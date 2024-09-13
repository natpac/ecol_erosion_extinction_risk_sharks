data {
  int<lower=0> T1;       // Number of time points
  int<lower=0> T2;       // Number of time points
  int<lower=0> T_pred;      // Number of future time points to predict
  vector[T1] y1;        // First observed dataset
  vector[T2] y2;        // Second observed dataset
}

parameters {
  real alpha1;            // Drift parameter
  real<lower=0> sigma1;   // Volatility parameter
  vector[T1] mu1;          // Random walk mean
  real alpha2;            // Drift parameter
  real<lower=0> sigma2;   // Volatility parameter
  vector[T2+T_pred] mu2;          // Random walk mean
}

model {
  // Prior distributions
  alpha1 ~ normal(0, 100000);    // Prior on drift
  sigma1 ~ cauchy(0, 500);    // Prior on volatility
  mu1[1] ~ normal(y1[1], 100000);   // Prior on initial value
  alpha2 ~ normal(0, 2000000);    // Prior on drift
  sigma2 ~ cauchy(0, 500);    // Prior on volatility
  mu2[1] ~ normal(y2[1], 1000000);   // Prior on initial value
  
  // Random walk model
  for (t in 2:T1) {
    mu1[t] ~ normal(mu1[t - 1] + alpha1, sigma1);  // Random walk with drift
  }
  for (t in 2:(T2+T_pred)) {
    mu2[t] ~ normal(mu2[t - 1] + alpha2, sigma2);  // Random walk with drift
  }
  
  // Likelihood
  y1 ~ normal(mu1, sigma1);  // Likelihood of the observed1 data
  y2 ~ normal(mu2[1:T2], sigma2);  // Likelihood of the observed2 data
}

generated quantities {
  vector[T1] y1_pred;  // Predicted values for the observed data
  vector[T2+T_pred] y2_pred;// Predicted values for the observed data
  vector[T1] CPUE;
  
  // Generate predicted values for the observed data
  for (t in 1:T1) {
    y1_pred[t] = normal_rng(mu1[t], sigma1);
  }
  
    
  // vector[T_pred] y2_future;  // Predictions for 9 additional observations
  
  // Generate predicted values for the observed data
  for (t in 1:T2+T_pred) {
    y2_pred[t] = normal_rng(mu2[t], sigma2);
  }
  
  // y2_future[1] = normal_rng(mu2[T2] + alpha2, sigma2);
  // // Generate predictions for 9 additional observations
  // for (t in 2:9) {
  //   y2_future[t] = normal_rng(y2_future[t - 1] + alpha2, sigma2);
  // }
  
  
  // Calculate CPUE as the ratio of y1 over y2
  for (t in 1:T1) {
    CPUE[t] = mu1[t] / mu2[t];
  }
}
