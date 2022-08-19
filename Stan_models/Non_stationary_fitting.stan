functions {
#include matern_ns_corr.stan
#include ns_utilities.stan
}
     data {
      int<lower=0,upper=1> fit_mu;
      int<lower=0,upper=1> fit_sigma;
      int<lower=0,upper=1> fit_lambda;
      int<lower=0, upper=1> fit_theta;
      int<lower=1> N;
      int<lower=0> P_mu;
      int<lower=0> P_sigma;
      int<lower=0> P_lambda;
      int<lower=0> P_theta;
      int<lower=1> nu;
      vector[N] Y;
      matrix[N, P_mu] X_mu[fit_mu];
      matrix[N, P_sigma] X_sigma[fit_sigma];
      matrix[N, P_lambda] X_lambda[fit_lambda];
      matrix[N, P_theta] X_theta[fit_theta];
      matrix[2,N] Locations;
      vector[P_mu] uB_mu[fit_mu];
      vector[P_sigma] uB_sigma[fit_sigma];
      vector[P_lambda] uB_lambda[fit_lambda];
      vector[P_theta] uB_theta[fit_theta];
      matrix[P_mu, P_mu] VB_mu[fit_mu];
      matrix[P_sigma, P_sigma] VB_sigma[fit_sigma];
      matrix[P_lambda, P_lambda] VB_lambda[fit_lambda];
      matrix[P_theta, P_theta] VB_theta[fit_theta];

  }

transformed data{
      matrix[P_mu, P_mu] LB_mu;
      matrix[P_sigma, P_sigma] LB_sigma;
      matrix[P_lambda, P_lambda] LB_lambda;
      matrix[P_theta, P_theta] LB_theta;
      if(fit_mu) LB_mu = cholesky_decompose(VB_mu[1]);
      if(fit_sigma) LB_sigma = cholesky_decompose(VB_sigma[1]);
      if(fit_lambda) LB_lambda = cholesky_decompose(VB_lambda[1]);
      if(fit_theta) LB_theta = cholesky_decompose(VB_theta[1]);
}

parameters {
  vector[P_mu] alfa_beta_mu[fit_mu];
  vector[P_sigma] alfa_beta_sigma[fit_sigma];
  vector[P_lambda] alfa_beta_lambda_1[fit_lambda];
  vector[P_lambda] alfa_beta_lambda_2[fit_lambda];
  vector[P_theta] alfa_beta_theta[fit_theta];
}

transformed parameters{
  vector[P_mu] beta_mu[fit_mu]; 
  vector[P_sigma] beta_sigma[fit_sigma];
  vector[P_lambda] beta_lambda_1[fit_lambda];
  vector[P_lambda] beta_lambda_2[fit_lambda]; 
  vector[P_theta] beta_theta[fit_theta]; 
  if(fit_mu)
    beta_mu[1]= uB_mu[1] + LB_mu*alfa_beta_mu[1];
  if(fit_sigma)
    beta_sigma[1] = uB_sigma[1] + LB_sigma*alfa_beta_sigma[1];
  if(fit_lambda){
    beta_lambda_1[1] = uB_lambda[1] + LB_lambda*alfa_beta_lambda_1[1];
    beta_lambda_2[1] = uB_lambda[1] + LB_lambda*alfa_beta_lambda_2[1];
  }
  if(fit_theta)
    beta_theta[1] = uB_theta[1] + LB_theta*alfa_beta_theta[1];
  
}

model {
  matrix[2,2] Aniso[N];
  vector[N] Dets;
  matrix[N,N] cho_cov;
  if (fit_mu)
    alfa_beta_mu[1] ~ std_normal();
  if (fit_sigma)
    alfa_beta_sigma[1] ~ std_normal();
  if (fit_lambda){
    alfa_beta_lambda_1[1] ~ std_normal();
    alfa_beta_lambda_2[1] ~ std_normal();
  }
  if (fit_theta)
    alfa_beta_theta[1] ~ std_normal();
  Aniso = Compute_Aniso(exp(X_lambda[1]*beta_lambda_1[1]),exp(X_lambda[1]*beta_lambda_2[1]),pi()/2*inv_logit(X_theta[1]*beta_theta[1]));
  Dets = Compute_Determinants(exp(X_lambda[1]*beta_lambda_1[1]),exp(X_lambda[1]*beta_lambda_2[1]));
  if (fit_sigma) cho_cov = cholesky_decompose(matern_ns_corr_mat(Locations,Aniso,Dets,nu, exp(X_sigma[1]*beta_sigma[1])));
  else           cho_cov = cholesky_decompose(matern_ns_corr_mat(Locations,Aniso,Dets,nu, rep_vector(1.,N)));
  if (fit_mu)
    Y ~ multi_normal_cholesky(X_mu[1]*beta_mu[1],cho_cov);
  else
    Y ~ multi_normal_cholesky(rep_vector(0,N),cho_cov);
}

