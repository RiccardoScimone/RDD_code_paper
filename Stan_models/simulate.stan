
functions {
#include /ns_utilities.stan
#include /matern_ns_corr.stan
}


data {
  int<lower=1> N;
  matrix[2,N] Locations;
  vector[N] lambda_1;
   vector[N]lambda_2;
   vector[N] theta;
   vector[N] sigma;
   vector[N] tau;
  int nu;
}

transformed data {
  matrix[2, 2] Aniso[N] =  Compute_Aniso(lambda_1,lambda_2,theta);
  vector[N] dets = Compute_Determinants(lambda_1,lambda_2);
  matrix[N,N] Cov = matern_ns_corr_mat(Locations, Aniso, dets, nu,sigma) + diag_matrix(tau);  
  matrix[N, N] L_cov = cholesky_decompose(Cov);
}

parameters {}
model {}

generated quantities {      
  vector[N] f = multi_normal_cholesky_rng(rep_vector(0, N), L_cov);
}