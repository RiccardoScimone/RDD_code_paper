functions{
#include nngp_matern52_lpdf.stan
}
data {
  int<lower=1> N_obs;
  vector[2] x_obs[N_obs];
  vector[N_obs] y_obs;
  int M;
  int NN_ind[N_obs-1,M];
}

parameters {
  real l_rho;
  real<lower=0> sigma;
  real<lower=0> tau;
  vector[N_obs] latent_functional;
}

transformed parameters{
  real rho = exp(l_rho);
}

model {
  latent_functional ~ nngp_matern52(rep_vector(0,N_obs),x_obs,NN_ind,N_obs,M,sigma,rho);
  y_obs ~ normal(latent_functional, tau);
  sigma ~ std_normal();
  l_rho ~ normal(-3,1);
  tau ~ std_normal();
}







