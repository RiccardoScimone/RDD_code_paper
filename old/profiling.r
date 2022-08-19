rm(list = ls())
library(BayesNSGP)
library(tidyverse)
library(rstan)
library(rstanarm)
library(shinystan)
library(cmdstanr)
source("NNMatrix.R")

### Defining spatial grid
sq_N = 40
a = -1
b = 1
grid = expand.grid(seq(a,b,length.out = sq_N),seq(a,b,length.out = sq_N))
names(grid) = c("x","y")
k = 5
distances = NNMatrix(coords = grid,n.neighbors = k,n.omp.threads = 6)

theta = (atan(grid$y/grid$x))*as.numeric(grid$x*grid$y >= 0) + (pi/2 + atan(grid$y/grid$x)) * as.numeric(grid$x*grid$y < 0)
range(theta)
lambda1 = 0.5 * as.numeric(grid$x*grid$y >= 0) + 0.3 * as.numeric(grid$x*grid$y < 0)
lambda2 = 0.5 * as.numeric(grid$x*grid$y < 0) + 0.3 * as.numeric(grid$x*grid$y >= 0)
res_theta = logit(2/pi * theta)
log_1 = log(lambda1)
log_2 = log(lambda2)
Sigma11_vec = inverseEigen(log_1,log_2,res_theta,1)
Sigma22_vec = inverseEigen(log_1,log_2,res_theta,2)
Sigma12_vec = inverseEigen(log_1,log_2,res_theta,3)

sigma = 1
tau = 1
nu = 3
dist_list = nsDist(grid)
Cor_mat <- nsCorr( dist1_sq = dist_list$dist1_sq, dist2_sq = dist_list$dist2_sq,
                   dist12 = dist_list$dist12, Sigma11 = Sigma11_vec,
                   Sigma22 = Sigma22_vec, Sigma12 = Sigma12_vec, nu = nu )
cholesky = chol(Cor_mat)
data = t(cholesky) %*% rnorm(sq_N^2)
grid$process = data

grid = grid[distances$ord,]

P = 1                  # number of regression coefficients
P_lambdas = 2
P_theta = 1
#uB = rep(0, P + 1)     # mean vector in the Gaussian prior of beta
VB_beta = diag(P + 1)*10  # covariance matrix in the Gaussian prior of beta
VB_lambdas = diag(P_lambdas + 1)*10
VB_theta =  diag(P_theta + 1)*10
#ss = 3 * sqrt(2)       # scale parameter in the normal prior of sigma 
#st = 3 * sqrt(0.1)     # scale parameter in the normal prior of tau     

### data
NN_ind = distances$NN_ind
N = sq_N^2
M = k
Locations = cbind(grid$x,grid$y)

Y = grid$process[1:N,]
X = matrix(data = 1, nrow = N, ncol = 1)
X = cbind(X,grid$x^2+grid$y^2)
X_lambdas = cbind(rep(1,N),as.numeric(grid$x*grid$y >= 0),as.numeric(grid$x*grid$y < 0))
X_theta = cbind(rep(1,N),res_theta)


model <- cmdstan_model("Profiling.stan")

data <- list(N = N, M = M, P = P,
             P_lambdas = P_lambdas, P_theta = P_theta,
             Y = Y, X = X, X_lambdas = X_lambdas,X_theta = X_theta, Locations = Locations,
             NN_ind = NN_ind, VB_beta = VB_beta, VB_lambdas = VB_lambdas, VB_theta = VB_theta, nu = nu ) 

parameters <- c("beta", "beta_lambda1", "beta_lambda2", "beta_theta")

fit = model$sample(data = data,chains = 1,iter_warmup = 300,iter_sampling = 300, adapt_delta = 0.999, max_treedepth = 6, refresh = 1,seed = 123)
save(fit, file = "Profiled.Rdata")
fit$cmdstan_diagnose()
fit$cmdstan_summary()
fit$inv_metric()

launch_shinystan(as.shinystan(fit))


expose_stan_functions(stanc("Profiling.stan"))


aniso_list = list()
aniso_matrix = NULL;
dets = NULL

for ( i in 1:N){
  aniso_list[[i]] = rbind(c(Sigma11_vec[i],Sigma12_vec[i]), c(Sigma12_vec[i],Sigma22_vec[i]))
  dets = c(dets, Sigma11_vec[i]*Sigma22_vec[i] - Sigma12_vec[i]^2)
  aniso_matrix = cbind(aniso_matrix,aniso_list[[i]]);
}



Cor = matrix(0,nrow = N,ncol = N)

system.time( 
  for ( i in 1:N)
    for( j in (i+1):N)
  { 
    #Cor[i,i] = 1
    if ( j < N){
        Cor[i,j] = matern_ns_corr(Locations[i,],Locations[j,],aniso_list[[i]],aniso_list[[j]],dets[i],dets[j],nu)
        Cor[j,i] = Cor[i,j]
    }
  }
)

Cor = matern_ns_corr_mat(t(Locations),aniso_list,dets,nu)
vec = matern_ns_corr_vec(Locations[1,],t(Locations[-1,]),aniso_list[[1]],aniso_list[2:N],dets[1],dets[-1],nu)

system.time(matern_ns_corr_mat(t(Locations),aniso_list,dets,nu))



nngp_lpdf(Y,rep(0,N),1,0,theta/2,lambda1,lambda2,t(Locations),split(NN_ind,1:(N-1)),N,M,nu)


