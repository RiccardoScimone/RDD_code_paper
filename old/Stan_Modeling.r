rm(list = ls())
library(BayesNSGP)
library(tidyverse)
library(rstan)
library(rstanarm)
library(tictoc)
source("NNMatrix.R")

### Defining spatial grid
sq_N = 80
a = -1
b = 1
grid = expand.grid(seq(a,b,length.out = sq_N),seq(a,b,length.out = sq_N))
names(grid) = c("x","y")
k = 5
distances = NNMatrix(coords = grid,n.neighbors = k,n.omp.threads = 6)



theta = (atan(grid$y/grid$x))*as.numeric(grid$x*grid$y >= 0) + (pi/2 + atan(grid$y/grid$x)) * as.numeric(grid$x*grid$y < 0)
range(theta)
lambda1 = 0.05 * as.numeric(grid$x*grid$y >= 0) + 0.3 * as.numeric(grid$x*grid$y < 0)
lambda2 = 0.05 * as.numeric(grid$x*grid$y < 0) + 0.3 * as.numeric(grid$x*grid$y >= 0)
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
tic()
Cor_mat <- nsCorr( dist1_sq = dist_list$dist1_sq, dist2_sq = dist_list$dist2_sq,
                   dist12 = dist_list$dist12, Sigma11 = Sigma11_vec,
                   Sigma22 = Sigma22_vec, Sigma12 = Sigma12_vec, nu = nu )
toc()
cholesky = chol(Cor_mat)
data = t(cholesky) %*% rnorm(sq_N^2)
grid$process = data

grid = grid[distances$ord,]

library(viridis)
x11()
p = ggplot(grid) + geom_tile(mapping = aes(x,y, fill = process)) + scale_fill_viridis(option = "inferno")
plot(p)


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
X_theta = cbind(rep(1,N),res_theta)[distances$ord,]

#------------------------------ response NNGP ---------------------------------#
options(mc.cores = parallel::detectCores())
data <- list(N = N, M = M, P = P,
             P_lambdas = P_lambdas, P_theta = P_theta,
             Y = Y, X = X, X_lambdas = X_lambdas,X_theta = X_theta, Locations = Locations,
             NN_ind = NN_ind, VB_beta = VB_beta, VB_lambdas = VB_lambdas, VB_theta = VB_theta, nu = nu ) 

parameters <- c("beta", "beta_lambda1", "beta_lambda2", "beta_theta")
rstan_options(auto_write = TRUE)
samples <- stan(
  file = "Profiling.stan",
  data = data,
  pars = parameters,
  iter = 600,
  chains = 6,
  warmup = 300,
  thin = 1,
  seed = 123,
  refresh = 1,
  verbose = TRUE,
  control = list(
  adapt_delta = 0.999,
  max_treedepth = 10)
)

save(samples, file = "stan_chains.Rdata")
