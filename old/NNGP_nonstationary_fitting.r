rm(list = ls())
library(rstan)
library(shinystan)
library(spNNGP)  
library(tidyverse)# Build neighbor index
library(viridis)
source("NNmatrix.R") 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
### Compile models
sim_model = stan_model("simulazioni.stan")
model = stan_model("Non_stationary_fitting_nngp.stan")

### Simulate data
sq_N = 100
N = sq_N^2
a = -1
b = 1
grid = expand.grid(seq(a,b,length.out = sq_N),seq(a,b,length.out = sq_N))
names(grid) = c("x","y")
M = 5

## compute neighbors
NN.matrix <- NNMatrix(coords = grid, n.neighbors = M, n.omp.threads = 12)

## Extract indexes and reorder data
grid = grid[NN.matrix$ord,]
NN_ind = NN.matrix$NN_ind

## Simulate data
lambda_1 = rep(0.006,N);
lambda_2 = rep(0.003,N);
theta = rep(pi/3,N);
sigma = rep(1.5,N);
#mu = rep(1,N)* (1 - sqrt(grid$x^2 + grid$y^2)/sqrt(2)) + rep(2,N) * sqrt(grid$x^2 + grid$y^2)/sqrt(2);
nu = 2

simu_data = list(N = N, Locations = as.matrix(t(grid)), lambda_1 = lambda_1, lambda_2 = lambda_2, theta = theta, nu = nu, sigma = sigma)

simu_fit <- sampling(object = sim_model, data=simu_data,
                     warmup=0, iter=100, chains=1, seed=18011996,
                     algorithm="Fixed_param", refresh=10)


realizations = matrix( data = unlist(simu_fit@sim$samples), ncol = 100,byrow = TRUE)
grid_toplot = data.frame(x = grid$x, y = grid$y, process = realizations[1:N,5])
library(viridis)
x11()
p = ggplot(grid_toplot) + geom_tile(mapping = aes(x,y, fill = process)) + scale_fill_viridis(option = "inferno")
plot(p)

### Prepare data for fitting

P_mu = 0
P_sigma = 1
P_lambda = 1
P_theta = 1
fit_mu = 0
fit_sigma = 1
fit_lambda = 1
fit_theta = 1
N = sq_N^2
nu = 2
Y = grid_toplot$process
X_mu = array(0,c(0,N,P_mu))
X_sigma = array(matrix(1, ncol = P_sigma, nrow = N), c(1,N,P_sigma))
X_lambda = X_sigma
X_theta = X_sigma
Locations = t(grid)
uB_mu = array(0,c(0,P_mu))
uB_sigma = array(0,c(1,P_sigma))
uB_lambda = array(0,c(1,P_lambda))
uB_theta = array(0,c(1,P_theta))
VB_mu = array(diag(1,P_mu,P_mu),c(0,P_mu,P_mu))
VB_sigma = array(diag(1,P_sigma,P_sigma),c(1,P_sigma,P_sigma) )
VB_lambda = array(diag(1,P_lambda,P_lambda), c(1,P_lambda,P_lambda))
VB_theta = array(diag(1,P_theta,P_theta), c(1, P_theta, P_theta))


data = list(
  fit_mu = fit_mu,
  fit_sigma = fit_sigma,
  fit_lambda = fit_lambda,
  fit_theta = fit_theta,
  N = N,
  P_mu = P_mu,
  P_sigma = P_sigma,
  P_lambda = P_lambda,
  P_theta = P_theta,
  nu = nu,
  Y = Y,
  X_mu = X_mu,
  X_sigma = X_sigma,
  X_lambda = X_lambda,
  X_theta = X_theta,
  Locations = Locations,
  uB_mu = uB_mu,
  uB_sigma = uB_sigma,
  uB_lambda = uB_lambda,
  uB_theta = uB_theta,
  VB_mu = VB_mu,
  VB_sigma = VB_sigma,
  VB_lambda = VB_lambda,
  VB_theta = VB_theta,
  M = M,
  NN_ind = NN_ind
)

fit = sampling(object = model, data=data,
               warmup=500, iter=1500, chains=6, seed=18011996, refresh=10,sample_file = "chain_data_nngp.csv")
saveRDS(fit, file = "fit1_nngp.rds")

#### Other Model

theta = pi/2* rstanarm::invlogit(grid$x^2 + grid$y^2)
range(theta)
lambda_1 = exp(-7.5  + sqrt(grid$x^2 + grid$y^2))
lambda_2 = exp(-5.7 + 2*sqrt(grid$x^2 + grid$y^2))
sigma = exp(grid$x^2)
mu = rep(1,N);

simu_data = list(N = N, Locations = as.matrix(t(grid)), lambda_1 = lambda_1, lambda_2 = lambda_2, theta = theta, nu = nu, sigma = sigma)

simu_fit <- sampling(object = sim_model, data=simu_data,
                     warmup=0, iter=100, chains=1, seed=18011996,
                     algorithm="Fixed_param", refresh=10)


realizations = matrix( data = unlist(simu_fit@sim$samples), ncol = 100,byrow = TRUE)
for ( i in 1:100) for(j in 1:N) realizations[j,i] = realizations[j,i] + mu[j];

grid_toplot = data.frame(x = grid$x, y = grid$y, process = realizations[1:N,5])
x11()
p = ggplot(grid_toplot) + geom_tile(mapping = aes(x,y, fill = process)) + scale_fill_viridis(option = "inferno")
plot(p)

P_mu = 2
P_sigma = 3
P_lambda = 3
P_theta = 3
fit_mu = 1
fit_sigma = 1
fit_lambda = 1
fit_theta = 1
N = sq_N^2
nu = 2
Y = grid_toplot$process
X_mu = array(cbind(rep(1,N),rnorm(N)),c(1,N,P_mu))
X_sigma = array(cbind(rep(1,N),grid$x^2,rnorm(N)), c(1,N,P_sigma))
X_lambda = array(cbind(rep(1,N),sqrt(grid$x^2 + grid$y^2),rnorm(N)), c(1,N,P_lambda))
X_theta = array(cbind(grid$x^2,grid$y^2,rnorm(N)), c(1,N,P_theta))
Locations = t(grid)
uB_mu = array(0,c(1,P_mu))
uB_sigma = array(0,c(1,P_sigma))
uB_lambda = array(c(-3,0),c(1,P_lambda))
uB_theta = array(0,c(1,P_theta))
VB_mu = array(diag(1,P_mu,P_mu),c(1,P_mu,P_mu))
VB_sigma = array(diag(1,P_sigma,P_sigma),c(1,P_sigma,P_sigma) )
VB_lambda = array(diag(1,P_lambda,P_lambda), c(1,P_lambda,P_lambda))
VB_theta = array(diag(1,P_theta,P_theta), c(1, P_theta, P_theta))


data = list(
  fit_mu = fit_mu,
  fit_sigma = fit_sigma,
  fit_lambda = fit_lambda,
  fit_theta = fit_theta,
  N = N,
  P_mu = P_mu,
  P_sigma = P_sigma,
  P_lambda = P_lambda,
  P_theta = P_theta,
  nu = nu,
  Y = Y,
  X_mu = X_mu,
  X_sigma = X_sigma,
  X_lambda = X_lambda,
  X_theta = X_theta,
  Locations = Locations,
  uB_mu = uB_mu,
  uB_sigma = uB_sigma,
  uB_lambda = uB_lambda,
  uB_theta = uB_theta,
  VB_mu = VB_mu,
  VB_sigma = VB_sigma,
  VB_lambda = VB_lambda,
  VB_theta = VB_theta,
  M = M,
  NN_ind = NN_ind
)

fit = sampling(object = model, data=data,
               warmup=200, iter=800, chains=6, seed=18011996, refresh=10,sample_file = "chain_data_nngp_complex.csv",control = list(adapt_delta = 0.8))
saveRDS(fit, file = "fit2_nngp_more.rds")


### no nuisances
P_mu = 1
P_sigma = 2
P_lambda = 2
P_theta = 2
fit_mu = 1
fit_sigma = 1
fit_lambda = 1
fit_theta = 1
N = sq_N^2
nu = 2
Y = grid_toplot$process
X_mu = array(cbind(rep(1,N)),c(1,N,P_mu))
X_sigma = array(cbind(rep(1,N),grid$x^2), c(1,N,P_sigma))
X_lambda = array(cbind(rep(1,N),sqrt(grid$x^2 + grid$y^2)), c(1,N,P_lambda))
X_theta = array(cbind(grid$x^2,grid$y^2), c(1,N,P_theta))
Locations = t(grid)
uB_mu = array(0,c(1,P_mu))
uB_sigma = array(0,c(1,P_sigma))
uB_lambda = array(c(-3,0),c(1,P_lambda))
uB_theta = array(0,c(1,P_theta))
VB_mu = array(diag(1,P_mu,P_mu),c(1,P_mu,P_mu))
VB_sigma = array(diag(1,P_sigma,P_sigma),c(1,P_sigma,P_sigma) )
VB_lambda = array(diag(1,P_lambda,P_lambda), c(1,P_lambda,P_lambda))
VB_theta = array(diag(1,P_theta,P_theta), c(1, P_theta, P_theta))


data = list(
  fit_mu = fit_mu,
  fit_sigma = fit_sigma,
  fit_lambda = fit_lambda,
  fit_theta = fit_theta,
  N = N,
  P_mu = P_mu,
  P_sigma = P_sigma,
  P_lambda = P_lambda,
  P_theta = P_theta,
  nu = nu,
  Y = Y,
  X_mu = X_mu,
  X_sigma = X_sigma,
  X_lambda = X_lambda,
  X_theta = X_theta,
  Locations = Locations,
  uB_mu = uB_mu,
  uB_sigma = uB_sigma,
  uB_lambda = uB_lambda,
  uB_theta = uB_theta,
  VB_mu = VB_mu,
  VB_sigma = VB_sigma,
  VB_lambda = VB_lambda,
  VB_theta = VB_theta,
  M = M,
  NN_ind = NN_ind
)

fit = sampling(object = model, data=data,
               warmup=200, iter=800, chains=6, seed=18011996, refresh=10,sample_file = "chain_data_nngp_complex_no_nuis.csv",control = list(adapt_delta = 0.8))
saveRDS(fit, file = "fit2_nngp.rds")



#### Other Model

theta = pi/2 * rstanarm::invlogit(grid$x^2 + grid$y^2)
range(theta)
lambda_1 = exp ( (-7 - grid$x)) 
lambda_2 = exp ( (-8 + grid$x))
sigma = exp(1/4*grid$x)
mu = grid$x + grid$y
theta_deg = theta/pi*180
nu = 2
simu_data = list(N = N, Locations = as.matrix(t(grid)), lambda_1 = lambda_1, lambda_2 = lambda_2, theta = theta, nu = nu, sigma = sigma)

simu_fit <- sampling(object = sim_model, data=simu_data,
                     warmup=0, iter=100, chains=1, seed=18011996,
                     algorithm="Fixed_param", refresh=10)

realizations = matrix( data = unlist(simu_fit@sim$samples), ncol = 100,byrow = TRUE)
for ( i in 1:100) for(j in 1:N) realizations[j,i] = realizations[j,i] + mu[j];

grid_toplot = data.frame(x = grid$x, y = grid$y, process = realizations[1:N,20])
x11()
p = ggplot(grid_toplot) + geom_tile(mapping = aes(x,y, fill = process)) + scale_fill_viridis(option = "inferno") + coord_fixed()
plot(p)

X = array(cbind(rep(1,N),grid$x,grid$y,grid$x^2,grid$y^2, grid$x*grid$y),c(1,N,6))
P_mu = 6
P_sigma = 6
P_lambda = 6
P_theta = 6
fit_mu = 1
fit_sigma = 1
fit_lambda = 1
fit_theta = 1
N = sq_N^2
nu = 2
Y = grid_toplot$process
X_mu = X
X_sigma = X
X_lambda = X
X_theta = X
Locations = t(grid)
uB_mu = array(0,c(1,P_mu))
uB_sigma = array(0,c(1,P_sigma))
uB_lambda = array(c(-3,0),c(1,P_lambda))
uB_theta = array(0,c(1,P_theta))
VB_mu = array(diag(1,P_mu,P_mu),c(1,P_mu,P_mu))
VB_sigma = array(diag(1,P_sigma,P_sigma),c(1,P_sigma,P_sigma) )
VB_lambda = array(diag(1,P_lambda,P_lambda), c(1,P_lambda,P_lambda))
VB_theta = array(diag(1,P_theta,P_theta), c(1, P_theta, P_theta))


data = list(
  fit_mu = fit_mu,
  fit_sigma = fit_sigma,
  fit_lambda = fit_lambda,
  fit_theta = fit_theta,
  N = N,
  P_mu = P_mu,
  P_sigma = P_sigma,
  P_lambda = P_lambda,
  P_theta = P_theta,
  nu = nu,
  Y = Y,
  X_mu = X_mu,
  X_sigma = X_sigma,
  X_lambda = X_lambda,
  X_theta = X_theta,
  Locations = Locations,
  uB_mu = uB_mu,
  uB_sigma = uB_sigma,
  uB_lambda = uB_lambda,
  uB_theta = uB_theta,
  VB_mu = VB_mu,
  VB_sigma = VB_sigma,
  VB_lambda = VB_lambda,
  VB_theta = VB_theta,
  M = M,
  NN_ind = NN_ind
)

fit = sampling(object = model, data=data,
               warmup=500, iter=1000, chains=6, seed=18011996, refresh=10,sample_file = "chain_data_easy.csv",control = list(adapt_delta = 0.9))
saveRDS(fit, file = "fit3_nngp.rds")

