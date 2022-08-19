### Non stationary simulations
rm(list = ls())
library(tidyverse)
library(rstan)
library(rstanarm)
library(tictoc)
library(cmdstanr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
graphics.off()
sq_N = 30
N = sq_N^2
a = -0.3
b = 0.3
grid = expand.grid(seq(a,b,length.out = sq_N),seq(a,b,length.out = sq_N))
names(grid) = c("x","y")
theta = (atan(grid$y/grid$x))*as.numeric(grid$x*grid$y >= 0) + (pi/2 + atan(grid$y/grid$x)) * as.numeric(grid$x*grid$y < 0)
range(theta)
lambda_1 = 0.0005 * as.numeric(grid$x*grid$y >= 0) + 0.003 * as.numeric(grid$x*grid$y < 0)
lambda_2 = 0.0005 * as.numeric(grid$x*grid$y < 0) + 0.003 * as.numeric(grid$x*grid$y >= 0)
sigma = rep(1,N)* (1 - sqrt(grid$x^2 + grid$y^2)/sqrt(2)) + rep(2,N) * sqrt(grid$x^2 + grid$y^2)/sqrt(2)
#lambda_1 = rep(0.0005,N);
#lambda_2 = rep(0.003,N);
#theta = rep(pi/3,N);
sigma = rep(1.5,N);
mu = rep(1,N)* (1 - sqrt(grid$x^2 + grid$y^2)/sqrt(2)) + rep(2,N) * sqrt(grid$x^2 + grid$y^2)/sqrt(2);
nu = 3
logit(2/pi*pi/3)
model = stan_model("simulazioni.stan",verbose = T)
cm_model = cmdstan_model("simulazioni.stan",cpp_options = list(stan_opencl = TRUE),include_paths = c("."))

#### If you want to debug/measure performance
#expose_stan_functions(model)
#tic()
#Aniso = Compute_Aniso(lambda_1,lambda_2,theta)
#Dets = Compute_Determinants(lambda_1,lambda_2)
#Cov = matern_ns_corr_mat(t(grid),Aniso,Dets,nu,sigma)
#toc()


data = list(N = N, Locations = as.matrix(t(grid)), lambda_1 = lambda_1, lambda_2 = lambda_2, theta = theta, nu = nu, sigma = sigma)

fit_cl <- cm_model$sample(data = data, chains = 1, parallel_chains = 1, fixed_param = T,iter_warmup = 0,iter_sampling = 100,opencl_ids = c(0, 0), refresh = 10)


simu_fit <- sampling(object = model, data=data,
                 warmup=0, iter=100, chains=1, seed=494838,
                 algorithm="Fixed_param", refresh=10)

realizations = matrix( data = unlist(simu_fit@sim$samples), ncol = 100,byrow = TRUE)

for ( i in 1:100) for(j in 1:N) realizations[j,i] = realizations[j,i] + mu[j];

grid_toplot = data.frame(x = grid$x, y = grid$y, process = realizations[1:N,5])

library(viridis)
x11()
p = ggplot(grid_toplot) + geom_tile(mapping = aes(x,y, fill = process)) + scale_fill_viridis(option = "inferno")
plot(p)




#### Now let's try with fitting. We start with zero mean stationary anisotropic.
model_fitting = stan_model("Non_stationary_fitting.stan")
cm_model_fitting = cmdstan_model("Non_stationary_fitting.stan",cpp_options = list(stan_opencl = TRUE),include_paths = c("."))

P_mu = 2
P_sigma = 1
P_lambda = 2
P_theta = 1
fit_mu = 1
fit_sigma = 1
fit_lambda = 1
fit_theta = 1
N = sq_N^2
nu = 2
Y = grid_toplot$process
X_mu = array(cbind(rep(1,N), sqrt(grid$x^2 + grid$y^2)),c(1,N,P_mu))
X_sigma = array(matrix(1, ncol = P_sigma, nrow = N), c(1,N,P_sigma))
X_lambda = array(cbind(as.numeric(grid$x*grid$y>=0),as.numeric(grid$x*grid$y<0)), c(1,N,P_lambda))
X_theta = array(cbind(invlogit(2/pi*theta)), c(1,N,P_theta))
Locations = t(grid)
uB_mu = array(0,c(1,P_mu))
uB_sigma = array(0,c(1,P_sigma))
uB_lambda = array(0,c(1,P_lambda))
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
VB_theta = VB_theta
)

#cm_fit_fit = cm_model_fitting$sample(data = data, chains = 3, parallel_chains = 3, fixed_param = F,iter_warmup = 1000,iter_sampling = 1000,opencl_ids = c(0, 0), refresh = 10)


fit_fit2 <- sampling(object = model_fitting, data=data,
                     warmup=1000, iter=2000, chains=4, seed=18011996, refresh=10)
saveRDS(fit_fit2,"fit_fit2.rds")

fit_fit2 = readRDS("fit_fit2.rds")
launch_shinystan(fit_fit2)

