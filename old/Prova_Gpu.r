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
sq_N = 100
N = sq_N^2
a = -0.1
b = 0.1
grid = expand.grid(seq(a,b,length.out = sq_N),seq(a,b,length.out = sq_N))
names(grid) = c("x","y")
theta = (atan(grid$y/grid$x))*as.numeric(grid$x*grid$y >= 0) + (pi/2 + atan(grid$y/grid$x)) * as.numeric(grid$x*grid$y < 0)
range(theta)
#lambda_1 = 0.0005 * as.numeric(grid$x*grid$y >= 0) + 0.003 * as.numeric(grid$x*grid$y < 0)
#lambda_2 = 0.0005 * as.numeric(grid$x*grid$y < 0) + 0.003 * as.numeric(grid$x*grid$y >= 0)
#sigma = rep(1,N)* (1 - sqrt(grid$x^2 + grid$y^2)/sqrt(2)) + rep(2,N) * sqrt(grid$x^2 + grid$y^2)/sqrt(2)
lambda_1 = rep(0.0005,N);
lambda_2 = rep(0.003,N);
theta = rep(pi/3,N);
sigma = rep(1.5,N);
nu = 3
logit(2/pi*pi/3)

# we recommend running this is a fresh R session or restarting your current session
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
path_to_opencl_lib <- "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.5/lib/x64"
cpp_options = list(
  "CXXFLAGS += -fpermissive -O3 -Wno-unused-variable -Wno-unused-function",
  "PRECOMPILED_HEADERS"=FALSE,
  paste0("LDFLAGS+= -L\"",path_to_opencl_lib,"\" -lOpenCL")
)

install_cmdstan(cores=12, overwrite = TRUE, cpp_options = cpp_options)

#Sys.setenv(PATH = paste(old_path, "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v11.5\\lib\\x64\\;C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v11.5\\include\\;C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v11.5\\include\\CL\\", sep = ";"))

model_cl <- cmdstan_model("simulazioni.stan",
                        cpp_options = list(stan_opencl = TRUE),include_paths = c("."), quiet = F, pedantic = T)

data = list(N = N, Locations = as.matrix(t(grid)), lambda_1 = lambda_1, lambda_2 = lambda_2, theta = theta, nu = nu, sigma = sigma)

fit_cl <- model_cl$sample(data = data, chains = 1, parallel_chains = 1, fixed_param = T,iter_warmup = 0,iter_sampling = 500,opencl_ids = c(0, 0), refresh = 0)

#### If you want to debug/measure performance
#expose_stan_functions(model)
#tic()
#Aniso = Compute_Aniso(lambda_1,lambda_2,theta)
#Dets = Compute_Determinants(lambda_1,lambda_2)
#Cov = matern_ns_corr_mat(t(grid),Aniso,Dets,nu,sigma)
#toc()

realizations = matrix( data = unlist(simu_fit@sim$samples), ncol = 10,byrow = TRUE)

grid_toplot = data.frame(x = grid$x, y = grid$y, process = realizations[1:N,5])

library(viridis)
x11()
p = ggplot(grid_toplot) + geom_tile(mapping = aes(x,y, fill = process)) + scale_fill_viridis(option = "inferno")
plot(p)









#### Debugging opencl

n<- 250000
k <- 20
X <- matrix(rnorm(n * k), ncol = k)
y <- rbinom(n, size = 1, prob = plogis(3 * X[,1] - 2 * X[,2] + 1))
mdata <- list(k = k, n = n, y = y, X = X)
mod_cl <- cmdstan_model("prova_gpu.stan",
                        cpp_options = list(stan_opencl = TRUE))

fit_cl <- mod_cl$sample(data = mdata, chains = 4, parallel_chains = 4,
                        opencl_ids = c(0, 0), refresh = 0)


N = 1500
x <- 22 * (0:(N - 1)) / (N - 1) - 11

alpha_true <- 3
rho_true <- 5.5
sigma_true = 1e-1
simu_data <- list(alpha=alpha_true, rho=rho_true, N=N, x=x)

model_cl <- cmdstan_model("simu.stan",
                          cpp_options = list(stan_opencl = TRUE),include_paths = c("."), quiet = F, pedantic = T)


fit_cl <- model_cl$sample(data = simu_data, chains = 1, fixed_param = T,iter_warmup = 0,iter_sampling = 4000,opencl_ids = c(0, 0), refresh = 4000)



simu_fit <- stan(file='simu.stan', data=simu_data,
                 warmup=0, iter=4000, chains=1, seed=494838,
                 algorithm="Fixed_param", refresh=4000,verbose = T)

realizations = matrix( data = unlist(simu_fit@sim$samples), ncol = 10,byrow = TRUE)

idxs = seq(1,1500, by = 1)
x_obs = x[idxs]
y_obs = realizations[idxs,1] + rnorm(1500,1e-2)
simu_data <- list(N_obs = 1500,x_obs = x_obs, y_obs = y_obs)

simu_data$alpha <- alpha_true
simu_data$rho <- rho_true
simu_data$sigma <- sigma_true
model_cl <- cmdstan_model("infer_gp.stan",
                          cpp_options = list(stan_opencl = TRUE),include_paths = c("."), quiet = F, pedantic = T)


fit_cl <- model_cl$sample(data = simu_data, chains = 4, parallel_chains = 4, fixed_param = F,iter_warmup = 1000,iter_sampling = 1000,opencl_ids = c(0, 0))

