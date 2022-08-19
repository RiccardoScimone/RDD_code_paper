rm(list = ls())
library(rstan)
library(plotly)
source("Scripts/NNMatrix.R")
rstan::rstan_options(auto_write = T)

options(mc.cores = parallel::detectCores(),Ncpus =  parallel::detectCores())

simulation_id = "simulation_1"
dir_path = paste0("Simulations/",simulation_id)
K = 16
est_spatial_pars = readRDS(paste0(dir_path,"/RDD_estimate",K,".rds"))
load(paste0(dir_path,"/simulated_process.Rdata"))

to_fit = "lambda_2"
rescale_log = T
M = 5
NN.matrix <- NNMatrix(coords = simulated_process[,1:2], n.neighbors = M, n.omp.threads = 12)

## Extract indexes and reorder data
simulated_process = simulated_process[NN.matrix$ord,]
NN_ind = NN.matrix$NN_ind

N_obs = 10000
x_obs = as.matrix(simulated_process[1:2])
y_obs = est_spatial_pars[[to_fit]][NN.matrix$ord]
if(rescale_log) y_obs = log(y_obs)

p = plot_ly(x = x_obs[,1], y = x_obs[,2], z = y_obs, type = "scatter3d", size = 1)

stan_data = list(N_obs = N_obs, x_obs = x_obs, y_obs = y_obs,NN_ind = NN_ind, M = M)
model = readRDS("Stan_models/Fit_mvariate_stationary_gp.rds")

fit = rstan::sampling(object = model, data=stan_data,
                      warmup=300, iter=700, chains=3, seed=18011996, 
                      refresh=10,sample_file = paste0(dir_path,"/chain_data_fitRDD",to_fit,".csv"),
                      control = list(adapt_delta = 0.8))

saveRDS(fit, file = paste0("fit_",to_fit,"_rdd_gp.rds"))
