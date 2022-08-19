rm(list = ls())
source("Scripts/Utilities.r")
source("Scripts/NNmatrix.R") 

##Specify simulation id

sim_id = "simulation_3"
dir = paste0("Simulations/",sim_id,"/")
load(paste0(dir,"simulated_process.Rdata"))

#Find neighbors
M = 5
NN.matrix <- NNMatrix(coords = simulated_process[,1:2], n.neighbors = M, n.omp.threads = 12)

## Extract indexes and reorder data
simulated_process = simulated_process[NN.matrix$ord,]
NN_ind = NN.matrix$NN_ind
saveRDS(NN.matrix$ord,paste0(dir,"ord.rds"))


## Create covariates 
X = model.matrix(data = simulated_process, object = process ~ 1 + x_1 + x_2 + I(x_1*x_2) + I(x_1^2) + I(x_2^2))

## Specify nu
nu = 2

## Create list for stan 
data = create_stan_list(covariate = X, spatial_data = simulated_process, nu = nu, NN_ind = NN_ind, M = M) 


## Fit the model
fitting_model = readRDS(file = "Stan_models/Non_stationary_fitting_nngp.rds")


fit = rstan::sampling(object = fitting_model, data=data,
               warmup=200, iter=600, chains=5, seed=18011996, 
               refresh=10,sample_file = paste0(dir,"chain_data_nngp.csv"),
               control = list(adapt_delta = 0.8))


saveRDS(fit, file = paste0(dir,"fitted.rds"))
saveRDS(data, file = paste0(dir,"fitting_input.rds"))

