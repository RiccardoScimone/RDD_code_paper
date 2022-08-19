rm(list = ls())
library(rstan)
options(mc.cores = parallel::detectCores(),Ncpus = parallel::detectCores())
rstan_options(auto_write = TRUE)
setwd("Stan_models")
Stan_simulate_ns_matern = stan_model("simulate.stan",includes = paste0('\n#include "', 
                                                                      file.path(getwd(), 'newbesselK.hpp'), '"\n'),allow_undefined = TRUE)
Stan_fit_linear_nngp_ns_matern = stan_model("Non_stationary_fitting_nngp.stan")
Stan_fit_linear_full_ns_matern = stan_model("Non_stationary_fitting.stan")
Stan_fit_linear_nngo_stationary_matern = stan_model("Fit_mvariate_stationary_gp.stan")


Stan_ns_Kriging = stan_model("nsKriging.stan",
                                  includes = paste0('\n#include "', file.path(getwd(), 'newbesselK.hpp'), '"\n'),allow_undefined = TRUE)

setwd("../")
