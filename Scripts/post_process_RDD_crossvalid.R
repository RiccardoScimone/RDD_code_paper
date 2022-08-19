### Analysis on the posterior distributions of spatial parameters
rm(list = ls())
graphics.off()
source("Scripts/Utilities.r")
expose_stan_functions("Stan_models/to_expose.stan")
library(SpatialTools)

dir_path = "GHCN"
K = 16

load(paste0(dir_path, "/RDD_fitting_Crossvalid",K,".Rdata"))

### Make predictions on folds

crossvalid_predictions = list()
mses = NULL
mse_nocv = NULL
obs = read.csv("GHCN/average_dataset_risser.csv") %>% select(LON, LAT, logprep)
names(obs)[1:2] = c("x_1","x_2")
for( i in 1:length(unique(folds))){
  idxs_fold = which(folds == i)
  idxs_nofold = which(folds != i)
  Locations = rbind ( obs$x_1, obs$x_2)
  lambda_1 =  RDD_estimate[[i]]$lambda_1
  lambda_2 =  RDD_estimate[[i]]$lambda_2
  theta =  RDD_estimate[[i]]$theta
  Anis = Compute_Aniso(lambda_1,lambda_2, theta)
  sigma = RDD_estimate[[i]]$sigma
  C = matern_ns_corr_mat(Locations, Anis, lambda_1*lambda_2, 5/2 , sigma)
  for ( j in 1:2311) C[j,j] = C[j,j] + 0.05
  krig = krige.sk(obs$logprep[idxs_nofold] - RDD_estimate[[i]]$mu[idxs_nofold], C[idxs_nofold,idxs_nofold],C[idxs_fold, idxs_fold], C[idxs_nofold,idxs_fold],m = 0)
  
  ret = list(krigpred = krig$pred + RDD_estimate[[i]]$mu[idxs_fold] )
  ret$Msecv = mean ( (ret$krigpred - obs$logprep[idxs_fold])^2)
  crossvalid_predictions[[i]] = ret
  mses = c(mses, ret$Msecv)
  mse_nocv = c(mse_nocv, ret$mse)
  print(i)
  }

