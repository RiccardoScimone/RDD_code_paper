rm(list = ls())
graphics.off()
library(BayesNSGP)
library(ggplot2)
library(ggnewscale)
library(viridis)

Colorado_data = read.csv("Colorado/Colorado_data.csv")
Colorado_grid = read.csv("Colorado/Colorado_grid.csv")

n_knot_grid = 8
knot_lon = seq(min(Colorado_data$Longitude), max(Colorado_data$Longitude), length.out = n_knot_grid+1)
knot_lat = seq(min(Colorado_data$Latitude), max(Colorado_data$Latitude), length.out = n_knot_grid+1)

knot_lon = (knot_lon[1:n_knot_grid] + knot_lon[2:(n_knot_grid+1)])/2
knot_lat = (knot_lat[1:n_knot_grid] + knot_lat[2:(n_knot_grid+1)])/2
knot_coord = expand.grid(knot_lon,knot_lat)
names(knot_coord) = c("Longitude", "Latitude")
x11()
ggplot() + geom_point(data = Colorado_data,
                      aes(x = Longitude, y = Latitude, col = logPrecip)) +scale_color_viridis(direction = -1)+ 
          geom_point(data = knot_coord, aes(x = Longitude, y = Latitude),size = 2, shape = 4)

knot_coord = as.matrix(knot_coord)

N = nrow(Colorado_data)

constants = list (Sigma_HP1 = c(10,10),Sigma_HP4 = c(10,20), nu = 2,
                  Sigma_knot_coords = knot_coord,Sigma_HP2 = c(5,5),Sigma_HP3 = c(3.85,3.85),mu_HP1 = 10, tau_HP1 = 10,maxAnisoRange = 16, k = 10)



Rmodel = nsgpModel(likelihood = "NNGP", tau_model =  "constant", sigma_model = "constant",
                      Sigma_model = "npApproxGP", mu_model = "constant", constants = constants, 
                   coords = as.matrix(Colorado_data[,c("Longitude","Latitude")]), data = Colorado_data$logPrecip, nu = 2)
conf = configureMCMC(Rmodel)

conf$printSamplers()
conf$removeSamplers(ind = c(10,11,12))


# Add block samplers
cut_lat = cut(knot_coord[,1],breaks = 4)
cut_lon = cut(knot_coord[,2],breaks = 2)
cut = as.factor(paste0(cut_lat,cut_lon))
table(cut)
x11()
plot(knot_coord[,1],knot_coord[,2], col = cut)

groups = list()
for ( i in 1:8) groups[[i]] = which(as.numeric(cut) == i)
for(g in 1:length(groups)){
  conf$addSampler(target = c(paste0("w1_Sigma[", groups[[g]], "]"),paste0("w2_Sigma[", groups[[g]], "]"),paste0("w3_Sigma[", groups[[g]], "]"))
                  , type = "RW_block", silent = TRUE)
}

conf$printSamplers()

Rmcmc = buildMCMC(conf)
Cmodel = compileNimble(Rmodel) # Compile the model in C++
Cmcmc = compileNimble(Rmcmc, project = Rmodel)

samples = runMCMC(Cmcmc, niter = 100000, nburnin = 50000, thin = 10,nchains = 1,progressBar = T)
save(samples,file = "Colorado/samples2.Rdata")
save(knot_coord, file = "Colorado/knots.Rdata")

rm(list = ls())
load("Colorado/samples2.Rdata")
load("Colorado/knots.Rdata")
Colorado_grid = read.csv("Colorado/Colorado_grid.csv")

### Computation for anisotropy
### retrieve eigenprocess on a grid

nu_latent = 5
number_of_knots = nrow(knot_coord)
gridplot = as.matrix(Colorado_grid[,c("longitude","latitude")])
rm(Colorado_grid)
knot_dist_pred = sqrt(nsCrossdist(knot_coord, gridplot, isotropic = TRUE)$dist1_sq)

library(foreach)
library(doParallel)
cl = makePSOCKcluster(12)
registerDoParallel(cl)
dir_path = "Colorado"
writeLines(c(""), paste0(dir_path,"/log.txt"))
ret = foreach(j = 1:nrow(samples),.packages = c("BayesNSGP"), .combine = "cbind")%dopar%{
  sink(paste0(dir_path,"/log.txt"),append = T)
  print(j)
  sink()
  w_lambda1_j = samples[j,paste("w1_Sigma[",1:number_of_knots,"]",sep = "")]
  w_lambda2_j = samples[j,paste("w2_Sigma[",1:number_of_knots,"]",sep = "")]
  w_theta_j =  samples[j,paste("w3_Sigma[",1:number_of_knots,"]",sep = "")]
  Pmat_theta_j = matern_corr(knot_dist_pred, samples[j,"SigmaGP_phi[2]"], nu_latent)
  Pmat_eig_j = matern_corr(knot_dist_pred, samples[j,"SigmaGP_phi[1]"], nu_latent)
  cbind(samples[j,"SigmaGP_mu[1]"] + samples[j,"SigmaGP_sigma[1]"] * Pmat_eig_j %*% w_lambda1_j,
  samples[j,"SigmaGP_mu[1]"] + samples[j,"SigmaGP_sigma[1]"] * Pmat_eig_j %*% w_lambda2_j,
  samples[j,"SigmaGP_mu[2]"] + samples[j,"SigmaGP_sigma[2]"] * Pmat_theta_j %*% w_theta_j)
}
stopCluster(cl)

lambda_1_samples = ret[,seq(1,14998, by = 3)]
lambda_2_samples = ret[,seq(2,14999, by = 3)]
theta_samples = ret[,seq(3,15000, by = 3)]


save(list = c("lambda_1_samples","lambda_2_samples","theta_samples"), file = "Colorado/Aniso_samples2.rdata")
load("Colorado/Aniso_samples.rdata")
source("Scripts/Utilities.r")

lambda_1_samples = exp(lambda_1_samples)
lambda_2_samples = exp(lambda_2_samples)
theta_samples = pi/2 * rstanarm::invlogit(theta_samples)

## Estimate
lambda_1_estimates = rowMeans(lambda_1_samples)
lambda_2_estimates = rowMeans(lambda_2_samples)
theta_estimates = rowMeans(theta_samples)


estimates = data.frame(x_1 = gridplot[,1], x_2 = gridplot[,2], 
                       lambda_1 = lambda_1_estimates, lambda_2 = lambda_2_estimates, theta = theta_estimates)


estimates$theta_deg = 90- 180/pi * estimates$theta
x11()
multiple_heatmaps(estimates)
x11()
plot_ellipses(data = estimates)


library(tidyverse)
## study just interesting samples

samples_f = data.frame(samples) %>% dplyr::select(!(starts_with("w")))

x11()
par(mfrow = c(2,5))
for ( i in 1:9)
  plot(samples[,i], type = "l", main = colnames(samples)[i], xlab = "Iteration", ylab = "")

