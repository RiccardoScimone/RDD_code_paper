rm(list = ls())
graphics.off()
library(BayesNSGP)
library(ggplot2)
library(ggnewscale)
library(viridis)
source("Scripts/Utilities.r")
a = -1
b = 1
n_knot_grid = 6

knot_lon = seq(a, b, length.out = n_knot_grid+1)
knot_lat = seq(a, b, length.out = n_knot_grid+1)

knot_lon = (knot_lon[1:n_knot_grid] + knot_lon[2:(n_knot_grid+1)])/2
knot_lat = (knot_lat[1:n_knot_grid] + knot_lat[2:(n_knot_grid+1)])/2
knot_coord = expand.grid(knot_lon,knot_lat)

names(knot_coord) = c("x_1", "x_2")
simulation_id = "simulation_7"
dir_path = paste0("Simulations/",simulation_id)
load(paste0(dir_path,"/simulated_process_list.Rdata"))


clusters = 5
cl = makePSOCKcluster(clusters)
registerDoParallel(cl)
writeLines(c(""), paste0(dir_path,"/log.txt"))


i = 1
constants = list (Sigma_HP1 = c(10,10),Sigma_HP4 = c(10,20), nu = 2,
                  Sigma_knot_coords = knot_coord,
                  tau_knot_coords = knot_coord,
                  sigma_knot_coords = knot_coord,
                  Sigma_HP2 = c(5,5),Sigma_HP3 = c(3.85,3.85),mu_HP1 = 10, tau_HP1 = 10,maxAnisoRange = 16, k = 10)


sample = sample(1:10000,3000)

Rmodel = nsgpModel(likelihood = "NNGP", tau_model =  "approxGP", sigma_model = "approxGP",
                   Sigma_model = "npApproxGP", mu_model = "constant", constants = constants, 
                   coords = as.matrix(simulated_process_list[[1]][sample,c("x_1","x_2")]), data = simulated_process_list[[i]]$process[sample], nu = 2)

conf = configureMCMC(Rmodel)


conf$printSamplers()
conf$removeSamplers(ind = c(14:18))


cut_lat = cut(knot_coord[,1],breaks = 2)
cut_lon = cut(knot_coord[,2],breaks = 2)
cut = as.factor(paste0(cut_lat,cut_lon))
table(cut)
x11()
plot(knot_coord[,1],knot_coord[,2], col = cut)

groups = list()
for ( i in 1:4) groups[[i]] = which(as.numeric(cut) == i)
for(g in 1:length(groups)){
  conf$addSampler(target = c(paste0("w1_Sigma[", groups[[g]], "]"),
                             paste0("w2_Sigma[", groups[[g]], "]"),
                             paste0("w3_Sigma[", groups[[g]], "]"),
                             paste0("w_sigma[", groups[[g]], "]"),
                             paste0("w_tau[", groups[[g]], "]"))
                  , type = "RW_block", silent = TRUE)
}

conf$printSamplers()

Rmcmc = buildMCMC(conf)
Cmodel = compileNimble(Rmodel) # Compile the model in C++
Cmcmc = compileNimble(Rmcmc, project = Rmodel)

samples = runMCMC(Cmcmc, niter = 10000, nburnin = 5000, thin = 2,nchains = 1,progressBar = T)



saveRDS(file = paste0(dir_path,"/NSGP_results.rds"),object = samples)
saveRDS(file = paste0(dir_path,"/NSGP_knots.rds"),object = samples)

number_of_knots = nrow(knot_coord)
gridplot = as.matrix(simulated_process_list[[1]][,c("x_1","x_2")])

knot_dist_pred = sqrt(nsCrossdist(knot_coord, gridplot, isotropic = TRUE)$dist1_sq)

nu_latent = 5
cl = makePSOCKcluster(12)
registerDoParallel(cl)
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



lambda_1_samples = ret[,seq(1,748, by = 3)]
lambda_2_samples = ret[,seq(2,749, by = 3)]
theta_samples = ret[,seq(3,750, by = 3)]


save(list = c("lambda_1_samples",
              "lambda_2_samples",
              "theta_samples"), file = paste0(dir_path,"/Aniso_samples.rdata"))


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
estimates$theta = NULL

x11()
multiple_heatmaps(estimates)