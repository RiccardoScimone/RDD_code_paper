rm(list = ls())
graphics.off()
source("Scripts/Utilities.r")
library(BayesNSGP)
simulation_id = "simulation_1"
dir_path = paste0("Simulations/",simulation_id)
load(paste0(dir_path,"/simulated_process.Rdata"))
sim_const = simulated_process
for ( N in c(1000,2000,5000,8000,10000)){
if(N < nrow(simulated_process)){
set.seed(18011996)
simulated_process = sim_const[sample(1:nrow(simulated_process),N),]
}
for ( K in c(1,4,8,16,32)){
B = 120
if (K == 1) B = 1
struct = c("K-Bessel", "Nugget Effect")
param = c("P=2")
lower = c("P=2")
upper = c("P=2")
dirs = seq(0,150, by = 30)
estimates = RDD_fit(data = simulated_process$process,coords = simulated_process[,1:2],center_grid = simulated_process[,1:2],
K = K, B = B, clusters = 12,dir_path = dir_path , struct = struct,dirs = dirs,
param = param,lower = lower, upper = upper, threshold = 10)
load(paste0(dir_path,"/simulated_process.Rdata"))
ret = RDD_predict(estimates,grid = sim_const[,1:2])
ret = simplify2array(ret)
RDD_estimate = extract_spatial_predictions_median(ret)
saveRDS(ret,paste0(dir_path,"/RDD_fitting",K,N,".rds"))
saveRDS(RDD_estimate,paste0(dir_path,"/RDD_estimate",K,N,".rds"))
print(paste0(N,"_",K))
}
}

### Bayes NSGP

a = -1
b = 1

n_knot_grid = 4
knot_lon = seq(a, b, length.out = n_knot_grid+1)
knot_lat = seq(a, b, length.out = n_knot_grid+1)

knot_lon = (knot_lon[1:n_knot_grid] + knot_lon[2:(n_knot_grid+1)])/2
knot_lat = (knot_lat[1:n_knot_grid] + knot_lat[2:(n_knot_grid+1)])/2
knot_coord = expand.grid(knot_lon,knot_lat)
names(knot_coord) = c("x_1", "x_2")

x11()
ggplot() + geom_tile(data = simulated_process,
                      aes(x = x_1, y = x_2, fill = process)) +scale_fill_viridis(option = "inferno")+ coord_fixed() +
  geom_point(data = knot_coord, aes(x = x_1, y = x_2),size = 2, shape = 4)

knot_coord = as.matrix(knot_coord)


constants = list (nu = 2, Sigma_knot_coords = knot_coord,k = 5,
                  tau_knot_coords = knot_coord,sigma_knot_coords = knot_coord)

library(BayesNSGP)
Rmodel = nsgpModel(likelihood = "NNGP", tau_model =  "approxGP", sigma_model = "approxGP",
                   Sigma_model = "npApproxGP", mu_model = "constant", constants = constants, 
                   coords = as.matrix(simulated_process[,1:2]), data = simulated_process$process, nu = 2)


conf = configureMCMC(Rmodel)

conf$printSamplers()

conf$removeSamplers(ind = 14:18)


cut_lat = cut(knot_coord[,1],breaks = 2)
cut_lon = cut(knot_coord[,2],breaks = 2)
cut = as.factor(paste0(cut_lat,cut_lon))
table(cut)
x11()
plot(knot_coord[,1],knot_coord[,2], col = cut)

groups = list()
for ( i in 1:4) groups[[i]] = which(as.numeric(cut) == i)
for(g in 1:length(groups)){
  conf$addSampler(target = c(paste0("w_tau[", groups[[g]], "]"),paste0("w_sigma[", groups[[g]], "]"),
    paste0("w1_Sigma[", groups[[g]], "]"),paste0("w2_Sigma[", groups[[g]], "]"),paste0("w3_Sigma[", groups[[g]], "]"))
                  , type = "RW_block", silent = TRUE)
}

conf$printSamplers()



Rmcmc = buildMCMC(conf)
Cmodel = compileNimble(Rmodel) # Compile the model in C++
Cmcmc = compileNimble(Rmcmc, project = Rmodel)

tictoc::tic()
samples = runMCMC(Cmcmc, niter = 100000, nburnin = 50000, thin = 10,nchains = 1,progressBar = T)
tictoc::toc()

save(samples,file = paste0(dir_path,"/samples.Rdata"))
save(knot_coord, file = paste0(dir_path,"/knots.Rdata"))

load(paste0(dir_path,"/samples.Rdata"))

nu_latent = 5
number_of_knots = nrow(knot_coord)
load(paste0(dir_path,"/simulated_process.Rdata"))
gridplot = as.matrix(simulated_process[,1:2])
knot_dist_pred = sqrt(nsCrossdist(knot_coord, gridplot, isotropic = TRUE)$dist1_sq)

library(foreach)
library(doParallel)
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
  w_tau_j =  samples[j,paste("w_tau[",1:number_of_knots,"]",sep = "")]
  w_sigma_j = samples[j,paste("w_sigma[",1:number_of_knots,"]",sep = "")]
  Pmat_theta_j = matern_corr(knot_dist_pred, samples[j,"SigmaGP_phi[2]"], nu_latent)
  Pmat_eig_j = matern_corr(knot_dist_pred, samples[j,"SigmaGP_phi[2]"], nu_latent)
  Pmat_tau_j = matern_corr(knot_dist_pred, samples[j,"SigmaGP_phi[2]"], nu_latent)
  Pmat_sigma_j = matern_corr(knot_dist_pred, samples[j,"SigmaGP_phi[2]"], nu_latent)
  cbind(samples[j,"SigmaGP_mu[1]"] + samples[j,"SigmaGP_sigma[1]"] * Pmat_eig_j %*% w_lambda1_j,
        samples[j,"SigmaGP_mu[1]"] + samples[j,"SigmaGP_sigma[1]"] * Pmat_eig_j %*% w_lambda2_j,
        samples[j,"SigmaGP_mu[2]"] + samples[j,"SigmaGP_sigma[2]"] * Pmat_theta_j %*% w_theta_j,
        samples[j,"tauGP_mu"] + samples[j,"tauGP_sigma"] * Pmat_tau_j %*% w_tau_j,
        samples[j,"sigmaGP_mu"] + samples[j,"sigmaGP_sigma"] * Pmat_sigma_j %*% w_sigma_j)
}
stopCluster(cl)

lambda_1_samples = ret[,seq(1,24996, by = 5)]
lambda_2_samples = ret[,seq(2,24997, by = 5)]
theta_samples = ret[,seq(3,24998, by = 5)]
tau_samples = ret[,seq(4,24999, by = 5)]
sigma_samples = ret[,seq(5,25000, by = 5)]

save(list = c("lambda_1_samples","lambda_2_samples","theta_samples","tau_samples","sigma_samples"), file = paste0(dir_path,"/Aniso_samples.Rdata"))
load(paste0(dir_path,"/Aniso_samples.Rdata"))
source("Scripts/Utilities.r")



#lambda_1_samples = exp(lambda_1_samples)
#lambda_2_samples = exp(lambda_2_samples)
theta_samples = pi/2 * rstanarm::invlogit(theta_samples)
#tau_samples = exp(tau_samples)
#sigma_samples = exp(sigma_samples)

## Estimate
lambda_1_estimates = (rowMeans(lambda_1_samples))
lambda_2_estimates = (rowMeans(lambda_2_samples))
theta_estimates = rowMeans(theta_samples)
tau_estimates = (rowMeans(tau_samples))
sigma_estimates = (rowMeans(sigma_samples))
estimates = data.frame(x_1 = gridplot[,1], x_2 = gridplot[,2], mu = rep(mean(samples[,"beta"]), nrow(samples)), sigma = sigma_estimates, tau = tau_estimates,
                       lambda_1 = lambda_1_estimates, lambda_2 = lambda_2_estimates, theta = theta_estimates)
#estimates$Sigma11 = Sigma11(lambda_2_estimates, lambda_1_estimates, pi/2-theta_estimates)
#estimates$Sigma22 = Sigma22(lambda_2_estimates, lambda_1_estimates, pi/2-theta_estimates)
#estimates$Sigma12 = Sigma12(lambda_2_estimates, lambda_1_estimates, pi/2-theta_estimates)
estimates$theta_deg = 90- 180/pi * estimates$theta
x11()
multiple_heatmaps(estimates)
x11()
plot_ellipses(data = estimates)