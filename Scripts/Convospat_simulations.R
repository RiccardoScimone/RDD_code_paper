rm(list = ls())
gc()
library(convoSPAT)
library(ellipse)
source("Scripts/Utilities.r")
a = -1
b = 1
n_knot_grid = 5

knot_lon = seq(a, b, length.out = n_knot_grid+1)
knot_lat = seq(a, b, length.out = n_knot_grid+1)

knot_lon = (knot_lon[1:n_knot_grid] + knot_lon[2:(n_knot_grid+1)])/2
knot_lat = (knot_lat[1:n_knot_grid] + knot_lat[2:(n_knot_grid+1)])/2
knot_coord = expand.grid(knot_lon,knot_lat)
names(knot_coord) = c("x_1", "x_2")
simulation_id = "simulation_4"
dir_path = paste0("Simulations/",simulation_id)
load(paste0(dir_path,"/simulated_process_list.Rdata"))
clusters = 5
cl = makePSOCKcluster(clusters)
registerDoParallel(cl)
writeLines(c(""), paste0(dir_path,"/log.txt"))


NSfit_model = foreach (realization = 1:10,.packages = c("convoSPAT")) %dopar% {
  sink(paste0(dir_path,"/log.txt"),append = T)
 # N = 1000
  set.seed(18011996)
  sim_const = simulated_process_list[[realization]]
  simulated_process = sim_const#[sample(1:nrow(sim_const),N),]
  at = NSconvo_fit(coords = simulated_process[,1:2], data = simulated_process[,3],
                          cov.model = "matern",fit.radius = 0.1, lambda.w = 0.5,
                          mc.locations = knot_coord, ns.mean = T, ns.nugget = T, ns.variance = T,
                          kappa = 2, fix.kappa = T)
  sink()
  at
#  x11()
#  plot(NSfit_model[[realization]],fit.radius = 0.35, xlim = 1.5*c(a,b), ylim = 1.5*c(a,b),asp = 1)

}
saveRDS(file = paste0(dir_path,"/convoSPAT_results.rds"),object = NSfit_model)
stopCluster(cl)


#### Post processing #####
rm(list = ls())
library(convoSPAT)
library(ellipse)
a = -1
b = 1
n_knot_grid = 10
knot_lon = seq(a, b, length.out = n_knot_grid+1)
knot_lat = seq(a, b, length.out = n_knot_grid+1)
knot_lon = (knot_lon[1:n_knot_grid] + knot_lon[2:(n_knot_grid+1)])/2
knot_lat = (knot_lat[1:n_knot_grid] + knot_lat[2:(n_knot_grid+1)])/2
knot_coord = expand.grid(knot_lon,knot_lat)
names(knot_coord) = c("x_1", "x_2")

simulation_id = "simulation_4"
dir_path = paste0("Simulations/",simulation_id)
NSfit_model = readRDS(file = paste0(dir_path,"/convoSPAT_results.rds"))
load(paste0(dir_path,"/simulated_process_list.Rdata"))

N = 10000
parameters_extracted = parameters_smoothed =  list()
for (realization in 1:10)
{
  #### build anisotropy estimates
  lambda_1 = lambda_2 = theta = rep(0,N)
  for ( i in 1:N)
  {
    eig = eigen(NSfit_model[[realization]]$kernel.ellipses[,,i])
    if (acos(eig$vectors[1,1]) > pi/2){
      theta[i] = 180/pi * (acos(eig$vectors[1,1]) - pi/2)
      lambda_1[i] = eig$values[2]
      lambda_2[i] = eig$values[1]
  }
  else
  {
    theta[i] = 180/pi * (acos(eig$vectors[1,1]))
    lambda_1[i] = eig$values[1]
    lambda_2[i] = eig$values[2]
  }
  }
  parameters_extracted[[realization]] = data.frame( 
    x_1 = NSfit_model[[realization]]$coords[,1],
    x_2 = NSfit_model[[realization]]$coords[,2],
    mu = NSfit_model[[realization]]$beta.est[,1],
    sigma = sqrt(NSfit_model[[realization]]$sigmasq.est),
    lambda_1 = lambda_1,
    lambda_2 = lambda_2,
    theta_deg = theta,
    tau = NSfit_model[[realization]]$tausq.est,
    phi = sqrt(lambda_1/lambda_2)
    )
                                                                                                      
}

saveRDS(file = paste0(dir_path,"/convoSPAT_results_estimates.rds"),object = parameters_extracted)
source("Scripts/Utilities.r")

plot.NSconvo(NSfit_model[[8]], xlim = 1.5*c(a,b), ylim = 1.5*c(a,b),asp = 1)

x11()
multiple_heatmaps(parameters_extracted[[6]][,1:6])


### average_analysis
real_pars = readRDS(paste0(dir_path,"/real_pars.rds"))
set.seed(18011996)
real_pars = real_pars[sample(1:nrow(real_pars),N),]
real_pars$theta = NULL
real_pars[,c(7,8)] = real_pars[,c(8,7)]
names(real_pars)[c(7,8)] = names(real_pars)[c(8,7)]
averaged_estimation = parameters_extracted[[1]][,-c(1,2)]
errors = parameters_extracted
for ( j in 1:10)
  errors[[j]] = abs (errors[[j]][,-c(1,2)] - real_pars[,-c(1,2)])

averaged_error_estimation = errors[[1]]

for(j in 2:10)
{
  averaged_estimation = parameters_extracted[[j]][,-c(1,2)] + averaged_estimation
  averaged_error_estimation = averaged_error_estimation + errors[[j]]
}

averaged_error_estimation = averaged_error_estimation/10
averaged_estimation = averaged_estimation/10

averaged_estimation = cbind(real_pars[,1:2],averaged_estimation)
averaged_error_estimation = cbind(real_pars[,1:2],averaged_error_estimation)

averaged_estimation$lambda_1 = log(averaged_estimation$lambda_1)
averaged_estimation$lambda_2 = log(averaged_estimation$lambda_2)

p = multiple_heatmaps(averaged_estimation)
ggsave(filename = paste0(dir_path,"/convoSPAT_multiple",".pdf"),plot = p,width = 30,height = 15,dpi = "retina")

spatial_averaged_errors = NULL

for ( j in 1:10)
  spatial_averaged_errors = rbind(spatial_averaged_errors, colMeans(errors[[j]]))
  
err_data = data.frame(spatial_averaged_errors)
err_data$lambda_1 = log(err_data$lambda_1)
err_data$lambda_2 = log(err_data$lambda_2)
err_data$sigma = log(err_data$sigma)

write_rds(err_data,file = paste0(dir_path,"/convoSPAT_errors.rds"))


for ( j in 1:10){
  pdf(paste0(dir_path,"/plots/convoSPAT_ellipse",j,".pdf"))
  par(cex = 1.5)
  plot.NSconvo(NSfit_model[[j]],fit.radius = 0.35, xlim = 1.5*c(a,b), ylim = 1.5*c(a,b),asp = 1, xlab = "x_1", ylab = "x_2")
  dev.off()
}
