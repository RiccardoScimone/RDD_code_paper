# Clean the environment
rm(list = ls())
graphics.off()
# Load the libraries
library(LocallyStationaryModels)
source("Scripts/Utilities.r")
n_knot_grid = 8

simulation_id = "simulation_7"
dir_path = paste0("Simulations/",simulation_id)
load(paste0(dir_path,"/simulated_process_list.Rdata"))
Fouedjofit_model = list()
Fouedjofit_anchors = list()


set.seed(18011996)
for (realization in 1:10){
#  N = 10000
#  set.seed(18011996)
  sim_const = simulated_process_list[[realization]]
  simulated_process = sim_const#[sample(1:nrow(sim_const),N),]
  d = cbind(simulated_process[,1],simulated_process[,2])
  a = find_anchorpoints.lsm(d,n = n_knot_grid,plot_output = F)
  z = as.matrix(simulated_process[,3])
  vario = variogram.lsm(z = z,d = d,a = a$anchorpoints,epsilon = 0.25,n_angles=6,n_intervals=16,dim=1,kernel_id="gaussian")
  solu = findsolutions.lsm(vario,lower.delta = 0.01, "exponentialnugget", c(0.01,0.01,0.4,3,0.3), 
                           upper.bound = c(0.5,0.5,pi/2,10,1),
                           remove_not_convergent = T)

  pars = smooth.lsm(solu,newpoints = as.matrix(simulated_process[,1:2]))
  pred = predict.lsm(solu,as.matrix(simulated_process[,1:2]),plot_output = F,predict_y = F)
  pars_df = cbind(simulated_process[,1:2],pars[[1]],pred$smoothed_means[,1])
 # names(pars_df) = c("x_1","x_2","lambda_1","lambda_2","theta","sigma","tau^2")
#  pars_df$lambda_1 = log(pars_df$lambda_1^2)
#  pars_df$lambda_2 = log(pars_df$lambda_2^2)
  Fouedjofit_model[[realization]] = pars_df
  Fouedjofit_anchors[[realization]] = solu
}

saveRDS(file = paste0(dir_path,"/Fouedjo_anchors.rds"),object = Fouedjofit_anchors)
saveRDS(file = paste0(dir_path,"/Fouedjo_results.rds"),object = Fouedjofit_model)

Fouedjofit_model = readRDS(paste0(dir_path,"/Fouedjo_results.rds"))
Fou_anchors = readRDS(paste0("Simulations/",simulation_id,"/Fouedjo_anchors.rds"))

for(j in 1:10) { 
  names(Fouedjofit_model[[j]]) = c("x_1","x_2","lambda_1","lambda_2","theta_deg","sigma","tau","mu")
  Fouedjofit_model[[j]]$theta_deg = 180/pi*Fouedjofit_model[[j]]$theta_deg
  Fouedjofit_model[[j]]$phi = Fouedjofit_model[[j]]$lambda_1/ Fouedjofit_model[[j]]$lambda_2
  Fouedjofit_model[[j]]$lambda_1 = Fouedjofit_model[[j]]$lambda_1^2
  Fouedjofit_model[[j]]$lambda_2 = Fouedjofit_model[[j]]$lambda_2^2
  
#  Fouedjofit_model[[j]]$lambda_1 = exp(Fouedjofit_model[[j]]$lambda_1 )
#  Fouedjofit_model[[j]]$lambda_2 = exp(Fouedjofit_model[[j]]$lambda_2 )
  Fouedjofit_model[[j]] = Fouedjofit_model[[j]][,c(1,2,8,6,3,4,5,7,9)]
}
saveRDS(file = paste0(dir_path,"/Fouedjo_results.rds"),object = Fouedjofit_model)


#x11()
#ggplot(simulated_process_list[[1]]) + geom_tile(mapping = aes(x_1,x_2, fill = process)) + coord_fixed() + theme_pubclean(base_size = 30) + 
#  scale_fill_viridis(option = "inferno") + theme(legend.text=element_text(size=15), legend.key.size = unit(1.5, 'cm'))

dir.create(paste0(dir_path,"/Fouedjo_multiple"))
for ( i in 1:10){
  temp_anchors_df = data.frame(cbind(Fou_anchors[[i]]$anchorpoints,Fou_anchors[[i]]$solutions[,1:3]))
  names(temp_anchors_df) = c("x_1","x_2","lambda_1","lambda_2","theta_deg")
  temp_anchors_df$theta_deg = temp_anchors_df$theta_deg * 180/pi
  temp_anchors_df$lambda_1 = temp_anchors_df$lambda_1^2
  temp_anchors_df$lambda_2 = temp_anchors_df$lambda_2^2
  
  p = plot_ellipses(Fouedjofit_model[[i]],l_out = 8,scale = 0.5)
  ggsave(filename = paste0(dir_path,"/Fouedjo_multiple/ellipses_",i,".pdf"),plot = p,width = 20,height = 15,dpi = "retina")
  p = plot_ellipses(temp_anchors_df,l_out = 8,scale = 0.5)
  ggsave(filename = paste0(dir_path,"/Fouedjo_multiple/ellipses_anchors_",i,".pdf"),plot = p,width = 20,height = 15,dpi = "retina")
}

for ( i in 1:10){
Fouedjofit_model[[i]]$lambda_1 = log(Fouedjofit_model[[i]]$lambda_1)
Fouedjofit_model[[i]]$lambda_2 = log(Fouedjofit_model[[i]]$lambda_2)
p = multiple_heatmaps(Fouedjofit_model[[i]])
ggsave(filename = paste0(dir_path,"/Fouedjo_multiple/estimates_",i,".pdf"),plot = p,width = 30,height = 15,dpi = "retina")
Fouedjofit_model[[i]]$lambda_1 = exp(Fouedjofit_model[[i]]$lambda_1)
Fouedjofit_model[[i]]$lambda_2 = exp(Fouedjofit_model[[i]]$lambda_2)
}




real_pars = readRDS(paste0(dir_path,"/real_pars.rds"))
real_pars$theta = NULL
real_pars[,c(7,8)] = real_pars[,c(8,7)]
names(real_pars)[c(7,8)] = names(real_pars)[c(8,7)]
p = plot_ellipses(real_pars)
ggsave(filename = paste0(dir_path,"/plots/real_ellipses",".pdf"),plot = p,width = 20,height = 15,dpi = "retina")
#real_pars$mu = NULL

averaged_estimation = Fouedjofit_model[[1]][,-c(1,2)]
errors = Fouedjofit_model
for ( j in 1:10)
  errors[[j]] = abs (errors[[j]][,-c(1,2)] - real_pars[,-c(1,2)])

averaged_error_estimation = errors[[1]]

for(j in 2:10)
{
  averaged_estimation = Fouedjofit_model[[j]][,-c(1,2)] + averaged_estimation
  averaged_error_estimation = averaged_error_estimation + errors[[j]]
}

averaged_error_estimation = averaged_error_estimation/10
averaged_estimation = averaged_estimation/10

averaged_estimation = cbind(Fouedjofit_model[[1]][,c(1,2)],averaged_estimation)
averaged_error_estimation = cbind(real_pars[,1:2],averaged_error_estimation)

averaged_estimation$lambda_1 = log(averaged_estimation$lambda_1)
averaged_estimation$lambda_2 = log(averaged_estimation$lambda_2)

p = multiple_heatmaps(averaged_estimation)
ggsave(filename = paste0(dir_path,"/Fouedjo_multiple",".pdf"),plot = p,width = 30,height = 15,dpi = "retina")

spatial_averaged_errors = NULL

for ( j in 1:10)
  spatial_averaged_errors = rbind(spatial_averaged_errors, colMeans(errors[[j]]))

err_data = data.frame(spatial_averaged_errors)
err_data$lambda_1 = log(err_data$lambda_1)
err_data$lambda_2 = log(err_data$lambda_2)
err_data$sigma = log(err_data$sigma)
err_data$tau = log(err_data$tau)

write_rds(err_data,file = paste0(dir_path,"/Fouedjo_multiple_errors.rds"))
