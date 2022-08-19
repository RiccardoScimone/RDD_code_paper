### Analysis on the posterior distributions of spatial parameters
rm(list = ls())
graphics.off()
source("Scripts/Utilities.r")
expose_stan_functions("Stan_models/to_expose.stan")
simulation = T
if (simulation){
  simulation_id = "simulation_4"
  dir_path = paste0("Simulations/",simulation_id)
}

if(!simulation)
  dir_path = "GHCN"

K = 32
N = 5000

dir.create(paste0(dir_path,"/",K,"_",N))


if(!simulation)
  predgrid = readRDS(paste0(dir_path,"/predgrid.rds"))

if(simulation){
  real_pars = readRDS(paste0(dir_path,"/real_pars.rds"))
  load(paste0(dir_path,"/simulated_process.Rdata"))
#  set.seed(18011996)
#  simulated_process = simulated_process[sample(1:nrow(simulated_process),N),]
}

fit = readRDS(paste0(dir_path,"/RDD_fitting",K,N,".rds"))
est_spatial_pars = readRDS(paste0(dir_path,"/RDD_estimate",K,N,".rds"))
dir_path = paste0(dir_path,"/",K,"_",N,"/")
dir.create(paste0(dir_path,"plots"))

if(!simulation){
  est_spatial_pars_obs_loc = est_spatial_pars[(nrow(predgrid)+1):nrow(est_spatial_pars),]
  est_spatial_pars = est_spatial_pars[1:nrow(predgrid),]
#  obs = read.csv("Scripts/Colorado_data.csv",) %>% select(Longitude, Latitude, logPrecip)
  obs = read.csv("GHCN/average_dataset_risser.csv") %>% select(LON, LAT, logprep)
  
  names(obs)[1:2] = c("x_1","x_2")
  est_spatial_pars_obs_loc = cbind(est_spatial_pars_obs_loc,obs[,1:2])
}


if(simulation)
  est_spatial_pars = cbind(simulated_process[,1:2],est_spatial_pars)

if(!simulation){
  est_spatial_pars = cbind(predgrid[,1:2],est_spatial_pars)
  names(est_spatial_pars)[1:2] = c("x_1","x_2")
}
est_spatial_pars$theta = NULL
if(simulation)
  real_pars$phi = sqrt(real_pars$lambda_1/real_pars$lambda_2)


est_spatial_pars$Sigma11 = Sigma11(est_spatial_pars$lambda_1,est_spatial_pars$lambda_2,pi/180*est_spatial_pars$theta_deg)
est_spatial_pars$Sigma22 = Sigma22(est_spatial_pars$lambda_1,est_spatial_pars$lambda_2,pi/180*est_spatial_pars$theta_deg)
est_spatial_pars$Sigma12 = Sigma12(est_spatial_pars$lambda_1,est_spatial_pars$lambda_2,pi/180*est_spatial_pars$theta_deg)

est_spatial_pars$lambda_1 = log(est_spatial_pars$lambda_1)
est_spatial_pars$lambda_2 = log(est_spatial_pars$lambda_2)

p = multiple_heatmaps(est_spatial_pars)
ggsave(filename = paste0(dir_path,"plots/est_pars_RDD.pdf"),plot = p,width = 30,height = 15,dpi = "retina")



## Some more meaningful map about elevation
colorado = FALSE
if(colorado){
  grid_plot = fread("Scripts/Colorado_grid.csv")
  p_elev = ggplot(grid_plot, aes(x = longitude, y = latitude, fill = elevation)) + geom_tile() +
   scale_fill_viridis(direction=1, option = "turbo") + theme_pubclean(base_size = 30)+ 
     theme(legend.text=element_text(size=15), legend.key.size = unit(1.5, 'cm'))+ coord_fixed() 
  ggsave(filename = paste0(dir_path,"plots/elevation.pdf"),plot = p_elev,width = 30,height = 15,dpi = "retina")
  }

ghcn = FALSE
if(ghcn){
  grid_plot = fread("GHCN/ghcn_grid.csv")
  p_elev = ggplot(grid_plot, aes(x = longitude, y = latitude, fill = elevation)) + geom_tile() +
    scale_fill_viridis(direction=1, option = "turbo") + theme_pubclean(base_size = 30)+ 
    theme(legend.text=element_text(size=15), legend.key.size = unit(1.5, 'cm'))+ coord_fixed() 
  ggsave(filename = paste0(dir_path,"plots/elevation.pdf"),plot = p_elev,width = 30,height = 15,dpi = "retina")

}

est_spatial_pars$lambda_1 = exp(est_spatial_pars$lambda_1)
est_spatial_pars$lambda_2 = exp(est_spatial_pars$lambda_2)
pdf(paste0(dir_path,"plots/EllipsesRDD.pdf"),width = 30,height = 15)
plot(plot_ellipses(est_spatial_pars))
dev.off()

if(!simulation){ # go with kriging
  Locations = rbind ( c (est_spatial_pars_obs_loc$x_1, est_spatial_pars$x_1) , c (est_spatial_pars_obs_loc$x_2, est_spatial_pars$x_2) )
  lambda_1 =  c (est_spatial_pars_obs_loc$lambda_1, est_spatial_pars$lambda_1)
  lambda_2 =  c (est_spatial_pars_obs_loc$lambda_2, est_spatial_pars$lambda_2)
  theta =  c (est_spatial_pars_obs_loc$theta, est_spatial_pars$theta_deg*pi/180)
  Anis = Compute_Aniso(lambda_1,lambda_2, theta)
  sigma = c (est_spatial_pars_obs_loc$sigma, est_spatial_pars$sigma)
  C = matern_ns_corr_mat(Locations, Anis, lambda_1*lambda_2, 5/2 , sigma)
  library(SpatialTools)
  N_obs = nrow(obs)
  for ( i in 1:N_obs) C[i,i] = C[i,i] + 0.125
  N_pred = nrow(est_spatial_pars)
#  krig = krige.ok(obs$logprep, C[1:N_obs,1:N_obs],C[(N_obs + 1) : (N_pred + N_obs), (N_obs + 1) : (N_pred + N_obs)], C[1:N_obs,(N_obs + 1) : (N_pred + N_obs)])
  krig = krige.sk(obs$logprep - est_spatial_pars_obs_loc$mu, C[1:N_obs,1:N_obs],C[(N_obs + 1) : (N_pred + N_obs), (N_obs + 1) : (N_pred + N_obs)], C[1:N_obs,(N_obs + 1) : (N_pred + N_obs)], m = 0)
  
  krig_to_plot = data.frame(x_1 = est_spatial_pars$x_1,x_2 = est_spatial_pars$x_2, ord = krig$pred + est_spatial_pars$mu, var = krig$mspe, nomean = krig$pred)
  krig_to_plot$preds = trim(krig_to_plot$ord,-2.4,2.5)
  p = ggplot(krig_to_plot) + geom_tile(aes(x = x_1,y = x_2, fill = preds)) + scale_fill_viridis(direction=-1,breaks = seq(-2.3,2.4,by = 0.7)) + theme_pubclean(base_size = 30)+ theme(legend.text=element_text(size=15), legend.key.size = unit(1.5, 'cm'))+ coord_fixed()
  ggsave(filename = paste0(dir_path,"plots/skRDD.pdf"),plot = p,width = 30,height = 15,dpi = "retina")
  ggsave(filename = paste0(dir_path,"plots/skRDDellipse.pdf"),plot = add_ellipses(est_spatial_pars,p,"red"),width = 30,height = 15,dpi = "retina")
  
  p = ggplot(krig_to_plot) + geom_tile(aes(x = x_1,y = x_2, fill = nomean)) + scale_fill_viridis(direction=-1) + theme_pubclean(base_size = 30)+ theme(legend.text=element_text(size=15), legend.key.size = unit(1.5, 'cm'))+ coord_fixed()
  ggsave(filename = paste0(dir_path,"plots/skRDDnomean.pdf"),plot = p,width = 30,height = 15,dpi = "retina")
  
  x11()
  ggplot(krig_to_plot) + geom_tile(aes(x = x_1,y = x_2, fill = sqrt(var))) + scale_fill_viridis(direction=-1, option = "inferno")
}




















if (simulation){
error_funs = list()
error_funs[["mu"]] = abs_err

for ( name in names(est_spatial_pars)[-(1:3)]) error_funs[[name]] = relative_err

log_rescale = c(FALSE,FALSE,TRUE,TRUE,FALSE,FALSE)

errs = compute_errors(data = est_spatial_pars,real = real_pars, error_funs = error_funs,log_rescale = log_rescale)
p = multiple_heatmaps(errs)
ggsave(filename = paste0(dir_path,"/plots/errs_RDD.pdf"),plot = p,width = 30,height = 15,dpi = "retina")



p = draw_pars_histograms(errs[3:8])
ggsave(filename = paste0(dir_path,"/plots/errs_RDD_hists.pdf"),plot = p,width = 30,height = 15,dpi = "retina")

#### 3D plotting
axx <- list(
  gridcolor='rgb(255, 255, 255)',
  zerolinecolor='rgb(255, 255, 255)',
  showbackground=TRUE,
  backgroundcolor='rgb(230, 230,230)'
)
ply_list = RDD_3d_plots(est_spatial_pars, real_pars)
plot = subplot(ply_list, nrows = 2)
plot = plot  %>% layout(title = "Fitted Parameters, Gaussian processes",
                 scene = list(domain=list(x=c(0,0.3),y=c(0.5,1)),
                              xaxis=axx, yaxis=axx, zaxis=axx,
                              aspectmode='cube'),
                 scene2 = list(domain=list(x=c(0.3,0.6),y=c(0.5,1)),
                               xaxis=axx, yaxis=axx, zaxis=axx,
                               aspectmode='cube'),
                 scene3 = list(domain=list(x=c(0.6,0.9),y=c(0.5,1)),
                               xaxis=axx, yaxis=axx, zaxis=axx,
                               aspectmode='cube'),
                 scene4 = list(domain=list(x=c(0,0.3),y=c(0,0.5)),
                              xaxis=axx, yaxis=axx, zaxis=axx,
                              aspectmode='cube'),
                 scene5 = list(domain=list(x=c(0.3,0.6),y=c(0,0.5)),
                               xaxis=axx, yaxis=axx, zaxis=axx,
                               aspectmode='cube'),
                 scene6 = list(domain=list(x=c(0.6,0.9),y=c(0,0.5)),
                               xaxis=axx, yaxis=axx, zaxis=axx,
                               aspectmode='cube')
                                 )
htmlwidgets::saveWidget(partial_bundle(plot), file = paste0(dir_path,"/plots/","fitted_parameters_RDD.html"))
}