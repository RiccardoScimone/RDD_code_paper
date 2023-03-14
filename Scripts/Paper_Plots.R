rm(list = ls())
graphics.off()
#### Figure 1 ###
source("Scripts/Utilities.r")
sim_vec = c("simulation_0","simulation_4","simulation_7")
realiz_list = p_list =  list()
for (i in 1:length(sim_vec))
{
  load(paste0("Simulations/",sim_vec[i],"/simulated_process.rdata"))
  real_pars = readRDS(paste0("Simulations/",sim_vec[i],"/real_pars.rds"))
  realiz_list[[i]] = simulated_process
  p_list[[i]] = ggplot(simulated_process) + geom_tile(mapping = aes(x_1,x_2, fill = process)) + coord_fixed() +
    theme_pubclean(base_size = 40) + 
    scale_fill_viridis(option = "inferno") + theme(legend.text=element_text(size=30), 
                                                   legend.key.size = unit(1.8, 'cm'),legend.title = element_blank(),
                                                   legend.position = "top",
                                                   axis.title = element_blank() )
#  if(i > 1)
#  p_list[[i]] = add_ellipses(real_pars,p_list[[i]])
}

plot = ggarrange(plotlist =  p_list, nrow = 1)
ggsave(filename = "Paper_plots/exprocess.pdf",plot = plot, width = 18,height = 8,dpi = "retina")


## Figure 2###
p = ggplot(simulated_process) + geom_tile(mapping = aes(x_1,x_2, fill = process)) + coord_fixed() + theme_pubclean(base_size = 40) + 
  scale_fill_viridis(option = "inferno", limits = c(-4.8,4.8), breaks = seq(-4.8,4.8,length.out = 5)) + theme(legend.text=element_text(size=30), 
                                                 legend.key.size = unit(1.8, 'cm'),legend.title = element_blank(),
                                                 legend.position = "right",
                                                 axis.title = element_blank() )
p = add_ellipses(real_pars,p)
ggsave(filename = "Paper_plots/nsprocessdetail.pdf",plot = p, width = 18,height = 8,dpi = "retina")
knitr::plot_crop("Paper_plots/nsprocessdetail.pdf")

### Plot realizations (Fig 3)

sim_vec = c("simulation_4","simulation_7")

for ( i in 1:length(sim_vec))
{
  load(paste0("Simulations/",sim_vec[i],"/simulated_process_list.rdata"))
  p_list =  list()
  lims = range(c(simulated_process_list[[1]]$process,simulated_process_list[[4]]$process,simulated_process_list[[8]]$process))
  for ( j in c(1,4,8) ){
    p_list[[j]] = ggplot(simulated_process_list[[j]]) + geom_tile(mapping = aes(x_1,x_2, fill = process)) + coord_fixed() + theme_pubclean(base_size = 40) + 
      scale_fill_viridis(option = "inferno", breaks = round(0.9*seq(lims[1],lims[2], length.out = 5),1)) + theme(legend.text=element_text(size=30), 
                                                     legend.key.size = unit(1.8, 'cm'),legend.title = element_blank(),
                                                     legend.position = "top",
                                                     axis.title = element_blank() )
  }
  p_list = p_list[-which(sapply(p_list, is.null))]
  plot = ggarrange(plotlist =  p_list, nrow = 1, common.legend = T)
  
  ggsave(filename = paste0("Paper_plots/realizations",sim_vec[i],".pdf"),plot = plot, width = 18,height = 8,dpi = "retina")
  
  
}


### Plot summary of real pars

sim_vec = c("simulation_7")

for (i in 1:length(sim_vec))
{
  real_pars = readRDS(paste0("Simulations/",sim_vec[i],"/real_pars.rds"))
  real_pars_short = real_pars %>% select(x_1,x_2,mu,sigma)
  p = multiple_heatmaps(real_pars_short)
  ellipses = plot_ellipses(real_pars,6)
  plot = ggarrange(p,ellipses, nrow = 1, widths = c(2,1),align = "h")
  ggsave(filename = paste0("Paper_plots/real_pars",sim_vec[i],".pdf"),plot = plot, width = 18,height = 8,dpi = "retina")
  knitr::plot_crop(paste0("Paper_plots/real_pars",sim_vec[i],".pdf"))
}



## Plot Fouedjo results for simulation 7

for (i in 1:length(sim_vec))
{
  Fou_est = readRDS(paste0("Simulations/",sim_vec[i],"/Fouedjo_results.rds"))
  Fou_anchors = readRDS(paste0("Simulations/",sim_vec[i],"/Fouedjo_anchors.rds"))
  p_list = ellipses = ellipses_anchors = list()
  for ( j in c(3,9,6) ){
    temp_anchors_df = data.frame(cbind(Fou_anchors[[j]]$anchorpoints,Fou_anchors[[j]]$solutions[,1:3]))
    names(temp_anchors_df) = c("x_1","x_2","lambda_1","lambda_2","theta_deg")
    temp_anchors_df$theta_deg = temp_anchors_df$theta_deg * 180/pi
    temp_anchors_df$lambda_1 = temp_anchors_df$lambda_1^2
    temp_anchors_df$lambda_2 = temp_anchors_df$lambda_2^2
    Fou_est_short = Fou_est[[j]] %>% select(x_1,x_2,sigma)
    p_list[[j]] = multiple_heatmaps(Fou_est_short)
    ellipses[[j]] = plot_ellipses(Fou_est[[j]],8,0.5) + ggtitle("")
    ellipses_anchors[[j]] = plot_ellipses(temp_anchors_df,8,0.5) + ggtitle("")
  }
  p_list = p_list[-which(sapply(p_list, is.null))]
  ellipses = ellipses[-which(sapply(ellipses, is.null))]
  ellipses_anchors = ellipses_anchors[-which(sapply(ellipses_anchors, is.null))]
  
  
  plot = ggarrange(plotlist = p_list, nrow = 1)
  ggsave(filename = paste0("Paper_plots/Fou_pars_sigma_ex",sim_vec[i],".pdf"),plot = plot, width = 18,height = 8,dpi = "retina")
  knitr::plot_crop(paste0("Paper_plots/Fou_pars_mu_sigma",sim_vec[i],".pdf"))
  
  plot = ggarrange(plotlist = ellipses, nrow = 1,align = "h")
  ggsave(filename = paste0("Paper_plots/Fou_ellipses_ex",sim_vec[i],".pdf"),plot = plot, width = 18,height = 8,dpi = "retina")
  knitr::plot_crop(paste0("Paper_plots/Fou_ellipses_ex",sim_vec[i],".pdf"))
  
  plot = ggarrange(plotlist = ellipses_anchors, nrow = 1,align = "h")
  ggsave(filename = paste0("Paper_plots/Fou_ellipses_exanchors",sim_vec[i],".pdf"),plot = plot, width = 18,height = 8,dpi = "retina")
  knitr::plot_crop(paste0("Paper_plots/Fou_ellipses_exanchors",sim_vec[i],".pdf"))
}


## Plot convospat results for simulation 7

for (i in 1:length(sim_vec))
{
  Fou_est = readRDS(file = paste0("Simulations/",sim_vec[i],"/convoSPAT_results_estimates.rds"))
  p_list = ellipses = list()
  for ( j in c(1,4,8) ){
    Fou_est_short = Fou_est[[j]] %>% select(x_1,x_2,sigma)
    p_list[[j]] = multiple_heatmaps(Fou_est_short)
    ellipses[[j]] = plot_ellipses(Fou_est[[j]],8,0.5) + ggtitle("")
  }
  p_list = p_list[-which(sapply(p_list, is.null))]
  ellipses = ellipses[-which(sapply(ellipses, is.null))]
  
  plot = ggarrange(plotlist = p_list, nrow = 1)
  ggsave(filename = paste0("Paper_plots/convo_pars_sigma_ex",sim_vec[i],".pdf"),plot = plot, width = 18,height = 8,dpi = "retina")
  knitr::plot_crop(paste0("Paper_plots/convo_pars_sigma_ex",sim_vec[i],".pdf"))
  plot = ggarrange(plotlist = ellipses, nrow = 1,align = "h")
  ggsave(filename = paste0("Paper_plots/convo_ellipses_ex",sim_vec[i],".pdf"),plot = plot, width = 18,height = 8,dpi = "retina")
  knitr::plot_crop(paste0("Paper_plots/convo_ellipses_ex",sim_vec[i],".pdf"))
}


## Convospat ellipses in anchorpoints
loc_ellips = list()
for (i in 1:length(sim_vec))
{
  convo = readRDS(file = paste0("Simulations/",sim_vec[i],"/convoSPAT_results.rds"))
  ### Build local elliplses
  N = 100
  for (realization in 1:10){
  lambda_1 = lambda_2 = theta = rep(0,N)
  for ( i in 1:N)
  {
    eig = eigen(convo[[realization]]$mc.kernels[,,i])
    if (acos(eig$vectors[1,1]) > pi/2)
      {
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
    loc_ellips[[realization]] = 
      data.frame(x_1 = convo[[realization]]$mc.locations$x_1,x_2 = convo[[realization]]$mc.locations$x_2,lambda_1 = lambda_1, lambda_2 = lambda_2, theta_deg = theta)
  }
  }
}

ellipses = list()

  for ( j in c(1,4,8) ){
    ellipses[[j]] = plot_ellipses(loc_ellips[[j]],10,0.2) + ggtitle("")
  }
  ellipses = ellipses[-which(sapply(ellipses, is.null))]
  
  plot = ggarrange(plotlist = ellipses, nrow = 1,align = "h")
  ggsave(filename = paste0("Paper_plots/convo_ellipses_ex_anchor",sim_vec[1],".pdf"),plot = plot, width = 18,height = 8,dpi = "retina")
  knitr::plot_crop(paste0("Paper_plots/convo_ellipses_ex_anchor",sim_vec[1],".pdf"))
  