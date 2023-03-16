### Analysis on the posterior distributions of spatial parameters
rm(list = ls())
graphics.off()
source("Scripts/Utilities.r")
expose_stan_functions("Stan_models/to_expose.stan")
simulation_id = "simulation_7"
dir_path = paste0("Simulations/",simulation_id)

K_vec =  c(1,4,8,16,32)
N_vec  = c(1000,2000,5000,8000,10000)
#K_vec =  c(8)
#N_vec  = c(5000)
real_pars = readRDS(paste0(dir_path,"/real_pars.rds"))
load(paste0(dir_path,"/simulated_process_list.Rdata"))


for ( N in N_vec){
  for (K in K_vec){
    if ((N != 1000 | K != 32) & !(N == 2000 & K == 32 ) ){
    multiple_realization_estimation = readRDS(paste0("Simulations/",simulation_id,"/RDD_estimate_mr_",K,N,".rds"))
    dir.create(paste0(dir_path,"/",K,"_",N))
    dir_path = paste0(dir_path,"/",K,"_",N,"/")
    dir.create(paste0(dir_path,"plots"))
    for (realization in 1:10){
      simulated_process = simulated_process_list[[realization]]
      est_spatial_pars = multiple_realization_estimation[[realization]]
      est_spatial_pars = cbind(simulated_process[,1:2],est_spatial_pars)
      est_spatial_pars$theta = NULL
      real_pars$phi = sqrt(real_pars$lambda_1/real_pars$lambda_2)
#    est_spatial_pars$Sigma11 = Sigma11(est_spatial_pars$lambda_1,est_spatial_pars$lambda_2,pi/180*est_spatial_pars$theta_deg)
#    est_spatial_pars$Sigma22 = Sigma22(est_spatial_pars$lambda_1,est_spatial_pars$lambda_2,pi/180*est_spatial_pars$theta_deg)
#    est_spatial_pars$Sigma12 = Sigma12(est_spatial_pars$lambda_1,est_spatial_pars$lambda_2,pi/180*est_spatial_pars$theta_deg)
     est_spatial_pars$lambda_1 = log(est_spatial_pars$lambda_1)
     est_spatial_pars$lambda_2 = log(est_spatial_pars$lambda_2)
    
     p = multiple_heatmaps(est_spatial_pars)
     ggsave(filename = paste0(dir_path,"/plots/est_pars_RDD",realization,".pdf"),plot = p,width = 30,height = 15,dpi = "retina")
     est_spatial_pars$lambda_1 = exp(est_spatial_pars$lambda_1)
     est_spatial_pars$lambda_2 = exp(est_spatial_pars$lambda_2)
     pdf(paste0(dir_path,"/plots/EllipsesRDD",realization,".pdf"))
     plot(plot_ellipses(est_spatial_pars))
     dev.off()
     print(realization)
    }
    dir_path = paste0("Simulations/",simulation_id)
    }
}
}

### Average Analysis with multiple realizations
### Load and process data
real_pars = readRDS(paste0(dir_path,"/real_pars.rds"))
load(paste0(dir_path,"/simulated_process.Rdata"))

dir_path = paste0("Simulations/",simulation_id,"/multiple")
dir.create(dir_path)
real_pars$phi = sqrt(real_pars$lambda_1/real_pars$lambda_2)

real_pars[,c(7,9)] = real_pars[,c(9,7)]
names(real_pars)[c(7,9)] = names(real_pars)[c(9,7)]

multiple_realization_estimation = list()
averaged_estimation = list()
errors = list()
averaged_error_estimation = list()
for ( N in as.character(N_vec)){
  multiple_realization_estimation[[N]] = list()
  averaged_estimation[[N]] = list()
  errors[[N]] = list()
  averaged_error_estimation[[N]] = list()
  for (K in as.character(K_vec)){
    if ((N != "1000" | K != "32")& !(N == "2000" & K=="32" )){
      multiple_realization_estimation[[N]][[K]] = readRDS(paste0("Simulations/",simulation_id,"/RDD_estimate_mr_",K,N,".rds"))
      errors[[N]][[K]] = multiple_realization_estimation[[N]][[K]]
      for ( j in 1:10)
        errors[[N]][[K]][[j]] = abs (errors[[N]][[K]][[j]] - real_pars[,-c(1,2)])
      
      averaged_estimation[[N]][[K]] = multiple_realization_estimation[[N]][[K]][[1]]
      averaged_error_estimation[[N]][[K]] = errors[[N]][[K]][[1]]
      
      for (j in 2:10){
        averaged_estimation[[N]][[K]] = averaged_estimation[[N]][[K]] + multiple_realization_estimation[[N]][[K]][[j]]
        averaged_error_estimation[[N]][[K]] = averaged_error_estimation[[N]][[K]] + errors[[N]][[K]][[j]]
      }
      
      
      averaged_estimation[[N]][[K]] = averaged_estimation[[N]][[K]]/10
      averaged_error_estimation[[N]][[K]] = averaged_error_estimation[[N]][[K]]/10
      averaged_estimation[[N]][[K]] = cbind(simulated_process[,1:2],averaged_estimation[[N]][[K]])
      averaged_error_estimation[[N]][[K]] = cbind(simulated_process[,1:2],averaged_error_estimation[[N]][[K]])
      
      averaged_estimation[[N]][[K]]$lambda_1 = log(averaged_estimation[[N]][[K]]$lambda_1)
      averaged_estimation[[N]][[K]]$lambda_2 = log(averaged_estimation[[N]][[K]]$lambda_2)
 #     p = multiple_heatmaps(averaged_estimation[[N]][[K]])
 #     ggsave(filename = paste0(dir_path,"/est_pars_RDD_multiple",N,"_",K,".pdf"),plot = p,width = 30,height = 15,dpi = "retina")
      
    }
  }
  print(N)
  
}

### Assemble error informations


spatial_averaged_errors = NULL
keys = NULL


for (N in as.character(N_vec))
  for(K in as.character(K_vec))
    if ((N != "1000" | K != "32")& !(N == "2000" & K=="32" ))
    for (j in 1:10)
    {
      keys = rbind(keys,c(N,K))
      spatial_averaged_errors = rbind(spatial_averaged_errors, colMeans(errors[[N]][[K]][[j]]))
    }


err_data = data.frame(keys)
err_data = cbind(err_data,spatial_averaged_errors)
names(err_data)[1:2] = c("samples","K")
#err_data$lambda_1 = log(err_data$lambda_1)
#err_data$lambda_2 = log(err_data$lambda_2)
#err_data$sigma = log(err_data$sigma)
#err_data$nugget = log(err_data$nugget)
err_data_pivoted = err_data %>% pivot_longer(cols = 3:10, names_to = "pars", values_to = "mean_spatial_errs")


err_data_summarises = err_data_pivoted %>% group_by(samples,K,pars) %>% summarise_all(list(average = mean,sd = sd))

err_data_summarises$K = as.numeric(err_data_summarises$K)
err_data_summarises$samples = factor(err_data_summarises$samples, levels = as.character(N_vec))
  
err_data_summarises = err_data_summarises %>% filter(pars != "theta")


errors_convo = read_rds("Simulations/simulation_7/convoSPAT_errors.rds")
#errors_convo$tau = log(errors_convo$tau)
errors_convo$nugget = errors_convo$tau
errors_convo = errors_convo %>% mutate_at(c("lambda_1","lambda_2","sigma"), exp)

errors_convo$tau = NULL
errors_convo = errors_convo %>% pivot_longer(cols = 1:7, names_to = "pars", values_to = "mean_spatial_errs")
errors_convo = errors_convo %>% group_by(pars) %>% summarise_all(list(average = mean, sd = sd))


errors_Fou = read_rds("Simulations/simulation_7/Fouedjo_multiple_errors.rds") 
errors_Fou$nugget = errors_Fou$tau
errors_Fou$tau = NULL
errors_Fou = errors_Fou %>% mutate_at(c("lambda_1","lambda_2","sigma","nugget"), exp)

errors_Fou = errors_Fou %>% pivot_longer(cols = 1:7, names_to = "pars", values_to = "mean_spatial_errs")
errors_Fou = errors_Fou %>% group_by(pars) %>% summarise_all(list(average = mean, sd = sd))

plot = ggplot(err_data_summarises, aes(x = K, y = log(average), colour = samples)) + geom_line() + geom_point() + 
   #coord_trans(y = "log") +
  scale_x_continuous(breaks = c(1,4,8,16,32)) + 
  geom_errorbar(aes(ymin = log( average - sd/sqrt(10)), ymax = log(average + sd/sqrt(10))), width = .2) +
  geom_hline(data = errors_convo,aes(yintercept = log(average) )) +
  geom_hline(data = errors_convo,aes(yintercept = log(average + sd/sqrt(10)) ),linetype="dashed") +  
  geom_hline(data = errors_convo,aes(yintercept = log(average - sd/sqrt(10)) ),linetype="dashed")+ 
  geom_hline(data = errors_Fou,aes(yintercept = log(average) ), col = "red") + geom_hline(data = errors_Fou,aes(yintercept = log(average + sd/sqrt(10))),linetype="dashed",col = "red") +  
  geom_hline(data = errors_Fou,aes(yintercept = log( average - sd/sqrt(10)) ),linetype="dashed", col = "red")+
  labs(y = "Logarithm of average error") + 
  facet_wrap(~pars,scales = "free") + theme_pubclean(base_size = 30)

ggsave(filename = paste0("Paper_plots/errs_curves",simulation_id,".pdf"),plot = plot, width = 22,height = 15,dpi = "retina")


