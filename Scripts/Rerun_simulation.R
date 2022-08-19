rm(list = ls())
source("Scripts/Utilities.r")

## Specify Simulation ID

sim_id = "simulation_1"
dir.create(paste0("Simulations/",sim_id))
dir.create(paste0("Simulations/",sim_id,"/plots"))


spatial_parameters = readRDS(paste0("Simulations/",sim_id,"/real_pars.rds"))
spatial_parameters$phi = sqrt(spatial_parameters$lambda_1/spatial_parameters$lambda_2)
#### Create plot with real parameter functions
spatial_parameters$theta = NULL
spatial_parameters = spatial_parameters[,c(1,2,6,5,3,4,7,8)]
p = multiple_heatmaps(spatial_parameters)
ggsave(filename = paste0("Simulations/",sim_id,"/plots/real_pars.pdf"),plot = p,width = 30,height = 15,dpi = "retina")
