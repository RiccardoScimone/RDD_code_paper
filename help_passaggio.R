rm(list = ls())
source("Scripts/Utilities.r")

rdd = read_rds("Simulations/simulation_7/RDD_fitting_mr_1610000.rds")

load(paste0("Simulations/simulation_7/simulated_process_list.Rdata"))

plist = list()
for ( j in 1:10)
{
  to_plot = data.frame(simulated_process_list[[1]][,1:2], sigma = rdd[[1]][,2,j])
  x11()
  plot(multiple_heatmaps(to_plot))
}


graphics.off()


rdd = read_rds("Simulations/simulation_7/RDD_estimate_mr_1610000.rds")
for ( j in 1:10)
{
  to_plot = data.frame(simulated_process_list[[1]][,1:2],rdd[[j]][,1:2])
  x11()
  plot(multiple_heatmaps(to_plot))
}