#devtools::install_github("giacomodecarlo/FunctionalLSM",force = T)
rm(list = ls())
graphics.off()
library(LocallyStationaryModels)
library(tidyverse)
library(ggpubr)
library(viridis)
simulation_id = "simulation_4"
dir_path = paste0("Simulations/",simulation_id)
load(paste0(dir_path,"/simulated_process_list.Rdata"))

realization = 4
N = 10000
set.seed(18011996)
sim_const = simulated_process_list[[realization]]
simulated_process = sim_const#[sample(1:nrow(sim_const),N),]

x11()
ggplot(simulated_process) + geom_tile(mapping = aes(x_1,x_2, fill = process)) + coord_fixed() + theme_pubclean(base_size = 30) + 
  scale_fill_viridis(option = "inferno") + theme(legend.text=element_text(size=15), legend.key.size = unit(1.5, 'cm'))


ancors = find_anchorpoints.lsm(dataset = as.matrix(simulated_process[,1:2]),4,plot_output =  T)

vario = variogram.lsm(z = simulated_process$process,d = as.matrix(simulated_process[,1:2]),
                      anchorpoints = ancors$anchorpoints, epsilon = 0.4 , n_angles = 3, n_intervals = 16, kernel_id = "gaussian")

x11()
plotvario(vario,pos = 6)

sol = findsolutions.lsm(vario = vario,id = "exponentialnugget",
                        initial.position = c(0.1,0.1,pi/4,0.3,0.01),remove_not_convergent = T)

x11()
par(mfrow = c(1,3))
res = plot.lsm(sol,ancors,simulated_process$process, as.matrix(simulated_process[,1:2]))


cv.lsm()


