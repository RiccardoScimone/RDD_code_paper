rm(list = ls())
graphics.off()
source("Scripts/Utilities.r")
source("Scripts/NNMatrix.R")
simulation_id = "simulation_1"
dir_path = paste0("Simulations/",simulation_id)
real_pars = readRDS(paste0(dir_path,"/real_pars.rds"))
real_pars$theta = NULL
real_pars = real_pars[,c(6,5,3,4,7)]
parameter_names = names(real_pars)
scaling_functions = NULL
scaling_functions[[1]] = function (x) {return(x)}
scaling_functions[[2]] = scaling_functions[[3]] = scaling_functions[[4]] =  exp
scaling_functions[[5]] = function(x) {return(180/pi*x)}

load(paste0(dir_path,"/simulated_process.Rdata"))
ord = readRDS(paste0(dir_path,"/ord.rds"))
grid = simulated_process[,1:2]

axx <- list(
  gridcolor='rgb(255, 255, 255)',
  zerolinecolor='rgb(255, 255, 255)',
  showbackground=TRUE,
  backgroundcolor='rgb(230, 230,230)'
)
x_1 = unique(grid[,1])
x_2 = unique(grid[,2])
plist = RDD_3d_plots_gp(parameter_names,dir_path,x_1,x_2,ord,real_pars,scaling_functions)
pphi = RDD_3d_plots_phi(dir_path,x_1,x_2,ord,real_pars)
plist[[6]] = pphi
plot = subplot(plist, nrows = 2)
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
htmlwidgets::saveWidget(partial_bundle(plot), file = paste0(dir_path,"/plots/","fitted_parameters_RDD_gp.html"))

