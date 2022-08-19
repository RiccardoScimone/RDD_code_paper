rm(list = ls())
source("Scripts/Utilities.r")

## Specify Simulation ID

sim_id = "simulation_7"
dir.create(paste0("Simulations/",sim_id))
dir.create(paste0("Simulations/",sim_id,"/plots"))
seed = 18011996

### Generate grid
sq_N = 100
N = sq_N^2
a = -1
b = 1
grid = expand.grid(seq(a,b,length.out = sq_N),seq(a,b,length.out = sq_N))
names(grid) = c("x_1","x_2")

## Set nu
nu = 2

### Specify spatial functions for parameters
lambda_1 = function(x,y) {return(exp ( (-8 + x)) )}
lambda_2 = function(x,y) {return(exp ( (-5)) )}
theta    = function(x,y) {return( 1/2*pi/3*(x+1))}
sigma    = function(x,y) {return( exp(1/4*y) )}
mu       = function(x,y) {return(x+y)}
tau = function(x,y) {return(0.1)}
funlist = list(mu = mu,sigma = sigma,lambda_1 = lambda_1,lambda_2 = lambda_2,theta = theta, tau = tau)
writeLines(c(""), paste0("Simulations/",sim_id,"/pars.txt"))
sink(paste0("Simulations/",sim_id,"/pars.txt"),append = T)
for (name in names(funlist))
{
  print(name)
  print(funlist[[name]])
  funlist[[name]] = Vectorize(funlist[[name]])
}
sink()

spatial_parameters = create_spatial_parameters(grid, funlist)
spatial_parameters$theta_deg = 180/pi*spatial_parameters$theta
spatial_parameters$phi = sqrt(spatial_parameters$lambda_1/spatial_parameters$lambda_2)
saveRDS(spatial_parameters,paste0("Simulations/",sim_id,"/real_pars.rds"))
#### Create plot with real parameter functions
p = multiple_heatmaps(spatial_parameters)
ggsave(filename = paste0("Simulations/",sim_id,"/plots/real_pars.pdf"),plot = p,width = 30,height = 15,dpi = "retina")

###Perform direct simulation via Cholesky decomposition in Stan (because it's fast)
sim_model = readRDS("Stan_models/simulate.rds")

data = list(N = N, Locations = as.matrix(t(grid)), lambda_1 = spatial_parameters$lambda_1, 
            lambda_2 = spatial_parameters$lambda_2, theta = spatial_parameters$theta, 
            nu = nu, sigma = spatial_parameters$sigma, tau = spatial_parameters$tau)

simu_fit = rstan::sampling(object = sim_model, data=data,
                     warmup=0, iter=10, chains=1, seed=494838,
                     algorithm="Fixed_param")

realizations = matrix( data = unlist(simu_fit@sim$samples), ncol = 10,byrow = TRUE)

set.seed(14061997)
selected_realization = sample(1:10,1)
for ( i in 1:10) for(j in 1:N) realizations[j,i] = realizations[j,i] + spatial_parameters$mu[j];

simulated_process = data.frame(x_1 = grid$x_1, x_2 = grid$x_2, process = realizations[1:N,selected_realization])
save(simulated_process, file =  paste0("Simulations/",sim_id,"/simulated_process.Rdata"))

simulated_process_list = list()
for ( i in 1:10) simulated_process_list[[i]] = data.frame(x_1 = grid$x_1, x_2 = grid$x_2, process = realizations[1:N,i])

save(simulated_process_list, file =  paste0("Simulations/",sim_id,"/simulated_process_list.Rdata"))


p = ggplot(simulated_process) + geom_tile(mapping = aes(x_1,x_2, fill = process)) + coord_fixed() + theme_pubclean(base_size = 30) + 
  scale_fill_viridis(option = "inferno") + theme(legend.text=element_text(size=15), legend.key.size = unit(1.5, 'cm'))

ggsave(filename = paste0("Simulations/",sim_id,"/plots/realization.pdf"),plot = p,width = 20,height = 15,dpi = "retina")

