rm(list = ls())
graphics.off()
library(BayesNSGP)
library(ggplot2)
library(ggnewscale)
library(viridis)

Colorado_data = read.csv("Colorado/Colorado_data.csv")
Colorado_grid = read.csv("Colorado/Colorado_grid.csv")

n_knot_grid = 8
knot_lon = seq(min(Colorado_data$Longitude), max(Colorado_data$Longitude), length.out = n_knot_grid+1)
knot_lat = seq(min(Colorado_data$Latitude), max(Colorado_data$Latitude), length.out = n_knot_grid+1)

knot_lon = (knot_lon[1:n_knot_grid] + knot_lon[2:(n_knot_grid+1)])/2
knot_lat = (knot_lat[1:n_knot_grid] + knot_lat[2:(n_knot_grid+1)])/2
knot_coord = expand.grid(knot_lon,knot_lat)
names(knot_coord) = c("Longitude", "Latitude")
x11()
ggplot() + geom_point(data = Colorado_data,
                      aes(x = Longitude, y = Latitude, col = logPrecip)) +scale_color_viridis(direction = -1)+ 
  geom_point(data = knot_coord, aes(x = Longitude, y = Latitude),size = 2, shape = 4)

knot_coord = as.matrix(knot_coord)

N = nrow(Colorado_data)
Xmat_mu = cbind(rep(1,N), rnorm(n = N,mean = 0,sd = 1))

constants = list (Sigma_HP1 = c(10,10),Sigma_HP4 = c(10,20), nu = 2,
                  Sigma_knot_coords = knot_coord,Sigma_HP2 = c(5,5),Sigma_HP3 = c(3.85,3.85), 
                  X_mu = Xmat_mu,mu_HP1 = 10, tau_HP1 = 10,maxAnisoRange = 16, k = 29)



Rmodel = nsgpModel(likelihood = "SGV", tau_model =  "constant", sigma_model = "constant",
                   Sigma_model = "npApproxGP", mu_model = "constant", constants = constants, 
                   coords = as.matrix(Colorado_data[,c("Longitude","Latitude")]), data = Colorado_data$logPrecip, nu = 2)
conf = configureMCMC(Rmodel)

conf$printSamplers()
conf$removeSamplers(ind = c(10,11,12))


# Add block samplers
cut_lat = cut(knot_coord[,1],breaks = 4)
cut_lon = cut(knot_coord[,2],breaks = 2)
cut = as.factor(paste0(cut_lat,cut_lon))
table(cut)
x11()
plot(knot_coord[,1],knot_coord[,2], col = cut)

groups = list()
for ( i in 1:8) groups[[i]] = which(as.numeric(cut) == i)
for(g in 1:length(groups)){
  conf$addSampler(target = c(paste0("w1_Sigma[", groups[[g]], "]"),paste0("w2_Sigma[", groups[[g]], "]"),paste0("w3_Sigma[", groups[[g]], "]"))
                  , type = "RW_block", silent = TRUE)
}

conf$printSamplers()

Rmcmc = buildMCMC(conf)
Cmodel = compileNimble(Rmodel) # Compile the model in C++
Cmcmc = compileNimble(Rmcmc, project = Rmodel)

samples = runMCMC(Cmcmc, niter = 100000, nburnin = 50000, thin = 10,nchains = 1,progressBar = T)