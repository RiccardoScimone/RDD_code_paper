rm(list = ls())
graphics.off()
library(RGeostats)
library(cluster)
library(tidyverse)
library(rstan)
library(viridis)
sq_N = 80
N = sq_N^2
a = 1
b = -1
grid = expand.grid(seq(a,b,length.out = sq_N),seq(a,b,length.out = sq_N))
names(grid) = c("x","y")
sim_model = stan_model("simulazioni.stan")



#theta = pi/2* rstanarm::invlogit(grid$x^2 + grid$y^2)
#range(theta)
#ambda_1 = exp(-7.5  + sqrt(grid$x^2 + grid$y^2))
#lambda_2 = exp(-5.7 + 2*sqrt(grid$x^2 + grid$y^2))
#sigma = exp(grid$x^2)
mu = rep(1,N);
nu = 2


lambda_1 = rep(0.001,N);
lambda_2 = rep(0.0005,N);
rho1 = sqrt(0.001)
rho2 = sqrt(0.0005)
theta = rep(pi/3,N);
sigma = rep(1,N);
simu_data = list(N = N, Locations = as.matrix(t(grid)), lambda_1 = lambda_1, lambda_2 = lambda_2, theta = theta, nu = nu, sigma = sigma)

simu_fit <- rstan::sampling(object = sim_model, data=simu_data,
                     warmup=0, iter=100, chains=1, seed=18011996,
                     algorithm="Fixed_param", refresh=10)


realizations = matrix( data = unlist(simu_fit@sim$samples), ncol = 100,byrow = TRUE)
mean(realizations)
for ( i in 1:100) for(j in 1:N) realizations[j,i] = realizations[j,i] + mu[j];

grid_toplot = data.frame(x = grid$x, y = grid$y, process = realizations[1:N,5])
x11()
p = ggplot(grid_toplot) + geom_tile(mapping = aes(x,y, fill = process)) + scale_fill_viridis(option = "inferno")
plot(p)


mean(realizations[1:N,5])


db = db.create(grid_toplot,flag.grid = F,ndim = 2,autoname = F)
print(db)
dirs = seq(0,150, by = 30)

vario = vario.calc(db,dir = dirs)
x11()
plot(vario)
x11()
fit = model.auto(vario,struct = c("K-Bessel"),param = c("V=3","R=0.01","P=2"),lower = c("V=0.1","R=1e-12","P=2"),upper = c("V=100","R=2","P=2"))
print(fit)

bigr = fit[1]
bigr

get_param = function(model, cost = 4.9)
{
  ret = list()
  ret$sigma = sqrt(model$sill)
  ret$lambda_1 = (model$range/cost)^2
  ret$lambda_2 = (model$range*model$aniso.coeffs[2]/cost)^2
  ret$theta = acos(model$aniso.rotmat[1,1])
  ret$theta_deg = ret$theta/pi*180
  return(ret)
}

get_param(bigr)
