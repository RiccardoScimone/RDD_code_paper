rm(list = ls())
graphics.off()
library(RGeostats)
library(cluster)
library(tidyverse)
library(rstan)
library(viridis)
library(foreach)
library(doParallel)
library(ggpubr)
source("NNmatrix.R") 

sq_N = 100
N = sq_N^2
a = -1
b = 1
grid = expand.grid(seq(a,b,length.out = sq_N),seq(a,b,length.out = sq_N))
names(grid) = c("x","y")
sim_model = stan_model("simulazioni.stan")
M = 5
NN.matrix <- NNMatrix(coords = grid, n.neighbors = M, n.omp.threads = 12)

## Extract indexes and reorder data
grid = grid[NN.matrix$ord,]
NN_ind = NN.matrix$NN_ind


#rotation = function(theta){return(cbind(c(cos(theta),sin(theta)),c(-sin(theta),cos(theta))))}
theta = pi/2 * rstanarm::invlogit(grid$x^2 + grid$y^2)
range(theta)
lambda_1 = exp ( (-7 - grid$x)) 
lambda_2 = exp ( (-8 + grid$x))
sigma = exp(1/4*grid$x)
mu = grid$x + grid$y
theta_deg = theta/pi*180
nu = 2
simu_data = list(N = N, Locations = as.matrix(t(grid)), lambda_1 = lambda_1, lambda_2 = lambda_2, theta = theta, nu = nu, sigma = sigma)

simu_fit <- sampling(object = sim_model, data=simu_data,
                     warmup=0, iter=100, chains=1, seed=18011996,
                     algorithm="Fixed_param", refresh=10)

realizations = matrix( data = unlist(simu_fit@sim$samples), ncol = 100,byrow = TRUE)
for ( i in 1:100) for(j in 1:N) realizations[j,i] = realizations[j,i] + mu[j];

reals = data.frame(x = grid$x, y = grid$y,mu = mu,sigma = sigma , lambda_1 = log(lambda_1), lambda_2 = log(lambda_2),theta_deg = theta_deg)
#mu = rep(1,N);
# lambda_1 = rep(0.0015,N);
# lambda_2 = rep(0.0005,N);
# rho1 = sqrt(0.001)
# rho2 = sqrt(0.0005)
# theta = rep(pi/4,N);
# theta_deg = theta/pi*180
# sigma = rep(1.5,N);
grid_toplot = data.frame(x = grid$x, y = grid$y, process = realizations[1:N,20])
x11()
p = ggplot(grid_toplot) + geom_tile(mapping = aes(x,y, fill = process)) + scale_fill_viridis(option = "inferno") + coord_fixed() + theme_pubclean(base_size = 16)
plot(p)

library(radiant.data)
assign = function(data,centers,Ncenters,Ndata){
  distmat = matrix(0,nrow = Ndata, ncol = Ncenters)
  for ( i in 1:Ncenters)
    distmat[,i] = (data[,1] - centers[i,1])^2 + (data[,2] - centers[i,2])^2
  return (which.pmin(distmat))
}

get_param = function(model, cost = 4.9)
{
  ret = list()
  ret$sigma = sqrt(model$sill)
  ret$lambda_1 = 0
  ret$lambda_2 = 0
  if(model$aniso.coeffs[1] < model$aniso.coeffs[2]){
    ret$lambda_2 = (model$range/cost)^2
    ret$lambda_1 = (model$range*min(model$aniso.coeffs[1],model$aniso.coeffs[2])/cost)^2
  }
  else
  {
    ret$lambda_1 = (model$range/cost)^2
    ret$lambda_2 = (model$range*min(model$aniso.coeffs[1],model$aniso.coeffs[2])/cost)^2
  }
  
  ret$theta_deg = acos(model$aniso.rotmat[1,1])/pi*180

  # if(ret$theta_deg > 90)
  # {
  #   ret$theta = ret$theta - pi/2
  #   ret$theta_deg = ret$theta_deg - 90
  #   t = ret$lambda_1
  #   ret$lambda_1 = ret$lambda_2
  #   ret$lambda_2 = t
  # }
  return(ret)
}

cl = makePSOCKcluster(12)
registerDoParallel(cl)
B = 360
K = 32

writeLines(c(""), "log.txt")

ret = foreach (i=1:B,.packages = c("RGeostats","radiant.data")) %dopar% {
  centers = matrix(runif(2*K,min = a,max = b), ncol = 2, nrow = K)
  assignment = assign(grid,centers,K,N)
  sink("log.txt",append = TRUE)
  while(min(table(assignment)) < 40){
    print("assign again")
    centers = matrix(runif(2*K,min = a,max = b), ncol = 2, nrow = K)
    assignment = assign(grid,centers,K,N)
  }
  mat_param = t(matrix(NA,N,5))
  for ( j in 1:K)
  { 
    idxs = which(assignment == j)
    db = db.create(grid_toplot[idxs,],flag.grid = F,ndim = 2,autoname = F)
    dirs = seq(0,160, by = 20)
    vario = vario.calc(db,dir = dirs)
    fit = model.auto(vario,struct = c("K-Bessel"),param = c("R=0.5","P=2","A=45"),lower = c("R=0","P=2","A=30"),upper = c("R=1","P=2","A = 85"),verbose = 0) 
    mat_param[2:5,idxs] = unlist(get_param(fit[1]))
    mat_param[1,idxs] = mean(grid_toplot[idxs,3])
  }
  print(i)
  sink()
  t(mat_param)
}
stopCluster(cl)
save(grid_toplot, file = "spatial_data.Rdata")
save(ret,file = "secondRDD.Rdata")
save(reals, file = "generators.Rdata")

load("secondRDD.Rdata")
load("spatial_data.Rdata")
load("generators.Rdata")
ret_array = simplify2array(ret)


# extract_par_loc = function(loc_id,array)
# {
#   temp = data.frame(t(as.matrix(array[loc,,])))
#   names(temp) = c("mu","sigma","lambda_1","lambda_2","theta_deg")
#   return(temp)
# }
# 
# 
# graphics.off()
# loc = 3000
# to_plot = extract_par_loc(loc,ret_array)
# to_plot = to_plot %>% pivot_longer(cols = 1:5, names_to = "parameter") %>% group_by(parameter) %>% mutate(mean = mean(value),med = median(value), true = eval(parse(text = parameter))[loc])
# 
# x11()
# ggplot(to_plot, mapping = aes(x = value, fill = parameter)) + geom_histogram()+ geom_vline(aes(xintercept = mean), color = "red") +
#   geom_vline(aes(xintercept = med), color = "black")+ geom_vline(aes(xintercept = true), color = "green") + 
#   facet_wrap(facets = "parameter",scales = "free") + theme_pubclean()


extract_spatial_predictions_mean = function(array)
{
  temp = array[,,1]
  B = dim(array)[3]
  for ( i in 2:B)
    temp = temp + array[,,i]
  temp = data.frame(round(temp/B, digits = 4))
  names(temp) = c("mu","sigma","lambda_1","lambda_2","theta_deg")
  return(temp)
}

extract_spatial_predictions_median = function(array)
{
  temp = matrix(nrow = N,ncol = npar)
  N  = dim(array)[1]
  B = dim(array)[3]
  npar = dim(array)[2]
  for( j in 1:N)
    for( i in 1:npar)
      temp[j,i] = median (array[j,i,])
  temp = data.frame(temp)
  names(temp) = c("mu","sigma","lambda_1","lambda_2","theta_deg")
  return(temp)
}


lambda_1 = log(lambda_1)
lambda_2 = log(lambda_2)
preds = extract_spatial_predictions_median(ret_array,B,N,5)
#preds = extract_spatial_predictions_mean(ret_array,B)

preds$lambda_1 = log(preds$lambda_1)
preds$lambda_2 = log(preds$lambda_2)

preds$x = grid$x
preds$y = grid$y
preds = preds %>% pivot_longer(cols = 1:5, names_to = "parameter")

x11()
plot_func = function(df,name)
{
  t = ggplot(df,mapping = aes(x = x, y = y, fill = value)) +  
    geom_tile() + coord_fixed() + theme_pubclean(base_size = 15) + scale_fill_viridis(option = "inferno", name = name) + theme(legend.text=element_text(size=14), legend.key.size = unit(1, 'cm'))
  return(t)
}


nested = preds %>% 
  group_by(parameter) %>% 
  nest() %>% 
  mutate(plots = map2(data,parameter,.f = plot_func)) 

gridExtra::grid.arrange(grobs = nested$plots,nrow = 2)

x11()
nested_reals = reals %>% pivot_longer(cols = 3:7, names_to = "parameter") %>% 
  group_by(parameter) %>% 
  nest() %>% 
  mutate(plots = map2(data,parameter,.f = plot_func)) 

gridExtra::grid.arrange(grobs = nested_reals$plots,nrow = 2)
x11()
plot_func = function(df,name)
{
  t = ggplot(df,mapping = aes(x = x, y = y, fill =  (value - eval(parse(text = name)))/(eval(parse(text = name)))   )) +  
    geom_tile() + coord_fixed() + theme_pubclean(base_size = 15) + scale_fill_viridis(option = "inferno", name = name) + theme(legend.text=element_text(size=12), legend.key.size = unit(1, 'cm'))
  if(name == "mu")
    t = ggplot(df,mapping = aes(x = x, y = y, fill =  (value - eval(parse(text = name))))) +  
      geom_tile() + coord_fixed() + theme_pubclean(base_size = 15) + scale_fill_viridis(option = "inferno", name = name) + theme(legend.text=element_text(size=12), legend.key.size = unit(1, 'cm'))
  return(t)
}

nested = preds %>% 
  group_by(parameter) %>% 
  nest() %>% 
  mutate(plots = map2(data,parameter,.f = plot_func)) 

gridExtra::grid.arrange(grobs = nested$plots,nrow = 2)


