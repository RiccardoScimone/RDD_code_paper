rm(list = ls())
graphics.off()
library(rstan)
library(rstanarm)
library(bayesplot)
library(ggpubr)
source("NNMatrix.R")
fit = read_rds("fit2_nngp_more.rds")
fit_plot = extract(fit)
fit_plot_2 = extract(fit, inc_warmup = F, permute = F)[,,15:28]

color_scheme_set("mix-blue-pink")
p <- mcmc_hist(fit_plot_2) + facet_text(size = 15) + theme_classic(base_size = 22)

x11()
plot(p)


samples_from_chain = function(filename) {
  stanfit = rstan::read_stan_csv(filename)
  samples = data.frame(stanfit@sim$samples)
  samples %>% 
    filter(lp__ != 0) %>%  # Remove non-sampled iterations
    mutate(
      iter = 1:nrow(.),
      chain = gsub('.csv', '', last(strsplit(filename, '_')[[1]]))  # Get chain number [filename_1.csv]
    )
}

chains = rbind(
  samples_from_chain('chain_data_easy_1.csv'),
  samples_from_chain('chain_data_easy_2.csv'),
  samples_from_chain('chain_data_easy_3.csv'),
  samples_from_chain('chain_data_easy_4.csv'),
  samples_from_chain('chain_data_easy_5.csv'),
  samples_from_chain('chain_data_easy_6.csv')
)
to_plot = chains[,c(31:60,62,63)]

to_plot = (to_plot %>% filter(iter > 500, chain != "1", chain != "6"))[,-c(31,32)]

to_plot = as.numeric(apply(to_plot,2,median))
sq_N = 100
N = sq_N^2
a = -1
b = 1
grid = expand.grid(seq(a,b,length.out = sq_N),seq(a,b,length.out = sq_N))
names(grid) = c("x","y")
nu = 2
M = 5
NN.matrix <- NNMatrix(coords = grid, n.neighbors = M, n.omp.threads = 12)

## Extract indexes and reorder data
grid = grid[NN.matrix$ord,]
NN_ind = NN.matrix$NN_ind

X = cbind(rep(1,N),grid$x,grid$y,grid$x^2,grid$y^2, grid$x*grid$y)
X_mu = X
X_sigma = X
X_lambda = X
X_theta = X

est_fun = data.frame(mu = X_mu%*%to_plot[1:6], sigma = exp(X_sigma%*%to_plot[7:12]), lambda_1 = X_lambda%*%to_plot[13:18], lambda_2 = X_lambda%*%to_plot[19:24], theta = pi/2*rstanarm::invlogit(X_theta%*%to_plot[25:30])/pi*180 )
est_fun$x = grid$x
est_fun$y = grid$y
save(est_fun, file = "Bayes_estimates.Rdata")
x11()
plot_func = function(df,name)
{
  t = ggplot(df,mapping = aes(x = x, y = y, fill = value)) +  
    geom_tile() + coord_fixed() + theme_pubclean(base_size = 15) + scale_fill_viridis(option = "inferno", name = name) + theme(legend.text=element_text(size=14), legend.key.size = unit(1, 'cm'))
  return(t)
}
nested_reals = est_fun %>% pivot_longer(cols = 1:5, names_to = "parameter") %>% 
  group_by(parameter) %>% 
  nest() %>% 
  mutate(plots = map2(data,parameter,.f = plot_func)) 

gridExtra::grid.arrange(grobs = nested_reals$plots,nrow = 2)

#rotation = function(theta){return(cbind(c(cos(theta),sin(theta)),c(-sin(theta),cos(theta))))}
theta = pi/2 * rstanarm::invlogit(grid$x^2 + grid$y^2)
range(theta)
lambda_1 = exp ( (-7 - grid$x)) 
lambda_2 = exp ( (-8 + grid$x))
sigma = exp(1/4*grid$x)
mu = grid$x + grid$y
theta_deg = theta/pi*180
lambda_1 = log(lambda_1)
lambda_2 = log(lambda_2)
theta = theta_deg

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

nested = est_fun %>% pivot_longer(cols = 1:5, names_to = "parameter") %>% 
  group_by(parameter) %>% 
  nest() %>% 
  mutate(plots = map2(data,parameter,.f = plot_func)) 

gridExtra::grid.arrange(grobs = nested$plots,nrow = 2)
