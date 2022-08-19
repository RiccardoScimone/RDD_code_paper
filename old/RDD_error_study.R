rm(list = ls())
graphics.off()
library(tidyverse)
library(ggpubr)
library(viridis)
load("secondRDD.Rdata")
load("Bayes_estimates.Rdata")
load("spatial_data.Rdata")
load("generators.Rdata")
attach(reals)
ret_array = simplify2array(ret)
extract_spatial_predictions_median = function(array,B,N, npar)
{
  temp = matrix(nrow = N,ncol = npar)
  for( j in 1:N)
    for( i in 1:npar)
      temp[j,i] = median (array[j,i,])
    temp = data.frame(temp)
    names(temp) = c("mu","sigma","lambda_1","lambda_2","theta_deg")
    return(temp)
}
plot_func_pred = function(df,name)
{
  t = ggplot(df,mapping = aes(x = x, y = y, fill = value)) +  
    geom_tile() + coord_fixed() + theme_pubclean(base_size = 15) + scale_fill_viridis(option = "inferno", name = name) + theme(legend.text=element_text(size=14), legend.key.size = unit(1, 'cm'))
  return(t)
}

B = 360
N = 10000
preds = extract_spatial_predictions_median(ret_array,B,N,5)
preds$lambda_1 = log(preds$lambda_1)
preds$lambda_2 = log(preds$lambda_2)
preds$x = grid_toplot$x
preds$y = grid_toplot$y
regression_from_RDD = preds

errors = data.frame(mu = rep(0,N))
errors$mu = (reals$mu - preds$mu)
errors[,2:5] = (reals[,4:7] - preds[,2:5])/reals[,4:7]
errors$x = preds$x
errors$y = preds$y

preds = preds %>% pivot_longer(cols = 1:5, names_to = "parameter")
errors = errors %>% pivot_longer(cols = 1:5, names_to = "parameter")

### Plot predictions
x11()
nested_predictions_plots = preds %>% 
  group_by(parameter) %>% 
  nest() %>% 
  mutate(plots = map2(data,parameter,.f = plot_func_pred)) 

gridExtra::grid.arrange(grobs = nested_predictions_plots$plots,nrow = 2)

### Errors
x11()
nested_errors_plots = errors %>% 
  group_by(parameter) %>% 
  nest() %>% 
  mutate(plots = map2(data,parameter,.f = plot_func_pred)) 

gridExtra::grid.arrange(grobs = nested_errors_plots$plots,nrow = 2)


### Errors distributions
x11()
hist_errs = ggplot(errors, mapping = aes(x = value, fill = parameter))+ geom_histogram(aes(y = ..density..))+ facet_wrap(facets = "parameter",scales = "free") + theme_pubclean(base_size = 16)
plot(hist_errs)

### Do the same with bayes
load("Bayes_estimates.Rdata")
errors = data.frame(mu = rep(0,N))
errors$mu = (reals$mu - est_fun$mu)
errors[,2:5] = (reals[,4:7] - est_fun[,2:5])/reals[,4:7]
errors$x = est_fun$x
errors$y = est_fun$y
est_fun = est_fun %>% pivot_longer(cols = 1:5, names_to = "parameter")
errors = errors %>% pivot_longer(cols = 1:5, names_to = "parameter")

### Plot predictions
x11()
nested_predictions_plots = est_fun %>% 
  group_by(parameter) %>% 
  nest() %>% 
  mutate(plots = map2(data,parameter,.f = plot_func_pred)) 

gridExtra::grid.arrange(grobs = nested_predictions_plots$plots,nrow = 2)

### Errors
x11()
nested_errors_plots = errors %>% 
  group_by(parameter) %>% 
  nest() %>% 
  mutate(plots = map2(data,parameter,.f = plot_func_pred)) 

gridExtra::grid.arrange(grobs = nested_errors_plots$plots,nrow = 2)


### Errors distributions
x11()
hist_errs = ggplot(errors, mapping = aes(x = value, fill = parameter))+ geom_histogram(aes(y = ..density..))+ facet_wrap(facets = "parameter",scales = "free") + theme_pubclean(base_size = 16)
plot(hist_errs)


### regression on RDD
X = cbind(rep(1,N),regression_from_RDD$x,regression_from_RDD$y,regression_from_RDD$x^2,regression_from_RDD$y^2, regression_from_RDD$x*regression_from_RDD$y)

for ( i in 1:5)
{
  temp = data.frame(y = regression_from_RDD[,i])
  temp = cbind(temp,X)
  model = lm(y~.,temp)
  regression_from_RDD[,i] = predict(model,temp)
}
errors = regression_from_RDD
errors$mu = (reals$mu - regression_from_RDD$mu)
errors[,2:5] = (reals[,4:7] - regression_from_RDD[,2:5])/reals[,4:7]
errors$x = regression_from_RDD$x
errors$y = regression_from_RDD$y
errors = errors %>% pivot_longer(cols = 1:5, names_to = "parameter")
regression_from_RDD =  regression_from_RDD %>% pivot_longer(cols = 1:5, names_to = "parameter")


x11()
nested_predictions_plots = regression_from_RDD %>% 
  group_by(parameter) %>% 
  nest() %>% 
  mutate(plots = map2(data,parameter,.f = plot_func_pred)) 

gridExtra::grid.arrange(grobs = nested_predictions_plots$plots,nrow = 2)

### Errors

x11()
nested_errors_plots = errors %>% 
  group_by(parameter) %>% 
  nest() %>% 
  mutate(plots = map2(data,parameter,.f = plot_func_pred)) 

gridExtra::grid.arrange(grobs = nested_errors_plots$plots,nrow = 2)

### Errors distributions
x11()
hist_errs = ggplot(errors, mapping = aes(x = value, fill = parameter))+ geom_histogram(aes(y = ..density..))+ facet_wrap(facets = "parameter",scales = "free") + theme_pubclean(base_size = 16)
plot(hist_errs)


