library(tidyverse)
library(purrr)
library(viridis)
library(ggpubr)
library(rstan)
library(foreach)
library(doParallel)
library(radiant.data)
library(RGeostats)
library(plotly)
library(SpatialTools)
options(mc.cores = parallel::detectCores(),Ncpus =  parallel::detectCores())
rstan_options(auto_write = TRUE)


create_spatial_parameters = function(grid, funlist)
{
  temp = grid
  coordnames = names(grid)
  par_names = names(funlist)
  for ( name in par_names)
    temp = cbind(temp, funlist[[name]](grid$x_1,grid$x_2))
  temp = data.frame(temp)
  names(temp) = c(coordnames,par_names)
  return(temp)
}

plot_func_1 = function(df,name)
{
  t = ggplot(df,mapping = aes(x = x_1, y = x_2, fill = value)) +  
    geom_tile()  + coord_fixed() + theme_pubclean(base_size = 30) + scale_fill_viridis(option = "viridis", direction = -1, name = name) + theme(legend.text=element_text(size=30), legend.key.size = unit(1.8, 'cm'), legend.position = "top",axis.title = element_blank())
  return(t)
}

multiple_heatmaps = function(data, plot_func = plot_func_1)
{ 
  data$theta = NULL
  names = names(data)
  par_names = names[3:length(names)]
  pivoted = data %>% pivot_longer( cols = 3:length(names), names_to = "parameter" , values_to = "value")
  nested = pivoted %>% 
    group_by(parameter) %>% 
    nest() %>% 
    mutate(plots = map2(data,parameter,.f = plot_func)) 
  
  lambdas_lim = c(min(c(data$lambda_1,data$lambda_2)) , max(c(data$lambda_1,data$lambda_2)))
  Sigmas_lim = c(min(c(data$Sigma11,data$Sigma22)) , max(c(data$Sigma11,data$Sigma22)))
  
  #nested$plots[[3]] = nested$plots[[3]] + scale_fill_viridis(option = "viridis", direction = -1, name = nested$parameter[3], limits = lambdas_lim)
  #nested$plots[[4]] = nested$plots[[4]] + scale_fill_viridis(option = "viridis", direction = -1, name = nested$parameter[4], limits = lambdas_lim)
  #nested$plots[[7]] = nested$plots[[7]] + scale_fill_viridis(option = "viridis", direction = -1,name = nested$parameter[7], limits = Sigmas_lim)
  #nested$plots[[8]] = nested$plots[[8]] + scale_fill_viridis(option = "viridis", direction = -1, name = nested$parameter[8], limits = Sigmas_lim)
  
  return(gridExtra::grid.arrange(grobs = nested$plots,nrow = 1))
}

create_stan_list = function(covariate, spatial_data, nu, NN_ind, M, fit_mu = 1, fit_sigma = 1, fit_lambda = 1, fit_theta = 1)
{
  P_mu = P_sigma = P_lambda = P_theta = ncol(covariate)
  Y = spatial_data$process
  N = length(Y)
  X_mu = array(covariate,c(1,N,P_mu))
  X_sigma = array(covariate, c(1,N,P_sigma))
  X_lambda = array(covariate, c(1,N,P_lambda))
  X_theta = array(covariate, c(1,N,P_theta))
  uB_mu = array(0,c(1,P_mu))
  uB_sigma = array(0,c(1,P_sigma))
  uB_lambda = array(c(0,0),c(1,P_lambda))
  uB_lambda[1,1] = -5
  uB_theta = array(0,c(1,P_theta))
  VB_mu = array(diag(1,P_mu,P_mu),c(1,P_mu,P_mu))
  VB_sigma = array(diag(1,P_sigma,P_sigma),c(1,P_sigma,P_sigma) )
  VB_lambda = array(diag(1,P_lambda,P_lambda), c(1,P_lambda,P_lambda))
  VB_theta = array(diag(1,P_theta,P_theta), c(1, P_theta, P_theta))
  Locations = t(spatial_data[,1:2])
  
  
  data = list(
    fit_mu = fit_mu,
    fit_sigma = fit_sigma,
    fit_lambda = fit_lambda,
    fit_theta = fit_theta,
    N = N,
    P_mu = P_mu,
    P_sigma = P_sigma,
    P_lambda = P_lambda,
    P_theta = P_theta,
    nu = nu,
    Y = Y,
    X_mu = X_mu,
    X_sigma = X_sigma,
    X_lambda = X_lambda,
    X_theta = X_theta,
    Locations = Locations,
    uB_mu = uB_mu,
    uB_sigma = uB_sigma,
    uB_lambda = uB_lambda,
    uB_theta = uB_theta,
    VB_mu = VB_mu,
    VB_sigma = VB_sigma,
    VB_lambda = VB_lambda,
    VB_theta = VB_theta,
    M = M,
    NN_ind = NN_ind
  )
  return(data)
}


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


from_chains_to_pars_distr = function (chains, fitting_data)
{ 
  aux_list_mat = fitting_data[c("X_mu","X_sigma","X_lambda","X_lambda","X_theta")]
  aux_list_P = fitting_data[c("P_mu","P_sigma","P_lambda","P_lambda","P_theta")]
  names = c("mu","sigma","lambda_1","lambda_2","theta")
  to_ret = list()
  acc = 0
  for( i in 1:5)
  {
    to_ret[[names[i]]] = aux_list_mat[[1]][1,,]%*%t(chains[,(acc + 1):ifelse(i == 1, aux_list_P[[i]] , aux_list_P[[i]] + acc)])
    acc = acc + aux_list_P[[i]]
  }
  to_ret$sigma = exp(to_ret$sigma)
  to_ret$lambda_1 = exp(to_ret$lambda_1)
  to_ret$lambda_2 = exp(to_ret$lambda_2)
  to_ret$theta = pi/2 * rstanarm::invlogit(to_ret$theta)
  to_ret$theta_deg = 180/pi*to_ret$theta
  to_ret$phi = sqrt(to_ret$lambda_1/to_ret$lambda_2)
  
  return(to_ret)
}
  
from_chains_to_pars_est = function(chains,fitting_data,est_fun)
{ 
  names = c("mu","sigma","lambda_1","lambda_2","theta")
  estimates = as.numeric(apply(chains,2,est_fun))
  aux_list_mat = fitting_data[c("X_mu","X_sigma","X_lambda","X_lambda","X_theta")]
  aux_list_P = fitting_data[c("P_mu","P_sigma","P_lambda","P_lambda","P_theta")]
  to_ret = matrix(nrow = fitting_data$N,ncol = 8)
  to_ret = data.frame(to_ret)
  names(to_ret) = c(c("x_1","x_2"),names,"theta_deg")
  to_ret$x_1 = fitting_data$Locations[1,]
  to_ret$x_2 = fitting_data$Locations[2,]
  acc = 0
  for( i in 1:5)
  {
    to_ret[[names[i]]] = aux_list_mat[[1]][1,,]%*%estimates[(acc + 1):ifelse(i == 1, aux_list_P[[i]] , aux_list_P[[i]] + acc)]
    acc = acc + aux_list_P[[i]]
  }
  to_ret$sigma = exp(to_ret$sigma)
  to_ret$lambda_1 = exp(to_ret$lambda_1)
  to_ret$lambda_2 = exp(to_ret$lambda_2)
  to_ret$theta = pi/2 * rstanarm::invlogit(to_ret$theta)
  to_ret$theta_deg = 180/pi*to_ret$theta
  to_ret$phi = sqrt(to_ret$lambda_1/to_ret$lambda_2)
  return(to_ret)
}
 


compute_errors = function (data,real,error_funs, log_rescale = rep(FALSE,6)){
  data$theta = NULL
  par_names = names(data)[-c(1,2)]
  error_to_plot = data
  for ( i in 1:length(par_names))
  { 
    if(!log_rescale[i])
      error_to_plot[[par_names[i]]] = error_funs[[par_names[i]]](data[[par_names[i]]],real[[par_names[i]]])
    else
      error_to_plot[[par_names[i]]] = error_funs[[par_names[i]]](log(data[[par_names[i]]]),log(real[[par_names[i]]]))
    
  }
  return(error_to_plot)
} 

draw_pars_histograms = function (data){
  data$theta = NULL
  order = names(data)
  temp = data %>% pivot_longer(cols = 1:ncol(data),names_to = "parameter") 
  temp$parameter = factor(temp$parameter,levels = order)
  p = ggplot(temp, mapping = aes(x = value, fill = parameter))+ geom_histogram(aes(y = ..density..))+ facet_wrap(facets = "parameter",scale = "free") + theme_pubclean(base_size = 30) + rremove("legend")
  return(p)
} 


relative_err = Vectorize(function(est,real){return((est - real)/real)})

abs_err = Vectorize(function(est,real){return((est - real))})

rpd_err = Vectorize(function(est,real){return((est - real)/( abs(real) + abs(est) ))})
                         
assign_centers = function(data,centers,Ncenters,Ndata){
  distmat = matrix(0,nrow = Ndata, ncol = Ncenters)
  data = as.matrix(data)
  centers = as.matrix(centers)
  for ( i in 1:Ncenters)
    distmat[,i] = (data[,1] - centers[i,1])^2 + (data[,2] - centers[i,2])^2
  return (which.pmin(distmat))
}

extract_spatial_predictions_mean = function(array)
{
  temp = array[,,1]
  B = dim(array)[3]
  
  for ( i in 2:B)
    temp = temp + array[,,i]
  temp = data.frame(round(temp/B, digits = 4))
  names(temp) = c("mu","sigma","lambda_1","lambda_2","theta_deg")
  temp$theta = pi/180 * temp$theta_deg
  temp = temp[,c(1,2,3,4,6,5)]
  return(temp)
}

extract_spatial_predictions_median = function(array,B)
{ 
  N = dim(array)[1]
  npar = dim(array)[2]
  B = dim(array)[3]
  
  temp = matrix(nrow = N,ncol = npar)
  for( j in 1:N)
    for( i in 1:npar)
      temp[j,i] = median (array[j,i,])
    temp = data.frame(temp)
    if (npar == 5)
      names(temp) = c("mu","sigma","lambda_1","lambda_2","theta_deg")
    else
      names(temp) = c("mu","sigma","lambda_1","lambda_2","theta_deg","nugget")
    temp$theta = pi/180 * temp$theta_deg
    temp$phi = sqrt(temp$lambda_1/temp$lambda_2)
    return(temp)
}


RDD_3d_plots = function(estimates, real)
{
  x_1 = unique(estimates[,1])
  x_2 = unique(estimates[,1])
  names = names(estimates)[3:length(names(estimates))]
  plist = list()
  for ( i in 1:length(names))
  { 
    z = matrix(data = real[[names[i]]], nrow = length(x_1),ncol = length(x_2),byrow = T)
    plist[[i]] = plot_ly(x = estimates[,1], y = estimates[,2], z = estimates[[names[i]]],title = names[i], type = "scatter3d", size = 0.1,color = "red",showlegend = F, scene = paste0("scene",i)) %>% add_surface(inherit = F, x = x_1, y = x_2, z = z,showscale=FALSE, scene = paste0("scene",i))  
  }   
  return(plist)
}



RDD_3d_plots_gp = function(parameter_names,dir_path, x_1,x_2, ord, real_pars, scaling_functions)
{
  nplots = length(parameter_names)
  plist = NULL
  for ( i in 1:nplots){
    if ( parameter_names[i] != "theta_deg")
      fit = rstan::extract(readRDS(paste0(dir_path,"/fit_",parameter_names[i],"_rdd_gp.rds")))$latent_functional
    else
      fit = rstan::extract(readRDS(paste0(dir_path,"/fit_","theta","_rdd_gp.rds")))$latent_functional
    fit = colMeans(data.frame(fit))
    fit = scaling_functions[[i]](fit)
    fit_rearrange = fit
    for ( j in 1:length(fit))
      fit_rearrange[ord[j]] = fit[j]
    z_real = matrix(data = real_pars[[parameter_names[i]]], nrow = length(x_1),ncol = length(x_2),byrow = T)
    fit_rearrange = matrix(data = fit_rearrange, nrow = length(x_1),ncol = length(x_2),byrow = T)
    plist[[i]] = plot_ly(x = x_1, y = x_2, z = fit_rearrange,scene = paste0("scene",i),type = "surface",showscale = F) %>% add_surface(x = x_1, y = x_2, z = z_real,scene = paste0("scene",i), color = "red",showscale = F) 
  }
  return(plist)
}

RDD_3d_plots_phi = function(dir_path, x_1,x_2, ord, real_pars)
{
  fit1 =  rstan::extract(readRDS(paste0(dir_path,"/fit_","lambda_1","_rdd_gp.rds")))$latent_functional
  fit2 =  rstan::extract(readRDS(paste0(dir_path,"/fit_","lambda_2","_rdd_gp.rds")))$latent_functional
  fit1 = exp(colMeans(fit1))
  fit2 = exp(colMeans(fit2))
  fit = sqrt(fit1/fit2)
  fit_rearrange = fit
  for ( j in 1:length(fit))
    fit_rearrange[ord[j]] = fit[j]
  z_real = sqrt(real_pars$lambda_1/real_pars$lambda_2)
  z_real = matrix(data = z_real, nrow = length(x_1),ncol = length(x_2),byrow = T)
  fit_rearrange = matrix(data = fit_rearrange, nrow = length(x_1),ncol = length(x_2),byrow = T)
  return(  plot_ly(x = x_1, y = x_2, z = fit_rearrange,scene = paste0("scene",6),type = "surface",showscale = F) %>% add_surface(x = x_1, y = x_2, z = z_real,scene = paste0("scene",6), color = "red",showscale = F) 
)
}

Sigma11 = function(lambda_1,lambda_2, theta){ return(sqrt(lambda_1 )* cos(theta)^2 + sqrt(lambda_2)*sin(theta)^2)}
Sigma22 = function(lambda_1,lambda_2, theta){ return(sqrt(lambda_2) * cos(theta)^2 + sqrt(lambda_1 )*sin(theta)^2)}
Sigma12 = function(lambda_1,lambda_2, theta){ return(sqrt(lambda_1 ) *sin(theta) * cos(theta) - sqrt(lambda_2)*sin(theta)*cos(theta))}

library(plotrix)

library(ggforce)
plot_ellipses = function(data, l_out = 10,scale = 1){
  seq_1 = unique(data$x_1)
  seq_2 = unique(data$x_2)
  seq_1 = seq_1[floor(seq(1,length(seq_1), length.out = l_out))]
  seq_2 = seq_2[floor(seq(1,length(seq_2), length.out = l_out))]
  data_t = data %>% dplyr::filter( x_1 %in% seq_1 , x_2 %in% seq_2)
  p = ggplot(data_t, aes(x0 = x_1, y0 = x_2, a = 3*scale*sqrt((lambda_1)), b = 3*scale*sqrt((lambda_2)), angle = theta_deg*pi/180)) +
     geom_ellipse() + geom_point(aes(x = x_1, y = x_2), size = 3) + theme_pubclean(base_size = 30) + coord_fixed() + ggtitle("Anisotropy Ellipses")+theme(axis.title = element_blank(),plot.title = element_text(margin = margin(60,0,0,0)))
  return(p)

}

add_ellipses = function(data, ggobj, color = "black",scale = 1){
  seq_1 = unique(data$x_1)
  seq_2 = unique(data$x_2)
  seq_1 = seq_1[floor(seq(1,length(seq_1), length.out = 6))]
  seq_2 = seq_2[floor(seq(1,length(seq_2), length.out = 6))]
  data_t = data %>% dplyr::filter( x_1 %in% seq_1 , x_2 %in% seq_2)
  p = ggobj + geom_ellipse(data = data_t, aes(x0 = x_1, y0 = x_2, a = 3*scale*sqrt(lambda_1), b = 3*scale*sqrt(lambda_2), angle = theta_deg*pi/180),inherit.aes = F, color = color, size = 0.01)+ geom_point(data = data_t, aes(x = x_1, y = x_2), size = 3) + coord_fixed() +theme(axis.title = element_blank())
  return(p)
  
}

#library(StanHeaders)
#stanFunction("modified_bessel_second_kind", v = 0.5, x = 10 )

Compute_Aniso = function(lambda_1,lambda_2,theta)
{
  N = length(lambda_1)
  Anis = list()
  for (i in 1:N){
    rot = cbind(c(cos(theta[i]), -sin(theta[i])) , c(sin(theta[i]), cos(theta[i])))
    eig = cbind( c(lambda_1[i],0) , c(0,lambda_2[i]))
    Anis[[i]] = rot%*%eig%*%t(rot);
  }
  return(Anis);
  }

matern_ns_corr = function(Locations, Aniso, dets, nu, sigma)
{ 
  N = ncol(Locations)
  sqrsqrdets = sqrt(sqrt(dets))
  norm_const = 2^(1-nu)/gamma(nu)
  C = matrix( nrow = N, ncol = N)
  for ( i in 1:N){
    C[i,i] = sigma[i]^2;
    for ( j in (i+1):N)
    { 
      if(j > N) break;
      metric = 0.5 * (Aniso[[i]] + Aniso[[j]])
      normdet = sqrt(1/(metric[1,1]*metric[2,2] - square(metric[1,2])));
      Q = sqrt( (Locations[,i] - Locations[,j]) %*% solve(metric, Locations[,i] - Locations[,j]));
        C[i,j] = sigma[i]*sigma[j]*sqrsqrdets[i]*sqrsqrdets[j] * normdet * norm_const * Q^nu * modified_bessel_second_kind(nu,Q);
        C[j,i] = C[i,j];
      }
    }
        return (C);
}

trim = Vectorize( function(x,lwr = -Inf, upr = +Inf)
{ if ( x < lwr) return(lwr)
  if(x > upr) return(upr)
  return(x)})

generate_centers = function (coords,center_grid, K , B, threshold = 10)
{
  N_grid = nrow(center_grid)
  N_obs = nrow(coords)
  center_list = list()
  for ( i in 1:B){
    centers = center_grid[sample(1:N_grid, size = K),]
    assignment = assign_centers(coords,centers,K,N_obs)
  while(min(table(assignment)) < threshold){
    print(sum(table(assignment)))
    centers = center_grid[sample(1:N_grid, size = K),] 
    assignment = assign_centers(coords,centers,K,N_obs)
  }
    center_list[[i]] = centers
  }
  return(center_list)
}

RDD_fit = function (data, coords, center_grid,
                    K, B, clusters = 12, dir_path, struct = c("K-Bessel","Nugget Effect"),
                    dirs = seq(0,150, by = 30), param, lower, upper,threshold)
{ 
  assign_centers = assign_centers
  get_param = get_param
  cl = makePSOCKcluster(clusters)
  registerDoParallel(cl)
  N_obs = length(data)
  N_grid = nrow(center_grid)
  spatial_data = cbind(coords,data)
  center_list = generate_centers(coords,center_grid,K,B,threshold)
  writeLines(c(""), paste0(dir_path,"/log.txt"))
  ret = foreach (i=1:B,.packages = c("RGeostats","radiant.data","seqinr")) %dopar% {
  
    centers = center_list[[i]]
    assignment_obs = assign_centers(coords,centers,K,N_obs)
    
    sink(paste0(dir_path,"/log.txt"),append = T)
    
    if(length(struct) == 1) 
      npar = 5 
    else 
      npar = 6 
      mat_param = t(matrix(NA,K,npar))
    for ( j in 1:K)
    { 
      idxs_obs = which(assignment_obs == j)
      db = db.create(spatial_data[idxs_obs,],flag.grid = F,ndim = 2,autoname = F)
      dirs = dirs
      vario = vario.calc(db,dir = dirs)
      fit = model.auto(vario,struct = struct,
                       param = param,lower = lower,
                       upper = upper,verbose = 0,flag.noreduce = T,flag.goulard = T,maxiter = 10000) 
      mat_param[2:npar,j] = unlist(get_param(fit))
      mat_param[1,j] = mean(spatial_data[idxs_obs,3])#global(db,fit,verbose = F,calcul = "mean")$zest
    }
    print(i)
#    print(mat_param)
    sink()
    t(mat_param)
  }
  stopCluster(cl)
  return(list(estimates =  ret , center_list = center_list) )
}

RDD_predict = function (RDD_output, grid)
{ 
  B = length(RDD_output$center_list)
  K = nrow(RDD_output$center_list[[1]])
  ret = list()
  for ( i in 1:B)
  {
    assign = assign_centers(grid,RDD_output$center_list[[i]],K,nrow(grid))
    ret[[i]] = RDD_output$estimates[[i]][assign,]
  }
  return(ret)
}
  
get_param = function(model, cost = 4.9)
{ 
#  print(model)
  ret = list()
  ret$sigma = sqrt(model@basics[[1]]$sill)
  ret$lambda_1 = 0
  ret$lambda_2 = 0
  ret$theta_deg = 0
  whichmin = which.min(c(model@basics[[1]]$aniso.coeffs[1],model@basics[[1]]$aniso.coeffs[2]))
  ratio = min(c(model@basics[[1]]$aniso.coeffs[1],model@basics[[1]]$aniso.coeffs[2]))
  if(whichmin == 1){
    ret$lambda_1 = (model@basics[[1]]$range * ratio / cost)^2
    ret$lambda_2 = (model@basics[[1]]$range / cost)^2
    
    if(model@basics[[1]]$aniso.rotmat[2,1] < 0)
    {
        t = ret$lambda_1 
        ret$lambda_1 = ret$lambda_2
        ret$lambda_2 = t
    }

    ret$theta_deg = acos(model@basics[[1]]$aniso.rotmat[1,1])/pi*180
    if(ret$theta_deg > 90){
         ret$theta_deg = ret$theta_deg - 90
         t = ret$lambda_1 
         ret$lambda_1 = ret$lambda_2
         ret$lambda_2 = t
     }
  }
  
  if(whichmin == 2){
    ret$lambda_2 = (model@basics[[1]]$range * ratio / cost)^2
    ret$lambda_1 = (model@basics[[1]]$range / cost)^2
    if(model@basics[[1]]$aniso.rotmat[2,1] <0)
    {
      t = ret$lambda_1 
      ret$lambda_1 = ret$lambda_2
      ret$lambda_2 = t
    }
    ret$theta_deg = acos(model@basics[[1]]$aniso.rotmat[1,1])/pi*180
     if(ret$theta_deg > 90){
       ret$theta_deg = ret$theta_deg - 90
       t = ret$lambda_1 
       ret$lambda_1 = ret$lambda_2
       ret$lambda_2 = t
     }
  }
  if(length(model@basics) == 2)
    ret$nugget = model@basics[[2]]$sill[1,1]
  #print(ret)
  return(ret)
}  
  

smooth_NSconvo = function(model, coords)
{
  weight_mat = as.matrix(proxy::dist(scale(coords),scale(model$mc.locations))^2)
  weight_mat = exp(-weight_mat/2*model$lambda.w)
  for ( i in 1:nrow(weight_mat)) weight_mat[i,] = weight_mat[i,]/sum(weight_mat[i,])
  weight_mat = as.matrix(weight_mat)
  knot_pars = model$MLEs.save
  knot_pars$n = NULL
  knot_pars$kappa = NULL
  knot_pars$theta = knot_pars$eta
  knot_pars$eta = NULL
  knot_pars$mu = model$beta.GLS[,1]
  knot_pars = as.matrix(knot_pars[,c(6,4,1,2,3,5)])
  par_est = weight_mat%*%knot_pars
  return(
    data.frame( 
    x_1 = coords[,1],
    x_2 = coords[,2],
    mu = par_est[,1],
    sigma = sqrt(par_est[,2]),
    lambda_1 = log(par_est[,3]),
    lambda_2 = log(par_est[,4]),
    nugget = par_est[,5],
    theta_deg = 180/pi * par_est[,6]
  )
  )
}
  
  
  
  
