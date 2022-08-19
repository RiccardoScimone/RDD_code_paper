##### Tries with BayesNSGP ####
#install.packages("BayesNSGP")
rm(list = ls())
#graphics.off()
library("BayesNSGP")
library(tidyverse)
library(viridis)
set.seed(18011996)
N = 90^2
x = y = seq(-1,1,length.out = sqrt(N))
grid = expand.grid(x,y)
names(grid) = c("x","y")
theta = (atan(grid$y/grid$x))*as.numeric(grid$x*grid$y >= 0) + (pi/2 + atan(grid$y/grid$x)) * as.numeric(grid$x*grid$y < 0)
#theta = rep(pi/4,N)
range(theta)
lambda1 = 0.0005 * as.numeric(grid$x*grid$y >= 0) + 0.003 * as.numeric(grid$x*grid$y < 0)
lambda2 = 0.0005 * as.numeric(grid$x*grid$y < 0) + 0.003 * as.numeric(grid$x*grid$y >= 0)


res_theta = logit(2/pi * theta)
log_1 = log(lambda1)
log_2 = log(lambda2)

Sigma11_vec = inverseEigen(log_1,log_2,res_theta,1)
Sigma22_vec = inverseEigen(log_1,log_2,res_theta,2)
Sigma12_vec = inverseEigen(log_1,log_2,res_theta,3)
sigma = 1
tau = 1
nu = 3
dist_list = nsDist(grid)
Cor_mat <- nsCorr( dist1_sq = dist_list$dist1_sq, dist2_sq = dist_list$dist2_sq,
                   dist12 = dist_list$dist12, Sigma11 = Sigma11_vec,
                   Sigma22 = Sigma22_vec, Sigma12 = Sigma12_vec, nu = nu )
cholesky = chol(Cor_mat)
data = t(cholesky) %*% rnorm(N)
grid$process = data

x11()
p = ggplot(grid) + geom_tile(mapping = aes(x,y, fill = process)) + scale_fill_viridis(option = "inferno")
plot(p)

X_Sigma = cbind(as.numeric(grid$x*grid$y >= 0),as.numeric(grid$x*grid$y < 0),res_theta )
model = nsgpModel(tau_model = "constant",sigma_model = "constant",Sigma_model = "compReg", mu_model = "zero", 
                  likelihood = "NNGP", coords = grid[,1:2], data = grid$process, constants = list(nu = 3, k = 10, 
                  X_Sigma = X_Sigma ))
conf <- configureMCMC(model)
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(model)
Cmcmc <- compileNimble(Rmcmc, project = model)
samples <- runMCMC(Cmcmc, niter = 300, nburnin = 100)

rm(cholesky, Cor_mat, dist_list)

save.image("env10.RData")

N = 500^2
x = y = seq(-1,1,length.out = sqrt(N))
predgrid = expand.grid(x,y)
names(predgrid) = c("x","y")
predtheta = (atan(predgrid$y/predgrid$x))*as.numeric(predgrid$x*predgrid$y >= 0) + (pi/2 + atan(predgrid$y/predgrid$x)) * as.numeric(predgrid$x*predgrid$y < 0)
postpred <- nsgpPredict( model = model, samples = samples, coords.predict = predgrid, PX_Sigma = cbind(as.numeric(predgrid$x*predgrid$y >= 0),as.numeric(predgrid$x*predgrid$y < 0), 
                                                                                                     logit( 2/pi*( predtheta))) )
predgrid$process = postpred$pred[200,]
save.image("env10.RData")
x11()
p = ggplot(predgrid) + geom_tile(mapping = aes(x,y, fill = process)) + scale_fill_viridis(option = "inferno")
plot(p)
############################################################################


rm(list = ls())

set.seed(18011996)
N = 90^2
x = y = seq(-1,1,length.out = sqrt(N))
grid = expand.grid(x,y)
names(grid) = c("x","y")
theta = (atan(grid$y/grid$x))*as.numeric(grid$x*grid$y >= 0) + (pi/2 + atan(grid$y/grid$x)) * as.numeric(grid$x*grid$y < 0)
#theta = rep(pi/4,N)
range(theta)
lambda1 = (0.0005 * as.numeric(grid$x*grid$y >= 0) + 0.003 * as.numeric(grid$x*grid$y < 0)) * (grid$x^2 + grid$y^2 + 1e-6)
lambda2 = (0.0005 * as.numeric(grid$x*grid$y < 0) + 0.003 * as.numeric(grid$x*grid$y >= 0)) * (grid$x^2 + grid$y^2 + 1e-6)


res_theta = logit(2/pi * theta)
log_1 = log(lambda1)
log_2 = log(lambda2)

Sigma11_vec = inverseEigen(log_1,log_2,res_theta,1)
Sigma22_vec = inverseEigen(log_1,log_2,res_theta,2)
Sigma12_vec = inverseEigen(log_1,log_2,res_theta,3)
sigma = 1
tau = 1
nu = 3
dist_list = nsDist(grid)
Cor_mat <- nsCorr( dist1_sq = dist_list$dist1_sq, dist2_sq = dist_list$dist2_sq,
                   dist12 = dist_list$dist12, Sigma11 = Sigma11_vec,
                   Sigma22 = Sigma22_vec, Sigma12 = Sigma12_vec, nu = nu )
cholesky = chol(Cor_mat)
data = t(cholesky) %*% rnorm(N)
grid$process = data

save(grid,file = "generated_data.Rdata")
load("generated_data.Rdata")
x11()
p = ggplot(grid) + geom_tile(mapping = aes(x,y, fill = process)) + scale_fill_viridis(option = "inferno")
plot(p)

#load("generated_data.Rdata")
X_Sigma = cbind(as.numeric(grid$x*grid$y >= 0), as.numeric(grid$x*grid$y < 0) , log(grid$x^2 + grid$y^2 + 1e-6) ,logit( 2/pi*theta ) )


model = nsgpModel(tau_model = "constant", sigma_model = "constant", Sigma_model = "compReg", mu_model = "zero", likelihood = "NNGP", coords = grid[,1:2], data = grid$process, 
                  constants = list(nu = 3, k = 10, p_Sigma = 3, X_Sigma = X_Sigma ) )

conf <- configureMCMC(model)
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(model)
Cmcmc <- compileNimble(Rmcmc, project = model)
samples <- runMCMC(Cmcmc, niter = 200, nburnin = 100)

rm(cholesky, Cor_mat, dist_list)

save.image("env_sim2.RData")

N = 500^2
x = y = seq(-1,1,length.out = sqrt(N))
predgrid = expand.grid(x,y)
names(predgrid) = c("x","y")
predtheta = (atan(predgrid$y/predgrid$x))*as.numeric(predgrid$x*predgrid$y >= 0) + (pi/2 + atan(predgrid$y/predgrid$x)) * as.numeric(predgrid$x*predgrid$y < 0)

PX_Sigma = cbind(as.numeric(predgrid$x*predgrid$y >= 0), as.numeric(predgrid$x*predgrid$y < 0) , log(predgrid$x^2 + predgrid$y^2 + 1e-6) ,logit( 2/pi*predtheta ) )


postpred <- nsgpPredict( model = model, samples = samples, coords.predict = predgrid, PX_Sigma = PX_Sigma )
save.image("env_sim2.RData")
predgrid$process = colMeans(postpred$pred)
x11()
p = ggplot(predgrid) + geom_tile(mapping = aes(x,y, fill = process)) + scale_fill_viridis(option = "inferno")
plot(p)
