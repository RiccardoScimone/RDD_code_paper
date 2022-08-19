rm(list = ls())
library(RandomFields)
library(tidyverse)
library(viridis)
library(rotations)
set.seed(18011996)
N = 500^2
x = y = seq(-1,1,length.out = sqrt(N))
grid = expand.grid(x,y)
names(grid) = c("x","y")

theta = pi/4
lambda_1 = 0.0006
lambda_2 = 0.0003
Lambda = diag(c(lambda_1,lambda_2))
nu = 3/2
sigma2 = 1

#rotation = function(theta) {return (rbind(c(sin(theta), - cos(theta)), c(cos(theta), sin(theta))))}
#Aniso = rotation(theta)%*%Lambda%*%t(rotation(theta))

model = RMmatern(nu = nu, var = sigma2, Aniso = RMangle(angle = theta, ratio = 2), scale = 0.2)
data = RFsimulate(model = RPgauss(model), x = grid$x, y= grid$y)
grid$process = data@data$variable1

x11()
p = ggplot(grid) + geom_tile(mapping = aes(x,y, fill = process)) + scale_fill_viridis(option = "inferno")
plot(p)