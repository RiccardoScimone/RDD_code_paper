rm(list = ls())
library(BayesNSGP)
library(tidyverse)
library(rstan)
library(rstanarm)
library(shinystan)
library(cmdstanr)

sq_N = 50
a = -1
b = 1
grid = expand.grid(seq(a,b,length.out = sq_N),seq(a,b,length.out = sq_N))
names(grid) = c("x","y")

### Simulate non stationary with STAN ###
### Begin stationary anisotropic
N = sq_N^2
sigma = rep(1,N)
theta = rep(pi/4,N)
lambda_1 = rep(0.05,N)
lambda_2 = rep(0.1,N)
nu = 2 







