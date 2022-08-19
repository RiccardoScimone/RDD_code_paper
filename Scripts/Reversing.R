lambda_1 = function(x,y) {return(exp ( (-7 - x)) )}
lambda_2 = function(x,y) {return(exp ( (-8 + x)) )}
theta    = function(x,y) {return( pi/2 * rstanarm::invlogit(x^2 + y^2) )}
sigma    = function(x,y) {return( exp(log(0.3) + 1/4*y) )}
mu       = function(x,y) {return( x+y )}
funlist = list(mu = mu,sigma = sigma,lambda_1 = lambda_1,lambda_2 = lambda_2,theta = theta)
for (name in names(funlist))
{
  funlist[[name]] = Vectorize(funlist[[name]])
}

pars = readRDS("Simulations/simulation_3/real_pars.rds")

m1 = lm(formula = I(rstanarm::logit(2/pi*theta)) ~ I(x_1^2) + I(x_2^2), data = pars)
summary(m1)        
m2 = lm(formula = I(log(sigma)) ~ 1 + x_1 + x_2, data = pars)
summary(m2)
m3 = lm(formula = I(log(lambda_1)) ~ 1 + x_1 + x_2, data = pars)
summary(m3)
m4 = lm(formula = I(log(lambda_2)) ~ 1 + x_1 + x_2, data = pars)
summary(m4)
m5 = lm(formula = mu ~ 1 + x_1 + x_2, data = pars)
summary(m5)

pars = readRDS("Simulations/simulation_1/real_pars.rds")

m1 = lm(formula = I(rstanarm::logit(2/pi*theta)) ~ I(x_1^2) + I(x_2^2), data = pars)
summary(m1)        
m2 = lm(formula = I(log(sigma)) ~ 1 + x_1 + x_2, data = pars)
summary(m2)
m3 = lm(formula = I(log(lambda_1)) ~ 1 + x_1 + x_2, data = pars)
summary(m3)
m4 = lm(formula = I(log(lambda_2)) ~ 1 + x_1 + x_2, data = pars)
summary(m4)
m5 = lm(formula = mu ~ 1 + x_1 + x_2, data = pars)
summary(m5)
