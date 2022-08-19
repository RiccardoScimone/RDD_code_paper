library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

util <- new.env()

par(family="CMU Serif", las=1, bty="l", cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 5))

N = 1101
x <- 22 * (0:(N - 1)) / (N - 1) - 11

alpha_true <- 3
rho_true <- 5.5

simu_data <- list(alpha=alpha_true, rho=rho_true, N=N, x=x)

simu_fit <- stan(file='simu.stan', data=simu_data,
                 warmup=0, iter=4000, chains=1, seed=494838,
                 algorithm="Fixed_param", refresh=4000)
source('gp_utility.R', local=util)
x11()
util$plot_gp_prior_realizations(simu_fit, x, "Realizations")

sigma_true <- 2
simu_data$sigma <- sigma_true

set.seed(2595)
simu_fit <- stan(file='simu_normal.stan', data=simu_data,
                 warmup=0, iter=4000, chains=1, seed=494838,
                 algorithm="Fixed_param", refresh=4000)



f <- extract(simu_fit)$f[1,]
y <- extract(simu_fit)$y[1,]

c_mid_teal="#487575"

plot(x, f, type="l", lwd=2, xlab="x", ylab="y",
     xlim=c(-11, 11), ylim=c(-10, 10))
points(x, y, col="white", pch=16, cex=0.6)
points(x, y, col=c_mid_teal, pch=16, cex=0.4)






