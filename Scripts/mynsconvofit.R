mynsconvofit = function (sp.SPDF = NULL, coords = NULL, data = NULL, cov.model = "exponential", 
          mean.model = data ~ 1, mc.locations = NULL, N.mc = NULL, 
          lambda.w = NULL, fixed.nugg2.var = NULL, mean.model.df = NULL, 
          mc.kernels = NULL, fit.radius = NULL, ns.nugget = FALSE, 
          ns.variance = FALSE, ns.mean = FALSE, local.aniso = TRUE, 
          fix.tausq = FALSE, tausq = 0, fix.kappa = FALSE, kappa = 0.5, 
          method = "reml", print.progress = TRUE, local.pars.LB = NULL, 
          local.pars.UB = NULL, global.pars.LB = NULL, global.pars.UB = NULL, 
          local.ini.pars = NULL, global.ini.pars = NULL) 
{
  if (is.null(fit.radius)) {
    cat("\nPlease specify a fitting radius.\n")
  }
  if (is.null(coords)) 
    stop("Please provide a Nx2 matrix of spatial coordinates.")
  if (is.null(data)) 
    stop("Please provide a vector of observed data values.")
  if (is.null(sp.SPDF) == FALSE) {
    if (class(sp.SPDF) != "SpatialPointsDataFrame") {
      stop("Please use a SpatialPointsDataFrame object for the 'sp.SPDF = ' input.")
    }
    coords <- sp.SPDF$coords
    data <- sp.SPDF$data
  }
  coords <- as.matrix(coords)
  N <- dim(coords)[1]
  data <- as.matrix(data, nrow = N)
  p <- dim(data)[2]
  if (cov.model != "matern" & cov.model != "gaussian" & cov.model != 
      "exponential") {
    stop("Please specify a valid covariance model (matern, gaussian, or exponential).")
  }
  if (ns.mean == TRUE) {
    if (ns.nugget == FALSE || ns.variance == FALSE) {
      stop("Cannot use ns.mean = TRUE and either ns.nugget = FALSE or ns.variance = FALSE (currently unsupported).")
    }
  }
  if (fix.tausq == TRUE) {
    if (ns.nugget == FALSE || ns.variance == FALSE) {
      stop("Cannot use fix.tausq == TRUE and either ns.nugget = FALSE or ns.variance = FALSE (currently unsupported).")
    }
  }
  if (is.null(mc.locations) == TRUE) {
    if (is.null(N.mc) == TRUE) {
      cat("Please enter the desired number of mixture component locations. \n")
    }
    lon_min <- min(coords[, 1])
    lon_max <- max(coords[, 1])
    lat_min <- min(coords[, 2])
    lat_max <- max(coords[, 2])
    mc_x <- seq(from = lon_min + 0.5 * (lon_max - lon_min)/floor(sqrt(N.mc)), 
                to = lon_max - 0.5 * (lon_max - lon_min)/floor(sqrt(N.mc)), 
                length = floor(sqrt(N.mc)))
    mc_y <- seq(from = lat_min + 0.5 * (lat_max - lat_min)/floor(sqrt(N.mc)), 
                to = lat_max - 0.5 * (lat_max - lat_min)/floor(sqrt(N.mc)), 
                length = floor(sqrt(N.mc)))
    mc.locations <- expand.grid(mc_x, mc_y)
    mc.locations <- matrix(c(mc.locations[, 1], mc.locations[, 
                                                             2]), ncol = 2, byrow = F)
  }
  K <- dim(mc.locations)[1]
  check.mc.locs <- mc_N(coords, mc.locations, fit.radius)
  cat("\n-----------------------------------------------------------\n")
  cat(paste("Fitting the nonstationary model: ", K, " local models with\nlocal sample sizes ranging between ", 
            min(check.mc.locs), " and ", max(check.mc.locs), ".", 
            sep = ""))
  if (ns.nugget == FALSE & ns.variance == FALSE) {
    cat("\nConstant nugget and constant variance.")
  }
  if (ns.nugget == FALSE & ns.variance == TRUE) {
    cat("\nConstant nugget and spatially-varying variance.")
  }
  if (ns.nugget == TRUE & ns.variance == FALSE) {
    cat("\nSpatially-varying nugget and constant variance.")
  }
  if (ns.nugget == TRUE & ns.variance == TRUE) {
    cat("\nSpatially-varying nugget and spatially-varying variance.")
  }
  if (local.aniso) {
    cat(paste("\nLocally anisotropic ", cov.model, " covariance.", 
              sep = ""))
  }
  if (!local.aniso) {
    cat(paste("\nLocally isotropic ", cov.model, " covariance.", 
              sep = ""))
  }
  if (min(check.mc.locs) < 5) {
    cat("\nWARNING: at least one of the mc locations has too few data points.\n")
  }
  cat("\n-----------------------------------------------------------\n")
  if (is.null(lambda.w) == TRUE) {
    lambda.w <- (0.5 * min(dist(mc.locations)))^2
  }
  if (is.null(fixed.nugg2.var) == TRUE) {
    fixed.nugg2.var <- matrix(0, N, N)
  }
  else {
    if (!is.matrix(fixed.nugg2.var)) {
      fixed.nugg2.var <- diag(fixed.nugg2.var)
    }
  }
  if (is.null(mean.model.df) == TRUE) {
    OLS.model <- lm(mean.model, x = TRUE)
    Xmat <- matrix(unname(OLS.model$x), nrow = N)
    beta.names <- colnames(OLS.model$x)
  }
  if (is.null(mean.model.df) == FALSE) {
    OLS.model <- lm(mean.model, x = TRUE, data = mean.model.df)
    Xmat <- matrix(unname(OLS.model$x), nrow = N)
    beta.names <- colnames(OLS.model$x)
  }
  q <- ncol(Xmat)
  lon_min <- min(coords[, 1])
  lon_max <- max(coords[, 1])
  lat_min <- min(coords[, 2])
  lat_max <- max(coords[, 2])
  max.distance <- sqrt(sum((c(lon_min, lat_min) - c(lon_max, 
                                                    lat_max))^2))
  if (p > 1) {
    ols.sigma <- NULL
    for (i in 1:length(names(summary(OLS.model)))) {
      ols.sigma <- c(ols.sigma, summary(OLS.model)[[i]]$sigma)
    }
    resid.var <- (max(ols.sigma))^2
  }
  if (p == 1) {
    resid.var <- summary(OLS.model)$sigma^2
  }
  if (is.null(local.pars.LB) == TRUE) {
    lam1.LB <- 1e-05
    lam2.LB <- 1e-05
    tausq.local.LB <- 1e-05
    sigmasq.local.LB <- 1e-05
    kappa.local.LB <- 1e-05
  }
  if (is.null(local.pars.LB) == FALSE) {
    lam1.LB <- local.pars.LB[1]
    lam2.LB <- local.pars.LB[2]
    tausq.local.LB <- local.pars.LB[3]
    sigmasq.local.LB <- local.pars.LB[4]
    kappa.local.LB <- local.pars.LB[5]
  }
  if (is.null(global.pars.LB) == TRUE) {
    tausq.global.LB <- 1e-05
    sigmasq.global.LB <- 1e-05
    kappa.global.LB <- 1e-05
  }
  if (is.null(global.pars.LB) == FALSE) {
    tausq.global.LB <- global.pars.LB[1]
    sigmasq.global.LB <- global.pars.LB[2]
    kappa.global.LB <- global.pars.LB[3]
  }
  if (is.null(local.pars.UB) == TRUE) {
    lam1.UB <- max.distance/4
    lam2.UB <- max.distance/4
    tausq.local.UB <- 4 * resid.var
    sigmasq.local.UB <- 4 * resid.var
    kappa.local.UB <- 30
  }
  if (is.null(local.pars.UB) == FALSE) {
    lam1.UB <- local.pars.UB[1]
    lam2.UB <- local.pars.UB[2]
    tausq.local.UB <- local.pars.UB[3]
    sigmasq.local.UB <- local.pars.UB[4]
    kappa.local.UB <- local.pars.UB[5]
  }
  if (is.null(global.pars.UB) == TRUE) {
    tausq.global.UB <- 4 * resid.var
    sigmasq.global.UB <- 4 * resid.var
    kappa.global.UB <- 30
  }
  if (is.null(global.pars.UB) == FALSE) {
    tausq.global.UB <- global.pars.UB[1]
    sigmasq.global.UB <- global.pars.UB[2]
    kappa.global.UB <- global.pars.UB[3]
  }
  if (is.null(local.ini.pars) == TRUE) {
    lam1.init <- max.distance/10
    lam2.init <- max.distance/10
    tausq.local.init <- 0.1 * resid.var
    sigmasq.local.init <- 0.9 * resid.var
    kappa.local.init <- 1
  }
  if (is.null(local.ini.pars) == FALSE) {
    lam1.init <- local.ini.pars[1]
    lam2.init <- local.ini.pars[2]
    tausq.local.init <- local.ini.pars[3]
    sigmasq.local.init <- local.ini.pars[4]
    kappa.local.init <- local.ini.pars[5]
  }
  if (is.null(global.ini.pars) == TRUE) {
    tausq.global.init <- 0.1 * resid.var
    sigmasq.global.init <- 0.9 * resid.var
    kappa.global.init <- 1
  }
  if (is.null(global.ini.pars) == FALSE) {
    tausq.global.init <- global.ini.pars[1]
    sigmasq.global.init <- global.ini.pars[2]
    kappa.global.init <- global.ini.pars[3]
  }
  if (is.null(mc.kernels) == TRUE) {
    mc.kernels <- array(NA, dim = c(2, 2, K))
    MLEs.save <- matrix(NA, K, 7)
    beta.GLS.save <- matrix(NA, K, ncol(Xmat))
    beta.cov.save <- array(NA, dim = c(ncol(Xmat), ncol(Xmat), 
                                       K))
    beta.coefs.save <- list()
    coords_mc_dist <- mahalanobis.dist(data.x = coords, 
                                       data.y = mc.locations, vc = diag(2))
    for (k in 1:K) {
      ind_local <- (coords_mc_dist[, k] <= fit.radius)
      temp.locations <- coords[ind_local, ]
      n.fit <- dim(temp.locations)[1]
      temp.data <- as.matrix(data[ind_local, ], nrow = n.fit)
      temp.nugg2.var <- fixed.nugg2.var[ind_local, ind_local]
      if (is.null(mean.model.df) == TRUE) {
        Xtemp <- as.matrix(Xmat[ind_local, ])
      }
      if (is.null(mean.model.df) == FALSE) {
        temp.mmdf <- mean.model.df[ind_local, ]
        Xtemp <- matrix(unname(lm(mean.model, x = TRUE, 
                                  data = temp.mmdf)$x), nrow = n.fit)
      }
      if (print.progress) {
        if (k == 1) {
          cat("Calculating the parameter set for:\n")
        }
        cat("mc location ", k, ", using ", n.fit, " observations...\n", 
            sep = "")
      }
      if (local.aniso) {
        if (cov.model == "matern" || cov.model == "cauchy") {
          if (!fix.kappa) {
            if (!fix.tausq) {
              if (method == "ml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 6 + q), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, lam2.init, 
                                      pi/4, tausq.local.init, sigmasq.local.init, 
                                      kappa.local.init, rep(0, q)), fn = f_loglik, 
                              method = "L-BFGS-B", lower = c(lam1.LB, 
                                                             lam2.LB, 0, tausq.local.LB, sigmasq.local.LB, 
                                                             kappa.local.LB, rep(-Inf, q)), upper = c(lam1.UB, 
                                                                                                      lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB, 
                                                                                                      kappa.local.UB, rep(Inf, q)))
                covMLE <- MLEs$par[1:6]
                betaMLE <- MLEs$par[-c(1:6)]
              }
              if (method == "reml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 6), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, lam2.init, 
                                      pi/4, tausq.local.init, sigmasq.local.init, 
                                      kappa.local.init), fn = f_loglik, 
                              method = "L-BFGS-B", lower = c(lam1.LB, 
                                                             lam2.LB, 0, tausq.local.LB, sigmasq.local.LB, 
                                                             kappa.local.LB), upper = c(lam1.UB, 
                                                                                        lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB, 
                                                                                        kappa.local.UB))
                covMLE <- MLEs$par
              }
            }
            if (fix.tausq) {
              if (method == "ml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 5 + q), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, lam2.init, 
                                      pi/4, sigmasq.local.init, kappa.local.init, 
                                      rep(0, q)), fn = f_loglik, method = "L-BFGS-B", 
                              lower = c(lam1.LB, lam2.LB, 0, sigmasq.local.LB, 
                                        kappa.local.LB, rep(-Inf, q)), upper = c(lam1.UB, 
                                                                                 lam2.UB, pi/2, sigmasq.local.UB, 
                                                                                 kappa.local.UB, rep(Inf, q)))
                covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4:5])
                betaMLE <- MLEs$par[-c(1:5)]
              }
              if (method == "reml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 5), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, lam2.init, 
                                      pi/4, sigmasq.local.init, kappa.local.init), 
                              fn = f_loglik, method = "L-BFGS-B", 
                              lower = c(lam1.LB, lam2.LB, 0, sigmasq.local.LB, 
                                        kappa.local.LB), upper = c(lam1.UB, 
                                                                   lam2.UB, pi/2, sigmasq.local.UB, 
                                                                   kappa.local.UB))
                covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4:5])
              }
            }
          }
          if (fix.kappa) {
            if (!fix.tausq) {
              if (method == "ml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 5 + q), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, lam2.init, 
                                      pi/4, tausq.local.init, sigmasq.local.init, 
                                      rep(0, q)), fn = f_loglik, method = "L-BFGS-B", 
                              lower = c(lam1.LB, lam2.LB, 0, tausq.local.LB, 
                                        sigmasq.local.LB, rep(-Inf, q)), 
                              upper = c(lam1.UB, lam2.UB, pi/2, 
                                        tausq.local.UB, sigmasq.local.UB, 
                                        rep(Inf, q)))
                covMLE <- c(MLEs$par[1:5], kappa)
                betaMLE <- MLEs$par[-c(1:5)]
              }
              if (method == "reml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 5), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, lam2.init, 
                                      pi/4, tausq.local.init, sigmasq.local.init), 
                              fn = f_loglik, method = "L-BFGS-B", 
                              lower = c(lam1.LB, lam2.LB, 0, tausq.local.LB, 
                                        sigmasq.local.LB), upper = c(lam1.UB, 
                                                                     lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB))
                covMLE <- c(MLEs$par[1:5], kappa)
              }
            }
            if (fix.tausq) {
              if (method == "ml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 4 + q), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, lam2.init, 
                                      pi/4, sigmasq.local.init, rep(0, q)), 
                              fn = f_loglik, method = "L-BFGS-B", 
                              lower = c(lam1.LB, lam2.LB, 0, sigmasq.local.LB, 
                                        rep(-Inf, q)), upper = c(lam1.UB, 
                                                                 lam2.UB, pi/2, sigmasq.local.UB, 
                                                                 rep(Inf, q)))
                covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4], 
                            kappa)
                betaMLE <- MLEs$par[-c(1:4)]
              }
              if (method == "reml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 4), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, lam2.init, 
                                      pi/4, sigmasq.local.init), fn = f_loglik, 
                              method = "L-BFGS-B", lower = c(lam1.LB, 
                                                             lam2.LB, 0, sigmasq.local.LB), upper = c(lam1.UB, 
                                                                                                      lam2.UB, pi/2, sigmasq.local.UB))
                covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4], 
                            kappa)
              }
            }
          }
        }
        if (cov.model != "matern" & cov.model != "cauchy") {
          if (!fix.tausq) {
            if (method == "ml") {
              f_loglik <- make_local_lik(locations = temp.locations, 
                                         cov.model = cov.model, data = temp.data, 
                                         Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                         fixed = rep(FALSE, 5 + q), method = method, 
                                         local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                         fix.kappa = fix.kappa, tausq = tausq, 
                                         kappa = kappa)
              MLEs <- optim(par = c(lam1.init, lam2.init, 
                                    pi/4, tausq.local.init, sigmasq.local.init, 
                                    rep(0, q)), fn = f_loglik, method = "L-BFGS-B", 
                            lower = c(lam1.LB, lam2.LB, 0, tausq.local.LB, 
                                      sigmasq.local.LB, rep(-Inf, q)), upper = c(lam1.UB, 
                                                                                 lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB, 
                                                                                 rep(Inf, q)))
              covMLE <- c(MLEs$par[1:5], NA)
              betaMLE <- MLEs$par[-c(1:5)]
            }
            if (method == "reml") {
              f_loglik <- make_local_lik(locations = temp.locations, 
                                         cov.model = cov.model, data = temp.data, 
                                         Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                         fixed = rep(FALSE, 5), method = method, 
                                         local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                         fix.kappa = fix.kappa, tausq = tausq, 
                                         kappa = kappa)
              MLEs <- optim(par = c(lam1.init, lam2.init, 
                                    pi/4, tausq.local.init, sigmasq.local.init), 
                            fn = f_loglik, method = "L-BFGS-B", 
                            lower = c(lam1.LB, lam2.LB, 0, tausq.local.LB, 
                                      sigmasq.local.LB), upper = c(lam1.UB, 
                                                                   lam2.UB, pi/2, tausq.local.UB, sigmasq.local.UB))
              covMLE <- c(MLEs$par[1:5], NA)
            }
          }
          if (fix.tausq) {
            if (method == "ml") {
              f_loglik <- make_local_lik(locations = temp.locations, 
                                         cov.model = cov.model, data = temp.data, 
                                         Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                         fixed = rep(FALSE, 4 + q), method = method, 
                                         local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                         fix.kappa = fix.kappa, tausq = tausq, 
                                         kappa = kappa)
              MLEs <- optim(par = c(lam1.init, lam2.init, 
                                    pi/4, sigmasq.local.init, rep(0, q)), 
                            fn = f_loglik, method = "L-BFGS-B", 
                            lower = c(lam1.LB, lam2.LB, 0, sigmasq.local.LB, 
                                      rep(-Inf, q)), upper = c(lam1.UB, 
                                                               lam2.UB, pi/2, sigmasq.local.UB, rep(Inf, 
                                                                                                    q)))
              covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4], 
                          NA)
              betaMLE <- MLEs$par[-c(1:4)]
            }
            if (method == "reml") {
              f_loglik <- make_local_lik(locations = temp.locations, 
                                         cov.model = cov.model, data = temp.data, 
                                         Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                         fixed = rep(FALSE, 4), method = method, 
                                         local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                         fix.kappa = fix.kappa, tausq = tausq, 
                                         kappa = kappa)
              MLEs <- optim(par = c(lam1.init, lam2.init, 
                                    pi/4, sigmasq.local.init), fn = f_loglik, 
                            method = "L-BFGS-B", lower = c(lam1.LB, 
                                                           lam2.LB, 0, sigmasq.local.LB), upper = c(lam1.UB, 
                                                                                                    lam2.UB, pi/2, sigmasq.local.UB))
              covMLE <- c(MLEs$par[1:3], tausq, MLEs$par[4], 
                          NA)
            }
          }
        }
      }
      if (!local.aniso) {
        if (cov.model == "matern" || cov.model == "cauchy") {
          if (!fix.kappa) {
            if (!fix.tausq) {
              if (method == "ml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 4 + q), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, tausq.local.init, 
                                      sigmasq.local.init, kappa.local.init, 
                                      rep(0, q)), fn = f_loglik, method = "L-BFGS-B", 
                              lower = c(lam1.LB, tausq.local.LB, 
                                        sigmasq.local.LB, kappa.local.LB, 
                                        rep(-Inf, q)), upper = c(lam1.UB, 
                                                                 tausq.local.UB, sigmasq.local.UB, 
                                                                 kappa.local.UB, rep(Inf, q)))
                covMLE <- c(MLEs$par[1], MLEs$par[1], 
                            0, MLEs$par[2:4])
                betaMLE <- MLEs$par[-c(1:4)]
              }
              if (method == "reml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 4), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, tausq.local.init, 
                                      sigmasq.local.init, kappa.local.init), 
                              fn = f_loglik, method = "L-BFGS-B", 
                              lower = c(lam1.LB, tausq.local.LB, 
                                        sigmasq.local.LB, kappa.local.LB), 
                              upper = c(lam1.UB, tausq.local.UB, 
                                        sigmasq.local.UB, kappa.local.UB))
                covMLE <- c(MLEs$par[1], MLEs$par[1], 
                            0, MLEs$par[2:4])
              }
            }
            if (fix.tausq) {
              if (method == "ml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 3 + q), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, sigmasq.local.init, 
                                      kappa.local.init, rep(0, q)), fn = f_loglik, 
                              method = "L-BFGS-B", lower = c(lam1.LB, 
                                                             sigmasq.local.LB, kappa.local.LB, 
                                                             rep(-Inf, q)), upper = c(lam1.UB, 
                                                                                      sigmasq.local.UB, kappa.local.UB, 
                                                                                      rep(Inf, q)))
                covMLE <- c(MLEs$par[1], MLEs$par[1], 
                            0, tausq, MLEs$par[2:3])
                betaMLE <- MLEs$par[-c(1:3)]
              }
              if (method == "reml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 3), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, sigmasq.local.init, 
                                      kappa.local.init), fn = f_loglik, 
                              method = "L-BFGS-B", lower = c(lam1.LB, 
                                                             sigmasq.local.LB, kappa.local.LB), 
                              upper = c(lam1.UB, sigmasq.local.UB, 
                                        kappa.local.UB))
                covMLE <- c(MLEs$par[1], MLEs$par[1], 
                            0, tausq, MLEs$par[2:3])
              }
            }
          }
          if (fix.kappa) {
            if (!fix.tausq) {
              if (method == "ml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 3 + q), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, tausq.local.init, 
                                      sigmasq.local.init, rep(0, q)), fn = f_loglik, 
                              method = "L-BFGS-B", lower = c(lam1.LB, 
                                                             tausq.local.LB, sigmasq.local.LB, 
                                                             rep(-Inf, q)), upper = c(lam1.UB, 
                                                                                      tausq.local.UB, sigmasq.local.UB, 
                                                                                      rep(Inf, q)))
                covMLE <- c(MLEs$par[1], MLEs$par[1], 
                            0, MLEs$par[2:3], kappa)
                betaMLE <- MLEs$par[-c(1:3)]
              }
              if (method == "reml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 3), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, tausq.local.init, 
                                      sigmasq.local.init), fn = f_loglik, 
                              method = "L-BFGS-B", lower = c(lam1.LB, 
                                                             tausq.local.LB, sigmasq.local.LB), 
                              upper = c(lam1.UB, tausq.local.UB, 
                                        sigmasq.local.UB))
                covMLE <- c(MLEs$par[1], MLEs$par[1], 
                            0, MLEs$par[2:3], kappa)
              }
            }
            if (fix.tausq) {
              if (method == "ml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 2 + q), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, sigmasq.local.init, 
                                      rep(0, q)), fn = f_loglik, method = "L-BFGS-B", 
                              lower = c(lam1.LB, sigmasq.local.LB, 
                                        rep(-Inf, q)), upper = c(lam1.UB, 
                                                                 sigmasq.local.UB, rep(Inf, q)))
                covMLE <- c(MLEs$par[1], MLEs$par[1], 
                            0, tausq, MLEs$par[2], kappa)
                betaMLE <- MLEs$par[-c(1:2)]
              }
              if (method == "reml") {
                f_loglik <- make_local_lik(locations = temp.locations, 
                                           cov.model = cov.model, data = temp.data, 
                                           Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                           fixed = rep(FALSE, 2), method = method, 
                                           local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                           fix.kappa = fix.kappa, tausq = tausq, 
                                           kappa = kappa)
                MLEs <- optim(par = c(lam1.init, sigmasq.local.init), 
                              fn = f_loglik, method = "L-BFGS-B", 
                              lower = c(lam1.LB, sigmasq.local.LB), 
                              upper = c(lam1.UB, sigmasq.local.UB))
                covMLE <- c(MLEs$par[1], MLEs$par[1], 
                            0, tausq, MLEs$par[2], kappa)
              }
            }
          }
        }
        if (cov.model != "matern" & cov.model != "cauchy") {
          if (!fix.tausq) {
            if (method == "ml") {
              f_loglik <- make_local_lik(locations = temp.locations, 
                                         cov.model = cov.model, data = temp.data, 
                                         Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                         fixed = rep(FALSE, 3 + q), method = method, 
                                         local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                         fix.kappa = fix.kappa, tausq = tausq, 
                                         kappa = kappa)
              MLEs <- optim(par = c(lam1.init, tausq.local.init, 
                                    sigmasq.local.init, rep(0, q)), fn = f_loglik, 
                            method = "L-BFGS-B", lower = c(lam1.LB, 
                                                           tausq.local.LB, sigmasq.local.LB, 
                                                           rep(-Inf, q)), upper = c(lam1.UB, 
                                                                                    tausq.local.UB, sigmasq.local.UB, 
                                                                                    rep(Inf, q)))
              covMLE <- c(MLEs$par[1], MLEs$par[1], 
                          0, MLEs$par[2:3], NA)
              betaMLE <- MLEs$par[-c(1:3)]
            }
            if (method == "reml") {
              f_loglik <- make_local_lik(locations = temp.locations, 
                                         cov.model = cov.model, data = temp.data, 
                                         Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                         fixed = rep(FALSE, 3), method = method, 
                                         local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                         fix.kappa = fix.kappa, tausq = tausq, 
                                         kappa = kappa)
              MLEs <- optim(par = c(lam1.init, tausq.local.init, 
                                    sigmasq.local.init), fn = f_loglik, 
                            method = "L-BFGS-B", lower = c(lam1.LB, 
                                                           tausq.local.LB, sigmasq.local.LB), 
                            upper = c(lam1.UB, tausq.local.UB, sigmasq.local.UB))
              covMLE <- c(MLEs$par[1], MLEs$par[1], 
                          0, MLEs$par[2:3], NA)
            }
          }
          if (fix.tausq) {
            if (method == "ml") {
              f_loglik <- make_local_lik(locations = temp.locations, 
                                         cov.model = cov.model, data = temp.data, 
                                         Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                         fixed = rep(FALSE, 2 + q), method = method, 
                                         local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                         fix.kappa = fix.kappa, tausq = tausq, 
                                         kappa = kappa)
              MLEs <- optim(par = c(lam1.init, sigmasq.local.init, 
                                    rep(0, q)), fn = f_loglik, method = "L-BFGS-B", 
                            lower = c(lam1.LB, sigmasq.local.LB, 
                                      rep(-Inf, q)), upper = c(lam1.UB, 
                                                               sigmasq.local.UB, rep(Inf, q)))
              covMLE <- c(MLEs$par[1], MLEs$par[1], 
                          0, tausq, MLEs$par[2], NA)
              betaMLE <- MLEs$par[-c(1:2)]
            }
            if (method == "reml") {
              f_loglik <- make_local_lik(locations = temp.locations, 
                                         cov.model = cov.model, data = temp.data, 
                                         Xmat = Xtemp, nugg2.var = temp.nugg2.var, 
                                         fixed = rep(FALSE, 2), method = method, 
                                         local.aniso = local.aniso, fix.tausq = fix.tausq, 
                                         fix.kappa = fix.kappa, tausq = tausq, 
                                         kappa = kappa)
              MLEs <- optim(par = c(lam1.init, sigmasq.local.init), 
                            fn = f_loglik, method = "L-BFGS-B", 
                            lower = c(lam1.LB, sigmasq.local.LB), 
                            upper = c(lam1.UB, sigmasq.local.UB))
              covMLE <- c(MLEs$par[1], MLEs$par[1], 
                          0, tausq, MLEs$par[2], NA)
            }
          }
        }
      }
      if (MLEs$convergence != 0) {
        if (MLEs$convergence == 52) {
          cat(paste("  There was a NON-FATAL error with optim(): \n  ", 
                    MLEs$convergence, "  ", MLEs$message, "\n", 
                    sep = ""))
        }
        else {
          cat(paste("  There was an error with optim(): \n  ", 
                    MLEs$convergence, "  ", MLEs$message, "\n", 
                    sep = ""))
        }
      }
      mc.kernels[, , k] <- kernel_cov(covMLE[1:3])
      if (ns.mean) {
        if (method == "ml") {
          beta.GLS.save[k, ] <- betaMLE
          beta.cov.save[, , k] <- matrix(NA, q, q)
          beta.coefs.save[[k]] <- NA
        }
        if (method == "reml") {
          dist.k <- mahalanobis.dist(data.x = temp.locations, 
                                     vc = mc.kernels[, , k])
          NS.cov.k <- covMLE[5] * cov_spatial(dist.k, 
                                              cov.model = cov.model, cov.pars = c(1, 1), 
                                              kappa = covMLE[6])
          Data.cov.k <- NS.cov.k + diag(rep(covMLE[4], 
                                            n.fit)) + temp.nugg2.var
          Data.chol.k <- chol(Data.cov.k)
          tX.Cinv.k <- t(backsolve(Data.chol.k, backsolve(Data.chol.k, 
                                                          Xtemp, transpose = TRUE)))
          beta.cov.k <- chol2inv(chol(tX.Cinv.k %*% 
                                        Xtemp))/p
          beta.GLS.k <- (p * beta.cov.k %*% tX.Cinv.k %*% 
                           temp.data)/p
          beta.GLS.save[k, ] <- beta.GLS.k
          beta.cov.save[, , k] <- beta.cov.k
          beta.coefs.save[[k]] <- data.frame(Estimate = beta.GLS.k, 
                                             Std.Error = sqrt(diag(beta.cov.k)), t.val = beta.GLS.k/sqrt(diag(beta.cov.k)))
        }
      }
      MLEs.save[k, ] <- c(n.fit, covMLE)
    }
    MLEs.save <- data.frame(MLEs.save)
    names(MLEs.save) <- c("n", "lam1", "lam2", "eta", "tausq", 
                          "sigmasq", "kappa")
  }
  weights.unnorm <- exp(-mahalanobis.dist(data.x = coords, 
                                          data.y = mc.locations, vc = diag(2))^2/(2 * lambda.w))
  weights <- t(apply(X = weights.unnorm, MARGIN = 1, FUN = function(x) {
    x/sum(x)
  }))
  obs.kernel11 <- rowSums(t(t(weights) * mc.kernels[1, 1, 
  ]))
  obs.kernel22 <- rowSums(t(t(weights) * mc.kernels[2, 2, 
  ]))
  obs.kernel12 <- rowSums(t(t(weights) * mc.kernels[1, 2, 
  ]))
  kernel.ellipses <- array(0, dim = c(2, 2, N))
  kernel.ellipses[1, 1, ] <- obs.kernel11
  kernel.ellipses[1, 2, ] <- obs.kernel12
  kernel.ellipses[2, 1, ] <- obs.kernel12
  kernel.ellipses[2, 2, ] <- obs.kernel22
  if (ns.nugget == TRUE) {
    obs.nuggets <- rowSums(t(t(weights) * MLEs.save$tausq))
  }
  if (ns.variance == TRUE) {
    obs.variance <- rowSums(t(t(weights) * MLEs.save$sigmasq))
  }
  if (ns.mean == TRUE) {
    obs.beta <- NULL
    for (t in 1:ncol(beta.GLS.save)) {
      obs.beta <- cbind(obs.beta, rowSums(t(t(weights) * 
                                              beta.GLS.save[, t])))
    }
  }
  cat("-----------------------------------------------------------\n")
  arg11 <- obs.kernel11
  arg22 <- obs.kernel22
  arg12 <- obs.kernel12
  det1 <- arg11 * arg22 - arg12^2
  mat11_1 <- matrix(arg11, nrow = N) %x% matrix(1, ncol = N)
  mat11_2 <- matrix(1, nrow = N) %x% matrix(arg11, ncol = N)
  mat22_1 <- matrix(arg22, nrow = N) %x% matrix(1, ncol = N)
  mat22_2 <- matrix(1, nrow = N) %x% matrix(arg22, ncol = N)
  mat12_1 <- matrix(arg12, nrow = N) %x% matrix(1, ncol = N)
  mat12_2 <- matrix(1, nrow = N) %x% matrix(arg12, ncol = N)
  mat11 <- 0.5 * (mat11_1 + mat11_2)
  mat22 <- 0.5 * (mat22_1 + mat22_2)
  mat12 <- 0.5 * (mat12_1 + mat12_2)
  det12 <- mat11 * mat22 - mat12^2
  Scale.mat <- diag(sqrt(sqrt(det1))) %*% sqrt(1/det12) %*% 
    diag(sqrt(sqrt(det1)))
  inv11 <- mat22/det12
  inv22 <- mat11/det12
  inv12 <- -mat12/det12
  dists1 <- as.matrix(dist(coords[, 1], upper = T, diag = T))
  dists2 <- as.matrix(dist(coords[, 2], upper = T, diag = T))
  temp1_1 <- matrix(coords[, 1], nrow = N) %x% matrix(1, ncol = N)
  temp1_2 <- matrix(1, nrow = N) %x% matrix(coords[, 1], ncol = N)
  temp2_1 <- matrix(coords[, 2], nrow = N) %x% matrix(1, ncol = N)
  temp2_2 <- matrix(1, nrow = N) %x% matrix(coords[, 2], ncol = N)
  sgn.mat1 <- (temp1_1 - temp1_2 >= 0)
  sgn.mat1[sgn.mat1 == FALSE] <- -1
  sgn.mat2 <- (temp2_1 - temp2_2 >= 0)
  sgn.mat2[sgn.mat2 == FALSE] <- -1
  dists1.sq <- dists1^2
  dists2.sq <- dists2^2
  dists12 <- sgn.mat1 * dists1 * sgn.mat2 * dists2
  Dist.mat <- sqrt(inv11 * dists1.sq + 2 * inv12 * dists12 + 
                     inv22 * dists2.sq)
  if (cov.model != "matern" & cov.model != "cauchy") {
    KAPPA <- NULL
    Unscl.corr <- cov_spatial(Dist.mat, cov.model = cov.model, 
                              cov.pars = c(1, 1), kappa = KAPPA)
    NS.corr <- Scale.mat * Unscl.corr
    if (ns.nugget == FALSE & ns.variance == FALSE) {
      if (print.progress) {
        cat("Calculating the variance parameter MLEs. \n")
      }
      overall.lik1 <- make_global_loglik1(data = data, 
                                          Xmat = Xmat, Corr = NS.corr, nugg2.var = fixed.nugg2.var)
      overall.MLEs <- optim(c(tausq.global.init, sigmasq.global.init), 
                            overall.lik1, method = "L-BFGS-B", lower = c(tausq.global.LB, 
                                                                         sigmasq.global.LB), upper = c(tausq.global.UB, 
                                                                                                       sigmasq.global.UB))
      if (print.progress) {
        if (overall.MLEs$convergence != 0) {
          if (overall.MLEs$convergence == 52) {
            cat(paste("  There was a NON-FATAL error with optim(): \n  ", 
                      overall.MLEs$convergence, "  ", overall.MLEs$message, 
                      "\n", sep = ""))
          }
          else {
            cat(paste("  There was an error with optim(): \n  ", 
                      overall.MLEs$convergence, "  ", overall.MLEs$message, 
                      "\n", sep = ""))
          }
        }
      }
      tausq.MLE <- overall.MLEs$par[1]
      sigmasq.MLE <- overall.MLEs$par[2]
      global.lik <- overall.MLEs$value
      ObsNuggMat <- diag(rep(tausq.MLE, N)) + fixed.nugg2.var
      ObsCov <- sigmasq.MLE * NS.corr
      obs.variance <- rep(sigmasq.MLE, N)
    }
    if (ns.nugget == TRUE & ns.variance == FALSE) {
      if (print.progress) {
        cat("Calculating the variance parameter MLEs. \n")
      }
      overall.lik2 <- make_global_loglik2(data = data, 
                                          Xmat = Xmat, Corr = NS.corr, obs.nuggets = obs.nuggets, 
                                          nugg2.var = fixed.nugg2.var)
      overall.MLEs <- optim(sigmasq.global.init, overall.lik2, 
                            method = "L-BFGS-B", lower = c(sigmasq.global.LB), 
                            upper = c(sigmasq.global.UB))
      if (print.progress) {
        if (overall.MLEs$convergence != 0) {
          if (overall.MLEs$convergence == 52) {
            cat(paste("  There was a NON-FATAL error with optim(): \n  ", 
                      overall.MLEs$convergence, "  ", overall.MLEs$message, 
                      "\n", sep = ""))
          }
          else {
            cat(paste("  There was an error with optim(): \n  ", 
                      overall.MLEs$convergence, "  ", overall.MLEs$message, 
                      "\n", sep = ""))
          }
        }
      }
      sigmasq.MLE <- overall.MLEs$par[1]
      global.lik <- overall.MLEs$value
      ObsNuggMat <- diag(obs.nuggets) + fixed.nugg2.var
      ObsCov <- sigmasq.MLE * NS.corr
      obs.variance <- rep(sigmasq.MLE, N)
    }
    if (ns.nugget == FALSE & ns.variance == TRUE) {
      if (print.progress) {
        cat("Calculating the variance parameter MLEs. \n")
      }
      overall.lik3 <- make_global_loglik3(data = data, 
                                          Xmat = Xmat, Corr = NS.corr, obs.variance = obs.variance, 
                                          nugg2.var = fixed.nugg2.var)
      overall.MLEs <- optim(tausq.global.init, overall.lik3, 
                            method = "L-BFGS-B", lower = c(tausq.global.LB), 
                            upper = c(tausq.global.UB))
      if (print.progress) {
        if (overall.MLEs$convergence != 0) {
          if (overall.MLEs$convergence == 52) {
            cat(paste("  There was a NON-FATAL error with optim(): \n  ", 
                      overall.MLEs$convergence, "  ", overall.MLEs$message, 
                      "\n", sep = ""))
          }
          else {
            cat(paste("  There was an error with optim(): \n  ", 
                      overall.MLEs$convergence, "  ", overall.MLEs$message, 
                      "\n", sep = ""))
          }
        }
      }
      tausq.MLE <- overall.MLEs$par[1]
      global.lik <- overall.MLEs$value
      ObsNuggMat <- diag(rep(tausq.MLE, N)) + fixed.nugg2.var
      ObsCov <- diag(sqrt(obs.variance)) %*% NS.corr %*% 
        diag(sqrt(obs.variance))
    }
    if (ns.nugget == TRUE & ns.variance == TRUE) {
      Cov <- diag(sqrt(obs.variance)) %*% NS.corr %*% 
        diag(sqrt(obs.variance))
      global.lik <- NA
      ObsNuggMat <- diag(obs.nuggets) + fixed.nugg2.var
      ObsCov <- Cov
    }
    kappa.MLE <- NA
  }
  if (cov.model == "matern" || cov.model == "cauchy") {
    if (!fix.kappa) {
      if (ns.nugget == FALSE & ns.variance == FALSE) {
        if (print.progress) {
          cat("Calculating the variance and smoothness parameter MLEs. \n")
        }
        overall.lik1.kappa <- make_global_loglik1_kappa(data = data, 
                                                        Xmat = Xmat, cov.model = cov.model, Scalemat = Scale.mat, 
                                                        Distmat = Dist.mat, nugg2.var = fixed.nugg2.var)
        overall.MLEs <- optim(c(tausq.global.init, sigmasq.global.init, 
                                kappa.global.init), overall.lik1.kappa, method = "L-BFGS-B", 
                              lower = c(tausq.global.LB, sigmasq.global.LB, 
                                        kappa.global.LB), upper = c(tausq.global.UB, 
                                                                    sigmasq.global.UB, kappa.global.UB))
        if (print.progress) {
          if (overall.MLEs$convergence != 0) {
            if (overall.MLEs$convergence == 52) {
              cat(paste("  There was a NON-FATAL error with optim(): \n  ", 
                        overall.MLEs$convergence, "  ", overall.MLEs$message, 
                        "\n", sep = ""))
            }
            else {
              cat(paste("  There was an error with optim(): \n  ", 
                        overall.MLEs$convergence, "  ", overall.MLEs$message, 
                        "\n", sep = ""))
            }
          }
        }
        tausq.MLE <- overall.MLEs$par[1]
        sigmasq.MLE <- overall.MLEs$par[2]
        kappa.MLE <- overall.MLEs$par[3]
        global.lik <- overall.MLEs$value
        ObsNuggMat <- diag(rep(tausq.MLE, N)) + fixed.nugg2.var
        ObsCov <- sigmasq.MLE * Scale.mat * cov_spatial(Dist.mat, 
                                                        cov.model = cov.model, cov.pars = c(1, 1), 
                                                        kappa = kappa.MLE)
        obs.variance <- rep(sigmasq.MLE, N)
      }
      if (ns.nugget == TRUE & ns.variance == FALSE) {
        if (print.progress) {
          cat("Calculating the variance and smoothness parameter MLEs. \n")
        }
        overall.lik2.kappa <- make_global_loglik2_kappa(data = data, 
                                                        Xmat = Xmat, cov.model = cov.model, Scalemat = Scale.mat, 
                                                        Distmat = Dist.mat, obs.nuggets = obs.nuggets, 
                                                        nugg2.var = fixed.nugg2.var)
        overall.MLEs <- optim(c(sigmasq.global.init, 
                                kappa.global.init), overall.lik2.kappa, method = "L-BFGS-B", 
                              lower = c(sigmasq.global.LB, kappa.global.LB), 
                              upper = c(sigmasq.global.UB, kappa.global.UB))
        if (print.progress) {
          if (overall.MLEs$convergence != 0) {
            if (overall.MLEs$convergence == 52) {
              cat(paste("  There was a NON-FATAL error with optim(): \n  ", 
                        overall.MLEs$convergence, "  ", overall.MLEs$message, 
                        "\n", sep = ""))
            }
            else {
              cat(paste("  There was an error with optim(): \n  ", 
                        overall.MLEs$convergence, "  ", overall.MLEs$message, 
                        "\n", sep = ""))
            }
          }
        }
        sigmasq.MLE <- overall.MLEs$par[1]
        kappa.MLE <- overall.MLEs$par[2]
        global.lik <- overall.MLEs$value
        ObsNuggMat <- diag(obs.nuggets) + fixed.nugg2.var
        ObsCov <- sigmasq.MLE * Scale.mat * cov_spatial(Dist.mat, 
                                                        cov.model = cov.model, cov.pars = c(1, 1), 
                                                        kappa = kappa.MLE)
        obs.variance <- rep(sigmasq.MLE, N)
      }
      if (ns.nugget == FALSE & ns.variance == TRUE) {
        if (print.progress) {
          cat("Calculating the variance and smoothness parameter MLEs. \n")
        }
        overall.lik3.kappa <- make_global_loglik3_kappa(data = data, 
                                                        Xmat = Xmat, cov.model = cov.model, Scalemat = Scale.mat, 
                                                        Distmat = Dist.mat, obs.variance = obs.variance, 
                                                        nugg2.var = fixed.nugg2.var)
        overall.MLEs <- optim(c(tausq.global.init, kappa.global.init), 
                              overall.lik3.kappa, method = "L-BFGS-B", lower = c(tausq.global.LB, 
                                                                                 kappa.global.LB), upper = c(tausq.global.UB, 
                                                                                                             kappa.global.UB))
        if (print.progress) {
          if (overall.MLEs$convergence != 0) {
            if (overall.MLEs$convergence == 52) {
              cat(paste("  There was a NON-FATAL error with optim(): \n  ", 
                        overall.MLEs$convergence, "  ", overall.MLEs$message, 
                        "\n", sep = ""))
            }
            else {
              cat(paste("  There was an error with optim(): \n  ", 
                        overall.MLEs$convergence, "  ", overall.MLEs$message, 
                        "\n", sep = ""))
            }
          }
        }
        tausq.MLE <- overall.MLEs$par[1]
        kappa.MLE <- overall.MLEs$par[2]
        global.lik <- overall.MLEs$value
        ObsNuggMat <- diag(rep(tausq.MLE, N)) + fixed.nugg2.var
        ObsCov <- diag(sqrt(obs.variance)) %*% (Scale.mat * 
                                                  cov_spatial(Dist.mat, cov.model = cov.model, 
                                                              cov.pars = c(1, 1), kappa = kappa.MLE)) %*% 
          diag(sqrt(obs.variance))
      }
      if (ns.nugget == TRUE & ns.variance == TRUE) {
        if (print.progress) {
          cat("Calculating the smoothness parameter MLE. \n")
        }
        overall.lik4.kappa <- make_global_loglik4_kappa(data = data, 
                                                        Xmat = Xmat, cov.model = cov.model, Scalemat = Scale.mat, 
                                                        Distmat = Dist.mat, obs.nuggets = obs.nuggets, 
                                                        obs.variance = obs.variance, nugg2.var = fixed.nugg2.var)
        overall.MLEs <- optim(kappa.global.init, overall.lik4.kappa, 
                              method = "L-BFGS-B", lower = c(kappa.global.LB), 
                              upper = c(kappa.global.UB))
        if (print.progress) {
          if (overall.MLEs$convergence != 0) {
            if (overall.MLEs$convergence == 52) {
              cat(paste("  There was a NON-FATAL error with optim(): \n  ", 
                        overall.MLEs$convergence, "  ", overall.MLEs$message, 
                        "\n", sep = ""))
            }
            else {
              cat(paste("  There was an error with optim(): \n  ", 
                        overall.MLEs$convergence, "  ", overall.MLEs$message, 
                        "\n", sep = ""))
            }
          }
        }
        kappa.MLE <- overall.MLEs$par[1]
        global.lik <- overall.MLEs$value
        # ObsNuggMat <- diag(obs.nuggets) + fixed.nugg2.var
        # ObsCov <- diag(sqrt(obs.variance)) %*% (Scale.mat * 
        #   cov_spatial(Dist.mat, cov.model = cov.model, 
        #     cov.pars = c(1, 1), kappa = kappa.MLE)) %*% 
        #   diag(sqrt(obs.variance))
      }
    }
    else {
      kappa.MLE <- kappa
      global.lik <- NA
      #  ObsNuggMat <- diag(obs.nuggets) + fixed.nugg2.var
      #  ObsCov <- diag(sqrt(obs.variance)) %*% (Scale.mat * 
      #    cov_spatial(Dist.mat, cov.model = cov.model, 
      #      cov.pars = c(1, 1), kappa = kappa.MLE)) %*% 
      #    diag(sqrt(obs.variance))
    }
  }
  #  Data.Cov <- ObsNuggMat + ObsCov
  #  Data.Cov.chol <- chol(Data.Cov)
  if (ns.mean == FALSE) {
    tX.Cinv <- t(backsolve(Data.Cov.chol, backsolve(Data.Cov.chol, 
                                                    Xmat, transpose = TRUE)))
    beta.cov <- chol2inv(chol(tX.Cinv %*% Xmat))/p
    beta.GLS <- (p * beta.cov %*% tX.Cinv %*% data)/p
    Mean.coefs <- data.frame(Estimate = beta.GLS, Std.Error = sqrt(diag(beta.cov)), 
                             t.val = beta.GLS/sqrt(diag(beta.cov)))
  }
  if (ns.mean == TRUE) {
    beta.cov <- beta.cov.save
    beta.GLS <- beta.GLS.save
    Mean.coefs <- beta.coefs.save
  }
  if (ns.nugget == TRUE) {
    tausq.out <- obs.nuggets
  }
  if (ns.nugget == FALSE) {
    tausq.out <- tausq.MLE
  }
  if (ns.variance == TRUE) {
    sigmasq.out <- obs.variance
  }
  if (ns.variance == FALSE) {
    sigmasq.out <- sigmasq.MLE
  }
  if (ns.mean == TRUE) {
    beta.out <- obs.beta
    names(beta.out) <- beta.names
  }
  if (ns.mean == FALSE) {
    beta.out <- beta.GLS
  }
  output <- list(mc.kernels = mc.kernels, mc.locations = mc.locations, 
                 MLEs.save = MLEs.save, kernel.ellipses = kernel.ellipses, 
                 data = data, beta.GLS = beta.GLS, beta.cov = beta.cov, 
                 Mean.coefs = Mean.coefs, tausq.est = tausq.out, sigmasq.est = sigmasq.out, 
                 beta.est = beta.out, kappa.MLE = kappa.MLE,# Cov.mat = Data.Cov, Cov.mat.chol = Data.Cov.chol, 
                 cov.model = cov.model, 
                 ns.nugget = ns.nugget, ns.variance = ns.variance, ns.mean = ns.mean, 
                 fixed.nugg2.var = fixed.nugg2.var, coords = coords, 
                 global.loglik = global.lik, Xmat = Xmat, lambda.w = lambda.w, 
                 fix.kappa = fix.kappa, kappa = kappa)
  class(output) <- "NSconvo"
  cat("Done.")
  cat("\n-----------------------------------------------------------\n")
  return(output)
}
