plot.NSconvo = function(x, plot.ellipses = TRUE,
         fit.radius = NULL, aniso.mat = NULL, true.mc = NULL,
         ref.loc = NULL, all.pred.locs = NULL, grid = TRUE,
         true.col = 1, aniso.col = 4,
         ns.col = 2, plot.mc.locs = TRUE, ... )
{
  if( !inherits(x, "NSconvo") ){
    warning("Object is not of type NSconvo.")
  }
  else{
    if( plot.ellipses == TRUE ){
      
      mc.locations <- x$mc.locations
      mc.kernels <- x$mc.kernels
      
      K <- dim(mc.locations)[1]
      
      plot(ellipse(mc.kernels[, , 1], centre = c(mc.locations[1,1],mc.locations[1,2]), level = 0.5),
           type = "l", col = ns.col, ... )
      
      if( plot.mc.locs ){ points(mc.locations[1,1], mc.locations[1,2], cex = 1, pch="+" ) }
      
      if (is.null(aniso.mat) == FALSE) {
        lines(ellipse(aniso.mat, centre = mc.locations[1, ],
                      level = 0.5), lty = "dashed", col = aniso.col )
      }
      if (is.null(true.mc) == FALSE) {
        lines(ellipse(true.mc[,,1], centre = mc.locations[1, ],
                      level = 0.5), col = true.col )
      }
      plotrix::draw.circle(mc.locations[1, 1], mc.locations[1, 2], fit.radius,
                           lty = "dashed" )
      for (k in 2:K) {
        
        lines(ellipse(mc.kernels[, , k], centre = c(mc.locations[k,1],mc.locations[k,2]), level = 0.5), col = ns.col )
        
        plotrix::draw.circle(mc.locations[k, 1], mc.locations[k,2], fit.radius, lty = "dashed")
        
        if( plot.mc.locs ){ points(mc.locations[k,1], mc.locations[k,2], cex = 1, pch="+" ) }
        if (is.null(aniso.mat) == FALSE) {
          lines(ellipse(aniso.mat, centre = mc.locations[k,], level = 0.5), lty = "dashed",
                col = aniso.col)
        }
        if (is.null(true.mc) == FALSE) {
          lines(ellipse(true.mc[,,k], centre = mc.locations[k, ],
                        level = 0.5), col = true.col)
        }
      }
    }
    
    if( plot.ellipses == FALSE ){
      
      M <- dim(all.pred.locs)[1]
      
      mc.kern <- x$mc.kernels
      mc.loc <- x$mc.locations
      K <- dim(mc.loc)[1]
      lambda.w <- x$lambda.w
      kappa <- x$kappa.MLE
      cov.model <- x$cov.model
      all.pred.weights <- matrix(NA, M, K)
      for (m in 1:M) {
        for (k in 1:K) {
          all.pred.weights[m, k] <- exp(-sum((all.pred.locs[m,] - mc.loc[k, ])^2)/(2 * lambda.w))
        }
        all.pred.weights[m, ] <- all.pred.weights[m, ]/sum(all.pred.weights[m,])
      }
      pred.weight <- rep(NA, K)
      for (k in 1:K) {
        pred.weight[k] <- exp(-sum((ref.loc - mc.loc[k, ])^2)/(2 * lambda.w))
      }
      pred.weight <- pred.weight/sum(pred.weight)
      pred.kernel.ellipses <- array(0, dim = c(2, 2, M))
      for (m in 1:M) {
        for (k in 1:K) {
          pred.kernel.ellipses[, , m] <- pred.kernel.ellipses[,,m] + all.pred.weights[m, k] * mc.kern[, ,
                                                                                                      k]
        }
      }
      pred.kernel.ellipse <- matrix(rep(0, 4), nrow = 2, ncol = 2)
      for (k in 1:K) {
        pred.kernel.ellipse <- pred.kernel.ellipse + pred.weight[k] *
          mc.kern[, , k]
      }
      Scale.cross <- rep(NA, M)
      Dist.cross <- rep(NA, M)
      Kerneli <- pred.kernel.ellipse
      det_i <- Kerneli[1, 1] * Kerneli[2, 2] - Kerneli[1, 2] *
        Kerneli[2, 1]
      Ui <- chol(Kerneli)
      for (j in 1:M) {
        Kernelj <- pred.kernel.ellipses[, , j]
        det_j <- Kernelj[1, 1] * Kernelj[2, 2] - Kernelj[1,
                                                         2] * Kernelj[2, 1]
        avg_ij <- 0.5 * (Kerneli + Kernelj)
        Uij <- chol(avg_ij)
        det_ij <- avg_ij[1, 1] * avg_ij[2, 2] - avg_ij[1,
                                                       2] * avg_ij[2, 1]
        vec_ij <- backsolve(Uij, (ref.loc - all.pred.locs[j,
        ]), transpose = TRUE)
        Scale.cross[j] <- sqrt(sqrt(det_i * det_j)/det_ij)
        Dist.cross[j] <- sqrt(sum(vec_ij^2))
      }
      Unscl.cross <- geoR::cov.spatial(Dist.cross, cov.model = cov.model,
                                       cov.pars = c(1, 1), kappa = kappa)
      Cov <- matrix(Scale.cross * Unscl.cross, ncol = 1)
      
      if (grid == TRUE) {
        grid_x <- unique(all.pred.locs[, 1])
        grid_y <- unique(all.pred.locs[, 2])
        image.plot(grid_x, grid_y, matrix(Cov, length(grid_x), length(grid_y)),
                   ... )
      }
      if (grid == FALSE) {
        quilt.plot(all.pred.locs, c(Cov), ... )
      }
    }
  }
}
