#' Computes EI for minimization
#' @title Expected Improvement criterion
#' @param x matrix of new designs, one point per row (size n x d)
#' @param model \code{homGP} or \code{hetGP} model, or their TP equivalents, including inverse matrices
#' @param cst optional plugin value used in the EI, see details
#' @param preds optional predictions at \code{x} to avoid recomputing if already done
#' @note This is a beta version at this point.
#' @details 
#' \code{cst} is classically the observed minimum in the deterministic case. 
#' In the noisy case, the min of the predictive mean works fine.
#'  
#' @references
#' Mockus, J.; Tiesis, V. & Zilinskas, A. (1978).
#' The application of Bayesian methods for seeking the extremum Towards Global Optimization, Amsterdam: Elsevier, 2, 2.\cr \cr
#' 
#' Vazquez E, Villemonteix J, Sidorkiewicz M, Walter E (2008). 
#' Global Optimization Based on Noisy Evaluations: An Empirical Study of Two Statistical Approaches, 
#' Journal of Physics: Conference Series, 135, IOP Publishing.\cr \cr
#' 
#' A. Shah, A. Wilson, Z. Ghahramani (2014), Student-t processes as alternatives to Gaussian processes, Artificial Intelligence and Statistics, 877--885.
#' @importFrom stats dt qt
#' @export
#' @examples 
#' ## Optimization example
#' set.seed(42)
#' 
#' 
#' ## Noise field via standard deviation
#' noiseFun <- function(x, coef = 1.1, scale = 1){
#' if(is.null(nrow(x)))
#'  x <- matrix(x, nrow = 1)
#'    return(scale*(coef + cos(x * 2 * pi)))
#' }
#' 
#' ## Test function defined in [0,1]
#' ftest <- function(x){
#' if(is.null(nrow(x)))
#' x <- matrix(x, ncol = 1)
#' return(f1d(x) + rnorm(nrow(x), mean = 0, sd = noiseFun(x)))
#' }
#' 
#' n_init <- 10 # number of unique designs
#' N_init <- 100 # total number of points
#' X <- seq(0, 1, length.out = n_init)
#' X <- matrix(X[sample(1:n_init, N_init, replace = TRUE)], ncol = 1)
#' Z <- ftest(X)
#' 
#' ## Predictive grid
#' ngrid <- 51
#' xgrid <- seq(0,1, length.out = ngrid)
#' Xgrid <- matrix(xgrid, ncol = 1)
#' 
#' model <- mleHetGP(X = X, Z = Z, lower = 0.001, upper = 1)
#' 
#' EIgrid <- crit_EI(Xgrid, model)
#' preds <- predict(x = Xgrid, model)
#' 
#' par(mar = c(3,3,2,3)+0.1)
#' plot(xgrid, f1d(xgrid), type = 'l', lwd = 1, col = "blue", lty = 3,
#' xlab = '', ylab = '', ylim = c(-8,16))
#' points(X, Z)
#' lines(Xgrid, preds$mean, col = 'red', lwd = 2)
#' lines(Xgrid, qnorm(0.05, preds$mean, sqrt(preds$sd2)), col = 2, lty = 2)
#' lines(Xgrid, qnorm(0.95, preds$mean, sqrt(preds$sd2)), col = 2, lty = 2)
#' lines(Xgrid, qnorm(0.05, preds$mean, sqrt(preds$sd2 + preds$nugs)), col = 3, lty = 2)
#' lines(Xgrid, qnorm(0.95, preds$mean, sqrt(preds$sd2 + preds$nugs)), col = 3, lty = 2)
#' par(new = TRUE)
#' plot(NA, NA, xlim = c(0, 1), ylim = c(0,max(EIgrid)), axes = FALSE, ylab = "", xlab = "")
#' lines(xgrid, EIgrid, lwd = 2, col = 'cyan')
#' axis(side = 4)
#' mtext(side = 4, line = 2, expression(EI(x)), cex = 0.8)
#' mtext(side = 2, line = 2, expression(f(x)), cex = 0.8)
crit_EI <- function(x, model, cst = NULL, preds = NULL){
  if(is.null(cst)) cst <- min(predict(model, x = model$X0)$mean)
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(preds)) preds <- predict(model, x = x)
  
  if(is(model, "homTP") || is(model, "hetTP")){
    gamma <- (cst - preds$mean)/sqrt(preds$sd2)
    res <- (cst - preds$mean) * pt(gamma, df = model$nu + length(model$Z))
    res <- res + sqrt(preds$sd2) * (1 + (gamma^2 - 1)/(model$nu + length(model$Z) - 1)) * dt(x = gamma, df = model$nu + length(model$Z))
    res[which(res < 1e-12)] <- 0 # for stability
    return(res)
  }
  
  xcr <- (cst - preds$mean)/sqrt(preds$sd2)
  res <- (cst - preds$mean) * pnorm(xcr)
  res <- res + sqrt(preds$sd2) * dnorm(xcr)
  res[which(preds$sd2 < sqrt(.Machine$double.eps) | res < 1e-12)] <- 0 # for stability
  return(res) 
}

#' Derivative of EI criterion for GP models
#' @param x matrix for the new design (size 1 x d)
#' @param model \code{homGP} or \code{hetGP} model
#' @param cst optional plugin value used in the EI, see details
#' @param preds optional predictions at \code{x} to avoid recomputing if already done
#' @export
#' @seealso \code{\link[hetGP]{crit_EI}} for the criterion
#' @references 
#' Ginsbourger, D. Multiples metamodeles pour l'approximation et l'optimisation de fonctions numeriques multivariables Ecole Nationale Superieure des Mines de Saint-Etienne, Ecole Nationale Superieure des Mines de Saint-Etienne, 2009 \cr \cr
#' Roustant, O., Ginsbourger, D., DiceKriging, DiceOptim: Two R packages for the analysis of computer experiments by kriging-based metamodeling and optimization, Journal of Statistical Software, 2012
deriv_crit_EI <- function(x, model, cst = NULL, preds = NULL){
  if(is.null(cst)) cst <- min(predict(model, x = model$X)$mean)
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(preds)) preds <- predict(model, x = x)
  
  pred_gr <- predict_gr(model, x)
  
  z <- (cst - preds$mean)/sqrt(preds$sd2)
  
  if(is(model, "homGP") || is(model, "hetGP")){
    res <- pred_gr$sd2 / (2 * sqrt(preds$sd2)) * dnorm(z) - pred_gr$mean * pnorm(z)
  }else{
    # dz = - dm/s - z ds/s = -dm/s - z * ds2/(2s2) 
    dz <- -pred_gr$mean/sqrt(preds$sd2) - z * pred_gr$sd2/(2 * matrix(preds$sd2, nrow = nrow(x), ncol(x)))
    
    # d( (cst - m(x)).pt(z(x))) = 
    p1 <- -pred_gr$mean * pt(z, df = model$nu + length(model$Z)) + (cst - preds$mean) * dz * dt(z, df = model$nu + length(model$Z))
    
    a <- model$nu + length(model$Z) - 1
    # d( s(x) (1 + (z^2-1)/(nu + N -1)) dt(z(x)) (in 2 lines)
    p2 <- (pred_gr$sd2/(2*sqrt(preds$sd2)) * (1 + (z^2 - 1)/a) + 2*sqrt(preds$sd2) * z * dz/a) * dt(z, df = model$nu + length(model$Z))
    p2 <- p2 + sqrt(preds$sd2) * (1 + (z^2 - 1)/a) * dz * dlambda(z, model$nu + length(model$Z))
    res <- p1 + p2
  }
  
  res[which(abs(res) < 1e-12)] <- 0 # for stability with optim
  
  return(res)
}

#' Prediction of the gradient of the predictive quantities
#' @param object GP/TP
#' @param x design location
#' @returns A list with elements:
#' \itemize{
#' \item \code{mean}: the gradient of the predictive mean;
#' \item \code{sd2}: the gradient of the predictive variance (not the variance of the gradient GP);
#' }
#' @seealso \code{\link[hetGP]{predict_derivgp}} for the prediction for the derivative GP itself
#' @noRd
predict_gr <- function(object, x){
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  kvec <-  cov_gen(X1 = object$X0, X2 = x, theta = object$theta, type = object$covtype)
  
  dm <- ds2 <- matrix(NA, nrow(x), ncol(x))
  
  for(i in 1:nrow(x)){
    dkvec <- matrix(NA, nrow(object$X0), ncol(x))
    for(j in 1:ncol(x)) dkvec[, j] <- drop(partial_cov_gen(X1 = x[i,,drop = F], X2 = object$X0, theta = object$theta, i1 = 1, i2 = j, arg = "X_i_j", type = object$covtype)) * kvec[,i]
    
    dm[i,] <- crossprod(object$Z0 - object$beta0, object$Ki) %*% dkvec
    if((is(object, "hetGP") || is(object, "homGP")) && object$trendtype == "OK") tmp <- drop(1 - colSums(object$Ki) %*% kvec[,i])/(sum(object$Ki)) * colSums(object$Ki) %*% dkvec else tmp <- 0
    ds2[i,] <- -2 * (crossprod(kvec[,i], object$Ki) %*% dkvec + tmp)
  }
  
  if(is(object, "hetGP") || is(object, "homGP")){
    return(list(mean = dm, sd2 = object$nu_hat * ds2))
  }else{
    return(list(mean = object$sigma2 * dm, sd2 =  (object$nu + object$psi - 2) / (object$nu + length(object$Z) - 2) * object$sigma2^2 * ds2)) 
  }
}


#' @param object GP/TP object
#' @param x design location
#' @returns A list with elements:
#' \itemize{
#' \item \code{mean}: the gradient gp mean
#' \item \code{sd2}: the gradient gp predictive variance;
#' }
#' @seealso \code{\link[hetGP]{predict_gr}} for the prediction for the derivative GP itself
#' @noRd
#' @note untested for TP objects. The cross covariance between x locations is not computed yet.
#' @references 
#' N. Wycoff, M. Binois, S. Wild (2021), Sequential Learning of Active Subspaces, JCGS \cr \cr
#' 
#' Wu, A.; Aoi, M. C. & Pillow, J. W. (2017), Exploiting gradients and Hessians in Bayesian optimization and Bayesian quadrature, arXiv preprint arXiv:1704.00060, 2017
predict_derivGP <- function(object, x){
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(object$trendtype != 'SK') warning("Only the Simple kriging case is treated.")
  
  kvec <-  cov_gen(X1 = object$X0, X2 = x, theta = object$theta, type = object$covtype)
  
  dm <- ds2 <- matrix(NA, nrow(x), ncol(x))
  
  
  for(i in 1:ncol(x)){
    if(object$covtype == "Gaussian") Km <- 2 / object$theta
    if(object$covtype == "Matern3_2") Km <- 3 / object$theta^2
    if(object$covtype == "Matern5_2") Km <- 5 / (3 * object$theta^2)
  }
  
  for(i in 1:nrow(x)){
    dkvec <- matrix(NA, nrow(object$X0), ncol(x))
    for(j in 1:ncol(x)) dkvec[, j] <- drop(partial_cov_gen(X1 = x[i,,drop = F], X2 = object$X0, theta = object$theta, i1 = 1, i2 = j, arg = "X_i_j", type = object$covtype)) * kvec[,i]
    
    dm[i,] <- crossprod(object$Z0 - object$beta0, object$Ki) %*% dkvec
    # if((is(object, "hetGP") || is(object, "homGP")) && object$trendtype == "OK") tmp <- drop(1 - colSums(object$Ki) %*% kvec[,i])/(sum(object$Ki)) * colSums(object$Ki) %*% dkvec else tmp <- 0
    tmp <- 0
    ds2[i,] <- as.vector(Km) - (fast_diag(crossprod(dkvec, object$Ki), dkvec) + tmp)
  }
  
  if(is(object, "hetGP") || is(object, "homGP")){
    return(list(mean = dm, sd2 = object$nu_hat * ds2))
  }else{
    return(list(mean = object$sigma2 * dm, sd2 =  (object$nu + object$psi - 2) / (object$nu + length(object$Z) - 2) * object$sigma2^2 * ds2)) 
  }  
}

#' Derivative of the student-t pdf
#' @param z input location
#' @param a degree of freedom parameter
#' @noRd
#' @examples
#' # grad(dt, x = 0.55, df = 3.6)
#' # dlambda(0.55, 3.6)
dlambda <- function(z, a){
  return(-(a + 1) * gamma((a + 1)/2)/(sqrt(pi * a) * a * gamma(a/2)) * z * ((a + z^2)/a)^(-(a +3)/2))
}

#' Fast approximated batch-Expected Improvement criterion (for minimization)
#' @title Parallel Expected improvement
#' @param x matrix of new designs representing the batch of q points,
#'  one point per row (size q x d)
#' @param model \code{homGP} or \code{hetGP} model, including inverse matrices.
#' @param cst optional plugin value used in the EI, see details
#' @param preds optional predictions at \code{x} to avoid recomputing if already done (must include the predictive covariance, i.e., the \code{cov} slot)
#' @details 
#' \code{cst} is classically the observed minimum in the deterministic case. 
#' In the noisy case, the min of the predictive mean works fine.
#' @note This is a beta version at this point. It may work for for TP models as well.
#' @references 
#' M. Binois (2015), Uncertainty quantification on Pareto fronts and high-dimensional strategies in Bayesian optimization, with applications in multi-objective automotive design.
#' Ecole Nationale Superieure des Mines de Saint-Etienne, PhD thesis.
#' @export
#' @importFrom stats cov2cor
#' @examples 
#' ## Optimization example (noiseless)
#' set.seed(42)
#' 
#' ## Test function defined in [0,1]
#' ftest <- f1d
#' 
#' n_init <- 5 # number of unique designs
#' X <- seq(0, 1, length.out = n_init)
#' X <- matrix(X, ncol = 1)
#' Z <- ftest(X)
#' 
#' ## Predictive grid
#' ngrid <- 51
#' xgrid <- seq(0,1, length.out = ngrid)
#' Xgrid <- matrix(xgrid, ncol = 1)
#' 
#' model <- mleHomGP(X = X, Z = Z, lower = 0.01, upper = 1, known = list(g = 2e-8))
#' 
#' # Regular EI function
#' cst <- min(model$Z0)
#' EIgrid <- crit_EI(Xgrid, model, cst = cst)
#' plot(xgrid, EIgrid, type = "l")
#' abline(v = X, lty = 2) # observations
#' 
#' # Create batch (based on regular EI peaks)
#' xbatch <- matrix(c(0.37, 0.17, 0.7), 3, 1)
#' abline(v = xbatch, col = "red")
#' fqEI <- crit_qEI(xbatch, model, cst = cst)
#' 
#' # Compare with Monte Carlo qEI
#' preds <- predict(model, xbatch, xprime = xbatch)
#' nsim <- 1e4
#' simus <- matrix(rnorm(3 * nsim), nsim) %*% chol(preds$cov)
#' simus <- simus + matrix(preds$mean, nrow = nsim, ncol = 3, byrow = TRUE)
#' MCqEI <- mean(apply(cst - simus, 1, function(x) max(c(x, 0))))
#' 
crit_qEI <- function(x, model, cst = NULL, preds = NULL){
  if(is.null(cst)) cst <- min(predict(model, x = model$X0)$mean)
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(preds) || is.null(preds$cov)) preds <- predict(model, x = x, xprime = x)
  
  cormat <- cov2cor(preds$cov)
  res <- qEI_cpp(mu = -preds$mean, s = sqrt(preds$sd2), cor = cormat, threshold = -cst)
  
  return(res) 
}


#' Search for best reduction in a criterion 
#' @title Criterion minimization
#' @param model \code{homGP} or \code{hetGP} model
#' @param crit considered criterion, one of \code{"crit_cSUR"}, \code{"crit_EI"}, \code{"crit_ICU"},
#'  \code{"crit_MCU"} and \code{"crit_tMSE"}. Note that \code{crit_IMSPE} has its dedicated method, see \code{\link[hetGP]{IMSPE_optim}}.
#' @param ... additional parameters of the criterion
#' @param replicate if \code{TRUE}, search only on existing designs
#' @param Xcand optional set of of candidates for discrete search
#' @param control list in case \code{Xcand == NULL}, with elements \code{multi.start},
#' to perform a multi-start optimization based on \code{\link[stats]{optim}}, with \code{maxit} iterations each.
#' Also, \code{tol_dist} defines the minimum distance to an existing design for a new point to be added, otherwise the closest existing design is chosen.
#' In a similar fashion, \code{tol_dist} is the minimum relative change of crit for adding a new design.
#' @param seed optional seed for the generation of designs with \code{\link[DiceDesign]{maximinSA_LHS}}
#' @param ncores number of CPU available (> 1 mean parallel TRUE), see \code{\link[parallel]{mclapply}}
#' @importFrom DiceDesign lhsDesign maximinSA_LHS
#' @return list with \code{par}, \code{value} elements, and additional slot \code{new} (boolean if it is or not a new design) and \code{id} giving the index of the duplicated design. 
#' @noRd
#' @examples 
#' ###############################################################################
#' ## Bi-variate example
#' ###############################################################################
#' 
#' nvar <- 2 
#' 
#' set.seed(42)
#' ftest <- function(x, coef = 0.1) return(sin(2*pi*sum(x)) + rnorm(1, sd = coef))
#' 
#' n <- 25 # must be a square
#' xgrid0 <- seq(0.1, 0.9, length.out = sqrt(n))
#' designs <- as.matrix(expand.grid(xgrid0, xgrid0))
#' X <- designs[rep(1:n, sample(1:10, n, replace = TRUE)),]
#' Z <- apply(X, 1, ftest)
#' 
#' model <- mleHomGP(X, Z, lower = rep(0.1, nvar), upper = rep(1, nvar))
#' 
#' ngrid <- 51
#' xgrid <- seq(0,1, length.out = ngrid)
#' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
#' 
#' preds <- predict(x = Xgrid, object =  model)
#' 
#' ## Initial plots
#' contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
#'         main = "Predicted mean", nlevels = 20)
#' points(model$X0, col = 'blue', pch = 20)
#' 
#' crit <- "crit_EI"
#' crit_grid <- apply(Xgrid, 1, crit, model = model)
#' filled.contour(x = xgrid, y = xgrid, matrix(crit_grid, ngrid),
#'                nlevels = 20, color.palette = terrain.colors, 
#'                main = "Initial criterion landscape",
#' plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
#' 
#' ## Sequential crit search
#' nsteps <- 1 # Increase for better results
#' 
#' for(i in 1:nsteps){
#'   res <- crit.search(model, crit = crit, control = list(multi.start = 100, maxit = 50))
#'   newX <- res$par
#'   newZ <- ftest(newX)
#'   model <- update(object = model, Xnew = newX, Znew = newZ)
#' }
#' 
#' ## Final plots
#' contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
#'         main = "Predicted mean", nlevels = 20)
#' points(model$X0, col = 'blue', pch = 20)
#' 
#' crit_grid <- apply(Xgrid, 1, crit, model = model)
#' filled.contour(x = xgrid, y = xgrid, matrix(crit_grid, ngrid),
#'                nlevels = 20, color.palette = terrain.colors, 
#'                main = "Final criterion landscape",
#' plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
#' 
crit.search <- function(model, crit, ..., replicate = FALSE, Xcand = NULL, 
                        control = list(tol_dist = 1e-6, tol_diff = 1e-6, multi.start = 20,
                                       maxit = 100, maximin = TRUE, Xstart = NULL), seed = NULL,
                        ncores = 1){
  # Only search on existing designs
  if(replicate){
    ## Discrete optimization
    res <- mclapply(1:nrow(model$X0), 
                    function(i) match.fun(crit)(x = model$X0[i,,drop = FALSE], model = model, ... = ...), mc.cores = ncores)
    res <- unlist(res)
    
    return(list(par = model$X0[which.max(res),,drop = FALSE], value = max(res), new = FALSE, id = which.max(res)))
  }
  
  if(is.null(control))
    control <- list(multi.start = 20, maxit = 100)
  
  if(is.null(control$multi.start))
    control$multi.start <- 20
  
  if(is.null(control$maxit))
    control$maxit <- 100
  
  if(is.null(control$maximin))
    control$maximin <- TRUE
  
  if(is.null(control$tol_dist)) control$tol_dist <- 1e-6
  if(is.null(control$tol_diff)) control$tol_diff <- 1e-6
  
  d <- ncol(model$X0)
  
  if(crit == "crit_EI") gr <- deriv_crit_EI else gr <- NULL
  crit <- match.fun(crit)
  
  ## Optimization
  if(is.null(Xcand)){
    ## Continuous optimization
    if(!is.null(control$Xstart)){
      Xstart <- control$Xstart
    }else{
      if(is.null(seed)) seed <- sample(1:2^15, 1) ## To be changed?
      if(control$maximin){
        if(d == 1){
          # perturbed 1d equidistant points
          Xstart <- matrix(seq(1/2, control$multi.start -1/2, length.out = control$multi.start) + runif(control$multi.start, min = -1/2, max = 1/2), ncol = 1)/control$multi.start
        }else{
          Xstart <- maximinSA_LHS(lhsDesign(control$multi.start, d, seed = seed)$design)$design
        }
      }else{
        Xstart <- lhsDesign(control$multi.start, d, seed = seed)$design
      }
    }
    
    res <- list(par = NA, value = -Inf, new = NA)
    
    local_opt_fun <- function(i){
      out <- try(optim(Xstart[i,, drop = FALSE], crit, ... = ..., method = "L-BFGS-B", gr = gr, 
                       lower = rep(0, d), upper = rep(1, d),
                       model = model, control = list(maxit = control$maxit, fnscale = -1)))
      if(is(out, "try-error")) return(NULL)
      return(out)
    }
    
    all_res <- mclapply(1:nrow(Xstart), local_opt_fun, mc.cores = ncores)
    res_max <- which.max(Reduce(c, lapply(all_res, function(x) x$value)))
    res <- list(par = apply(all_res[[res_max]]$par, c(1, 2), function(x) max(min(x , 1), 0)),
                value = all_res[[res_max]]$value, new = TRUE, id = NULL)
    
    # for(i in 1:nrow(Xstart)){
    #   out <- try(optim(Xstart[i,, drop = FALSE], crit, ... = ..., method = "L-BFGS-B", gr = gr, 
    #                    lower = rep(0, d), upper = rep(1, d),
    #                    model = model, control = list(maxit = control$maxit, fnscale = -1)))
    #   if(class(out) != "try-error"){
    #     if(out$value > res$value)
    #       res <- list(par = out$par, value = out$value, new = TRUE, id = NULL)
    #   }
    # }
    
    if(control$tol_dist > 0 || control$tol_diff > 0){
      ## Check if new design is not to close to existing design
      dists <- sqrt(distance_cpp(res$par, model$X0))
      if(min(dists) < control$tol_dist){
        res <- list(par = model$X0[which.min(dists),,drop = FALSE],
                    value = crit(x = model$X0[which.min(dists),, drop = F], model = model, ... = ...),
                    new = FALSE, id = which.min(dists))
      }else{
        ## Check if crit difference between replication and new design is significative
        id_closest <- which.min(dists) # closest point to new design
        crit_rep <- crit(model$X0[which.min(dists),,drop = FALSE], model = model, ...=...)
        
        ## EI can be 0, in which case it is better to replicate even if crit_rep$value is also 0
        if(res$value == 0 || (res$value - crit_rep)/res$value < control$tol_diff){
          res <- list(par = model$X0[which.min(dists),,drop = FALSE],
                      value = crit_rep,
                      new = FALSE, id = which.min(dists))
        }
      }
    }
    
    
    return(res)
    
    
  }else{
    ## Discrete optimization
    res <- mclapply(1:nrow(model$X0), 
                    function(i) match.fun(crit)(x = Xcand[i,,drop = FALSE], model = model, ... = ...), mc.cores = ncores)
    res <- unlist(res)
    
    tmp <- which(duplicated(rbind(model$X0, Xcand[which.max(res),,drop = FALSE]), fromLast = TRUE))
    if(length(tmp) > 0) return(list(par = Xcand[which.max(res),,drop = FALSE], value = max(res), new = FALSE, id = tmp))
    return(list(par = Xcand[which.max(res),,drop = FALSE], value = max(res), new = TRUE, id = NULL))
  }
}


#' Search for the best value of available criterion, possibly using a h-steps lookahead strategy to favor designs with replication
#' @title Criterion optimization
#' @param model \code{homGP} or \code{hetGP} model
#' @param crit considered criterion, one of \code{"crit_cSUR"}, \code{"crit_EI"}, \code{"crit_ICU"},
#'  \code{"crit_MCU"} and \code{"crit_tMSE"}. Note that \code{crit_IMSPE} has its dedicated method, see \code{\link[hetGP]{IMSPE_optim}}.
#' @param ... additional parameters of the criterion
##' @param Xcand optional discrete set of candidates (otherwise a maximin LHS is used to initialize continuous search)
#' @param control list in case \code{Xcand == NULL}, with elements \code{multi.start},
#' to perform a multi-start optimization based on \code{\link[stats]{optim}}, with \code{maxit} iterations each.
#' Also, \code{tol_dist} defines the minimum distance to an existing design for a new point to be added, otherwise the closest existing design is chosen.
#' In a similar fashion, \code{tol_dist} is the minimum relative change of crit for adding a new design.
#' @param seed optional seed for the generation of LHS designs with \code{\link[DiceDesign]{maximinSA_LHS}}
#' @param h horizon for multi-step ahead framework.
#' The decision is made between:
#' \itemize{
#'  \item sequential crit search starting by a new design (optimized first) then adding \code{h} replicates
#'  \item sequential crit searches starting by \code{1} to \code{h} replicates before adding a new point
#' }
#' Use \code{h = 0} for the myopic criterion, i.e., not looking ahead.
#' @param ncores number of CPU available (> 1 mean parallel TRUE), see \code{\link[parallel]{mclapply}}
#' @details 
#' When looking ahead, the kriging believer heuristic is used,
#'  meaning that the non-observed value is replaced by the mean prediction in the update.  
#' @return list with elements:
#' \itemize{
#' \item \code{par}: best first design,
#' \item \code{value}: criterion h-steps ahead starting from adding \code{par},
#' \item \code{path}: list of elements list(\code{par}, \code{value}, \code{new}) at each step \code{h}
#' }
#' @references
#' M. Binois, J. Huang, R. B. Gramacy, M. Ludkovski (2019), 
#' Replication or exploration? Sequential design for stochastic simulation experiments,
#' Technometrics, 61(1), 7-23.\cr 
#' Preprint available on arXiv:1710.03206.
#' @export 
#' @examples
#' ###############################################################################
#' ## Bi-variate example (myopic version)
#' ###############################################################################
#' 
#' nvar <- 2 
#' 
#' set.seed(42)
#' ftest <- function(x, coef = 0.1) return(sin(2*pi*sum(x)) + rnorm(1, sd = coef))
#' 
#' n <- 25 # must be a square
#' xgrid0 <- seq(0.1, 0.9, length.out = sqrt(n))
#' designs <- as.matrix(expand.grid(xgrid0, xgrid0))
#' X <- designs[rep(1:n, sample(1:10, n, replace = TRUE)),]
#' Z <- apply(X, 1, ftest)
#' 
#' model <- mleHomGP(X, Z, lower = rep(0.1, nvar), upper = rep(1, nvar))
#' 
#' ngrid <- 51
#' xgrid <- seq(0,1, length.out = ngrid)
#' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
#' 
#' preds <- predict(x = Xgrid, object =  model)
#' 
#' ## Initial plots
#' contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
#'         main = "Predicted mean", nlevels = 20)
#' points(model$X0, col = 'blue', pch = 20)
#' 
#' crit <- "crit_EI"
#' crit_grid <- apply(Xgrid, 1, crit, model = model)
#' filled.contour(x = xgrid, y = xgrid, matrix(crit_grid, ngrid),
#'                nlevels = 20, color.palette = terrain.colors, 
#'                main = "Initial criterion landscape",
#' plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
#' 
#' ## Sequential crit search
#' nsteps <- 1 # Increase for better results
#' 
#' for(i in 1:nsteps){
#'   res <- crit_optim(model, crit = crit, h = 0, control = list(multi.start = 50, maxit = 30))
#'   newX <- res$par
#'   newZ <- ftest(newX)
#'   model <- update(object = model, Xnew = newX, Znew = newZ)
#' }
#' 
#' ## Final plots
#' contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
#'         main = "Predicted mean", nlevels = 20)
#' points(model$X0, col = 'blue', pch = 20)
#' 
#' crit_grid <- apply(Xgrid, 1, crit, model = model)
#' filled.contour(x = xgrid, y = xgrid, matrix(crit_grid, ngrid),
#'                nlevels = 20, color.palette = terrain.colors, 
#'                main = "Final criterion landscape",
#' plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
#' 
#' ###############################################################################
#' ## Bi-variate example (look-ahead version)
#' ###############################################################################
#' \dontrun{  
#' nvar <- 2 
#' 
#' set.seed(42)
#' ftest <- function(x, coef = 0.1) return(sin(2*pi*sum(x)) + rnorm(1, sd = coef))
#' 
#' n <- 25 # must be a square
#' xgrid0 <- seq(0.1, 0.9, length.out = sqrt(n))
#' designs <- as.matrix(expand.grid(xgrid0, xgrid0))
#' X <- designs[rep(1:n, sample(1:10, n, replace = TRUE)),]
#' Z <- apply(X, 1, ftest)
#' 
#' model <- mleHomGP(X, Z, lower = rep(0.1, nvar), upper = rep(1, nvar))
#' 
#' ngrid <- 51
#' xgrid <- seq(0,1, length.out = ngrid)
#' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
#' 
#' nsteps <- 5 # Increase for more steps
#' crit <- "crit_EI"
#' 
#' # To use parallel computation (turn off on Windows)
#' library(parallel)
#' parallel <- FALSE #TRUE #
#' if(parallel) ncores <- detectCores() else ncores <- 1
#' 
#' for(i in 1:nsteps){
#'   res <- crit_optim(model, h = 3, crit = crit, ncores = ncores,
#'                     control = list(multi.start = 100, maxit = 50))
#'   
#'   # If a replicate is selected
#'   if(!res$path[[1]]$new) print("Add replicate")
#'   
#'   newX <- res$par
#'   newZ <- ftest(newX)
#'   model <- update(object = model, Xnew = newX, Znew = newZ)
#'   
#'   ## Plots 
#'   preds <- predict(x = Xgrid, object =  model)
#'   contour(x = xgrid,  y = xgrid, z = matrix(preds$mean, ngrid),
#'           main = "Predicted mean", nlevels = 20)
#'   points(model$X0, col = 'blue', pch = 20)
#'   points(newX, col = "red", pch = 20)
#'   
#'   crit_grid <- apply(Xgrid, 1, crit, model = model)
#'   filled.contour(x = xgrid, y = xgrid, matrix(crit_grid, ngrid),
#'                  nlevels = 20, color.palette = terrain.colors,
#'   plot.axes = {axis(1); axis(2); points(model$X0, pch = 20)})
#' }
#' }
crit_optim <- function(model, crit, ..., h = 2, Xcand = NULL, control = list(multi.start = 10, maxit = 100), seed = NULL, ncores = 1){
  d <- ncol(model$X0)
  
  if(crit == "crit_IMSPE") stop("crit_IMSPE is intended to be optimized by IMSPE_optim")
  
  ## A) Setting to beat: first new point then replicate h times
  crit_A <- crit.search(model = model, crit = crit, ... = ..., control = control, Xcand = Xcand, seed = seed, ncores = ncores)
  new_designA <- crit_A$par ## store first considered design to be added
  path_A <- list(crit_A)
  
  if(h > 0){
    newmodelA <- model
    for(i in 1:h){
      ZnewA <- predict(newmodelA, crit_A$par)$mean
      newmodelA <- update(object = newmodelA, Xnew = crit_A$par, Znew = ZnewA, maxit = 0)
      crit_A <- crit.search(model = newmodelA, crit = crit, ... = ..., replicate = TRUE, control = control, seed = seed, ncores = ncores)
      path_A <- c(path_A, list(crit_A))
    }
  }
  
  if(h == -1) return(list(par = new_designA, value = crit_A$value, path = path_A)) 
  
  ## B) Now compare with waiting to add new point
  newmodelB <- model
  
  if(h == 0){
    crit_B <- crit.search(model = newmodelB, crit = crit, ... =..., replicate = TRUE, control = control, seed = seed, ncores = ncores)
    new_designB <- crit_B$par ## store considered design to be added
    
    # search from best replicate
    if(is.null(Xcand)){
      crit_C <- crit.search(model = newmodelB, crit = crit,
                            control = list(Xstart = crit_B$par, maxit = control$maxit,
                                           tol_dist = control$tol_dist, tol_diff = control$tol_diff), seed = seed, ncores = ncores)
    }else{
      crit_C <- crit_B
    }
    
    if(crit_C$value > max(crit_A$value, crit_B$value)) return(list(par = crit_C$par, value = crit_C$value, path = list(crit_C)))
    
    if(crit_B$value > crit_A$value){
      return(list(par = crit_B$par, value = crit_B$value, path = list(crit_B)))
    } 
  }else{
    for(i in 1:h){
      ## Add new replicate
      crit_B <- crit.search(model = newmodelB, crit = crit, ... = ..., replicate = TRUE, control = control, seed = seed, ncores = ncores)
      
      if(i == 1){
        new_designB <- matrix(crit_B$par, nrow = 1) ##store first considered design to add
        path_B <- list()
      } 
      
      path_B <- c(path_B, list(crit_B))
      ZnewB <- predict(newmodelB, crit_B$par)$mean
      newmodelB <- update(object = newmodelB, Xnew = crit_B$par, Znew = ZnewB, maxit = 0)
      
      ## Add new design
      crit_C <- crit.search(model = newmodelB, crit = crit, ... = ..., control = control, Xcand = Xcand, seed = seed, ncores = ncores)
      path_C <- list(crit_C)
      
      if(i < h){
        newmodelC <- newmodelB
        
        for(j in i:(h-1)){
          ## Add remaining replicates
          ZnewC <- predict(newmodelC, crit_C$par)$mean
          newmodelC <- update(object = newmodelC, Xnew = crit_C$par, Znew = ZnewC, maxit = 0)
          crit_C <- crit.search(model = newmodelC, crit = crit, ... = ..., replicate = TRUE, control = control, seed = seed, ncores = ncores)
          path_C <- c(path_C, list(crit_C))
        }
      }
      
      if(crit_C$value < crit_A$value) return(list(par = new_designB, value = crit_C$value, path = c(path_B, path_C)))
    }
  }
  
  return(list(par = new_designA, value = crit_A$value, path = path_A))
  
}


#' log1mexp function
#' @references Maechler, Martin (2012). Accurately Computing log(1-exp(-|a|)). Assessed from the Rmpfr package.
#' @noRd
log1mexp <- function(x){
  if (is.na(x)){
    warnings("x is NA")
    return(NA)
  }
  if(any(x < 0)){
    warnings("x < 0")
    return(NA)
  }
  
  ifelse(x <= log(2), log(-expm1(-x)), log1p(-exp(-x)))
}

erfc <- function(x){
  2 * pnorm(x * sqrt(2), lower.tail = FALSE)
} 

erfcx <- function (x) 
{
  exp(x^2) * erfc(x)
}

#' log_h function from 
#' @references Ament, S., Daulton, S., Eriksson, D., Balandat, M., & Bakshy, E. (2024). Unexpected improvements to expected improvement for bayesian optimization. Advances in Neural Information Processing Systems, 36.
#' @noRd
log_h <- function(z, eps = .Machine$double.eps){
  c1 <- log(2*pi)/2
  c2 <- log(pi/2)/2
  if(z > -1) return(log(dnorm(z) + z * pnorm(z)))
  if(z < -1 / sqrt(eps)) return(-z^2/2 - c1 - 2 * log(abs(z)))
  res <- -z^2/2 - c1 + log1mexp(-(log(erfcx(-z/sqrt(2)) * abs(z)) + c2))
  if(is.na(res)) res <- -750
  return(res)
}

#' Computes log of EI for minimization, with improved stability with respect to EI
#' @title Logarithm of Expected Improvement criterion
#' @param x matrix of new designs, one point per row (size n x d)
#' @param model \code{homGP} or \code{hetGP} model, or their TP equivalents, including inverse matrices. For TP models, the computation is using the one from regular EI.
#' @param cst optional plugin value used in the EI, see details
#' @param preds optional predictions at \code{x} to avoid recomputing if already done
#' @note This is a beta version at this point.
#' @seealso \code{\link[hetGP]{crit_EI}} for the regular EI criterion and compare the outcomes
#' @details 
#' \code{cst} is classically the observed minimum in the deterministic case. 
#' In the noisy case, the min of the predictive mean works fine.
#'  
#' @references
#' Ament, S., Daulton, S., Eriksson, D., Balandat, M., & Bakshy, E. (2024). Unexpected improvements to expected improvement for Bayesian optimization. Advances in Neural Information Processing Systems, 36.
#' @importFrom stats dt qt
#' @export
#' @examples 
#' ## Optimization example
#' set.seed(42)
#' 
#' 
#' ## Noise field via standard deviation
#' noiseFun <- function(x, coef = 1.1, scale = 1){
#' if(is.null(nrow(x)))
#'  x <- matrix(x, nrow = 1)
#'    return(scale*(coef + cos(x * 2 * pi)))
#' }
#' 
#' ## Test function defined in [0,1]
#' ftest <- function(x){
#' if(is.null(nrow(x)))
#' x <- matrix(x, ncol = 1)
#' return(f1d(x) + rnorm(nrow(x), mean = 0, sd = noiseFun(x)))
#' }
#' 
#' n_init <- 10 # number of unique designs
#' N_init <- 100 # total number of points
#' X <- seq(0, 1, length.out = n_init)
#' X <- matrix(X[sample(1:n_init, N_init, replace = TRUE)], ncol = 1)
#' Z <- ftest(X)
#' 
#' ## Predictive grid
#' ngrid <- 51
#' xgrid <- seq(0,1, length.out = ngrid)
#' Xgrid <- matrix(xgrid, ncol = 1)
#' 
#' model <- mleHetGP(X = X, Z = Z, lower = 0.001, upper = 1)
#' 
#' logEIgrid <- crit_logEI(Xgrid, model)
#' preds <- predict(x = Xgrid, model)
#' 
#' par(mar = c(3,3,2,3)+0.1)
#' plot(xgrid, f1d(xgrid), type = 'l', lwd = 1, col = "blue", lty = 3,
#' xlab = '', ylab = '', ylim = c(-8,16))
#' points(X, Z)
#' lines(Xgrid, preds$mean, col = 'red', lwd = 2)
#' lines(Xgrid, qnorm(0.05, preds$mean, sqrt(preds$sd2)), col = 2, lty = 2)
#' lines(Xgrid, qnorm(0.95, preds$mean, sqrt(preds$sd2)), col = 2, lty = 2)
#' lines(Xgrid, qnorm(0.05, preds$mean, sqrt(preds$sd2 + preds$nugs)), col = 3, lty = 2)
#' lines(Xgrid, qnorm(0.95, preds$mean, sqrt(preds$sd2 + preds$nugs)), col = 3, lty = 2)
#' par(new = TRUE)
#' plot(NA, NA, xlim = c(0, 1), ylim = range(logEIgrid), axes = FALSE, ylab = "", xlab = "")
#' lines(xgrid, logEIgrid, lwd = 2, col = 'cyan')
#' axis(side = 4)
#' mtext(side = 4, line = 2, expression(logEI(x)), cex = 0.8)
#' mtext(side = 2, line = 2, expression(f(x)), cex = 0.8)
crit_logEI <- function(x, model, cst = NULL, preds = NULL){
  if(is.null(cst)) cst <- min(predict(model, x = model$X0)$mean)
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(preds)) preds <- predict(model, x = x)
  
  if(is(model, "homTP") || is(model, "hetTP")){
    gamma <- (cst - preds$mean)/sqrt(preds$sd2)
    res <- (cst - preds$mean) * pt(gamma, df = model$nu + length(model$Z))
    res <- res + sqrt(preds$sd2) * (1 + (gamma^2 - 1)/(model$nu + length(model$Z) - 1)) * dt(x = gamma, df = model$nu + length(model$Z))
    res[which(res < 1e-12)] <- 1e-14 # for stability
    return(log(res))
  }
  
  xcr <- (cst - preds$mean)/sqrt(preds$sd2)
  res <- sapply(xcr, log_h) + log(sqrt(preds$sd2))
  res[preds$sd2 == 0] <- -800
  # res <- (cst - preds$mean) * pnorm(xcr)
  # res <- res + sqrt(preds$sd2) * dnorm(xcr)
  # res[which(preds$sd2 < sqrt(.Machine$double.eps) | res < 1e-12)] <- 0 # for stability
  return(res) 
}


#' derivative of log_h function from 
#' @references Ament, S., Daulton, S., Eriksson, D., Balandat, M., & Bakshy, E. (2024). Unexpected improvements to expected improvement for bayesian optimization. Advances in Neural Information Processing Systems, 36.
#' @noRd
#' @examples
#' # Verification
#' z <- 1.5
#' eps <- 1e-4
#' hetGP:::dlog_h(z)
#' (hetGP:::log_h(z+1e-4) - hetGP:::log_h(z))*1e4
#' 
#' z <- -110.3
#' hetGP:::dlog_h(z, eps = eps)
#' (hetGP:::log_h(z+1e-4, eps = eps) - hetGP:::log_h(z, eps = eps))*1e4
#' 
#' # d(erfcx(z)) = 2 * z * exp(z^2)*erfc(z) - 2/sqrt(pi)
#' z <- -2.1
#' 2 * z * exp(z^2) * hetGP:::erfc(z) - 2/sqrt(pi)
#' (hetGP:::erfcx(z+1e-4) - hetGP:::erfcx(z))*1e4
#' 
#' # d(erfcx(-z))
#' -2 * z * exp(z^2) * hetGP:::erfc(-z) - 2/sqrt(pi)
#' (hetGP:::erfcx(-z+1e-4) - hetGP:::erfcx(-z))*1e4
#' 
#' # d(erfcx(-z/sqrt(2)))
#' z * exp(z^2/2) * (2 - hetGP:::erfc(z/sqrt(2))) + sqrt(2/pi)
#' (hetGP:::erfcx(-(z + 1e-4)/sqrt(2)) - hetGP:::erfcx(-z/sqrt(2)))*1e4
#' 
#' # d(log(erfcx(-z/sqrt(2)) * -z))
#' ((z^2 + 1) * exp(z^2 / 2) * (hetGP:::erfc(z / sqrt(2)) - 2) - (sqrt(2) * z) / sqrt(pi))/(-hetGP:::erfcx(-z/sqrt(2)) * z)
#' ft <- function(z) return(log(hetGP:::erfcx(-z/sqrt(2)) * abs(z)))
#' (ft(z + 1e-4) - ft(z))*1e4
#' 
#' # d logm1exp(...)
#' ft <- function(z) return(-z^2/2 + hetGP:::log1mexp(-(log(hetGP:::erfcx(-z/sqrt(2)) * abs(z)) + log(pi/2)/2)))
#' (ft(z + 1e-4) - ft(z))*1e4
#' -z + (exp(log(pi/2)/2) * (sqrt(pi) * (z^2 + 1) * exp(z^2 / 2) * (hetGP:::erfc(z / sqrt(2)) - 2) - sqrt(2) * z)) / (sqrt(pi) * (z * exp((z^2 + 2 * log(pi/2)/2) / 2) * (hetGP:::erfc(z / sqrt(2)) - 2) - 1))
#' 
#' # Now more stable
#' (ft(z + 1e-4) - ft(z))*1e4
#' -z + (exp(log(pi/2)/2) * (-sqrt(pi) * (z^2 + 1) * hetGP:::erfcx(-z/sqrt(2)) - sqrt(2) * z)) / (sqrt(pi) * (z * exp(log(pi/2)/2) * -hetGP:::erfcx(-z / sqrt(2)) - 1))
#' 
#' hetGP:::dlog_h(z, eps = eps)
#' (hetGP:::log_h(z+1e-4, eps = eps) - hetGP:::log_h(z, eps = eps))*1e4
dlog_h <- function(z, eps = .Machine$double.eps){
  c2 <- log(pi/2)/2
  if(z > -1){
    # d log(dnorm(z) + z * pnorm(z)) = (pnorm(z))/(dnorm(z) + z * pnorm(z))
    return(pnorm(z)/(dnorm(z) + z * pnorm(z)))
  } 
  if(z < -1 / sqrt(eps)){
    return(-z - 2 / z)
  }
  # d(-z^2/2 -c1 + log1mexp(-(log(erfcx(-z/sqrt(2)) * abs(z)) + c2))) = 
  res <- -z + (exp(c2) * (-sqrt(pi) * (z^2 + 1) * erfcx(-z/sqrt(2)) - sqrt(2) * z)) / (sqrt(pi) * (z * exp(c2) * -erfcx(-z / sqrt(2)) - 1))
  
  if(is.na(res)) res <- -z - 2 / z # Use more stable regime
  return(res)
}


#' Derivative of logEI criterion for GP models
#' @param x matrix for the new design (size 1 x d)
#' @param model \code{homGP} or \code{hetGP} model
#' @param cst threshold for contour criteria
#' @param preds pre-computed preds for contour criteria
#' @export
#' @seealso \code{\link[hetGP]{crit_logEI}} for the criterion
#' @references 
#' Ginsbourger, D. Multiples metamodeles pour l'approximation et l'optimisation de fonctions numeriques multivariables Ecole Nationale Superieure des Mines de Saint-Etienne, Ecole Nationale Superieure des Mines de Saint-Etienne, 2009 \cr \cr
#' Roustant, O., Ginsbourger, D., DiceKriging, DiceOptim: Two R packages for the analysis of computer experiments by kriging-based metamodeling and optimization, Journal of Statistical Software, 2012
#' Ament, S., Daulton, S., Eriksson, D., Balandat, M., & Bakshy, E. (2024). Unexpected improvements to expected improvement for Bayesian optimization. Advances in Neural Information Processing Systems, 36.
deriv_crit_logEI <- function(x, model, cst = NULL, preds = NULL){
  if(is.null(cst)) cst <- min(predict(model, x = model$X)$mean)
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(preds)) preds <- predict(model, x = x)
  
  pred_gr <- predict_gr(model, x)
  
  z <- (cst - preds$mean)/sqrt(preds$sd2)
  
  if(is(model, "homGP") || is(model, "hetGP")){
    
    # (sqrt(predict(model, x = x + c(1e-4,0))$sd2) - sqrt(predict(model, x = x)$sd2))*1e4
    # (sqrt(predict(model, x = x + c(0,1e-4))$sd2) - sqrt(predict(model, x = x)$sd2))*1e4
    # ((cst - predict(model, x = x +  c(1e-4,0))$mean)/sqrt(predict(model, x = x + 1e-4)$sd2) - z)*1e4
    # ((cst - predict(model, x = x +  c(0,1e-4))$mean)/sqrt(predict(model, x = x + 1e-4)$sd2) - z)*1e4
    # 
    # grad(function(x){sqrt(predict(model, x = x)$sd2)}, x = x)
    # grad(function(x){predict(model, x = x)$sd2}, x = x)
    # grad(function(x){(cst - predict(model, x = x)$mean)/sqrt(predict(model, x = x)$sd2)}, x = x)
    
    # dz = - dm/s - z ds/s = -dm/s - z * ds2/(2s2) 
    ds <- pred_gr$sd2/(2 * matrix(sqrt(preds$sd2), nrow = nrow(x), ncol(x)))
    dz <- -pred_gr$mean/sqrt(preds$sd2) - z * ds/sqrt(preds$sd2)
    res <- dz * dlog_h(z) + ds / sqrt(preds$sd2)
    
  }else{
    # dz = - dm/s - z ds/s = -dm/s - z * ds2/(2s2) 
    dz <- -pred_gr$mean/sqrt(preds$sd2) - z * pred_gr$sd2/(2 * matrix(preds$sd2, nrow = nrow(x), ncol(x)))
    
    # d( (cst - m(x)).pt(z(x))) = 
    p1 <- -pred_gr$mean * pt(z, df = model$nu + length(model$Z)) + (cst - preds$mean) * dz * dt(z, df = model$nu + length(model$Z))
    
    a <- model$nu + length(model$Z) - 1
    # d( s(x) (1 + (z^2-1)/(nu + N -1)) dt(z(x)) (in 2 lines)
    p2 <- (pred_gr$sd2/(2*sqrt(preds$sd2)) * (1 + (z^2 - 1)/a) + 2*sqrt(preds$sd2) * z * dz/a) * dt(z, df = model$nu + length(model$Z))
    p2 <- p2 + sqrt(preds$sd2) * (1 + (z^2 - 1)/a) * dz * dlambda(z, model$nu + length(model$Z))
    res <- p1 + p2
    res[which(abs(res) < 1e-12)] <- 0 # for stability with optim
    
    eitmp <- (cst - preds$mean) * pt(z, df = model$nu + length(model$Z))
    eitmp <- eitmp + sqrt(preds$sd2) * (1 + (z^2 - 1)/(model$nu + length(model$Z) - 1)) * dt(x = z, df = model$nu + length(model$Z))
    
    res <- res / eitmp
  }
  
  return(res)
}

#' Bayesian optimization loop with parallel EI starting from initial observations
#' @title BO loop with qEI
#' @param X0 initial design of experiments matrix
#' @param Y0 initial vector of responses at \code{X0}
#' @param model \code{homGP} or \code{hetGP} model
#' @param q batch size
#' @param nrep number of replicates at the q points, default to 1
#' @param fun test function to minimize
#' @param budget optimization budget
#' @param lower,upper domain bounds
#' @param control list with parameters 
#' \itemize{
#' \item nunif: number of uniformly sampled candidates points for acquisition function optimization
#' \item ncb: number of initial candidate batches
#' \item maxit: passed to \code{\link[stats]{optim}}
#' \item trackbest: if \code{TRUE}, the best estimated solution will be saved at each iteration
#' \item maxtime: alternative stopping criterion
#' }
#' @return  A list with components:
#' \itemize{
#' \item \code{par}: all points evaluated,
#' \item \code{value}: the matrix of objective values at the points given in \code{par},
#' \item \code{model}: the last kriging models fitted.
#' \item \code{membest}: a matrix of best estimated designs at each iteration.
#' \item \code{estbest}: corresponding predicted mean values.
#' }
#' @export
#' @importFrom methods is
#' @examples
#' d <- 2
#' n <- 10*d
#' N <- 5*n
#' budget <- 130 # Increase for better results
#' 
#' ## Noise field via standard deviation
#' noiseFun <- function(x){
#'   if(is.null(nrow(x)))
#'     x <- matrix(x, nrow = 1)
#'   return(1/5*(3*(2 + 2*sin(x[,1]*pi)*cos(x[,2]*3*pi) + 5*rowSums(x^2))))
#' }
#' 
#' ## Branin redefined in [0,1]^2
#' branin <- function(x){
#'  if(is.null(nrow(x))) x <- matrix(x, nrow = 1)
#'   x1 <- x[,1] * 15 - 5
#'   x2 <- x[,2] * 15
#'   return((x2 - 5/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(x1) + 10)
#'  }
#' 
#' ## data generating function combining mean and noise fields
#' ftest <- function(x){
#'   if(is.null(nrow(x))) x <- matrix(x, nrow = 1)
#'   return(branin(x) + rnorm(nrow(x), mean = 0, sd = noiseFun(x)))
#' }
#'
#'
#' ngrid <- 51
#' Xgrid <- as.matrix(expand.grid(seq(0,1,length.out = ngrid), seq(0,1,length.out = ngrid)))
#' Ygrid <- branin(Xgrid)
#' Ygrid <- cbind(Ygrid, noiseFun(Xgrid)^2)
#' x1 <- c(0.9616520, 0.15); x2 <- c(0.1238946, 0.8166644); x3 <- c(0.5427730, 0.15)
#' fstar <- 0.3978874
#'
#'
#' X0 <- matrix(runif(n*d),n, d)
#' X0 <- X0[sample(1:n, size = N, replace = TRUE),]
#' Y0 <- ftest(X0)
#' 
#' mod <- mleHetGP(X0, Y0, covtype = "Matern5_2", known = list(beta0 = 0))
#' 
#' opt <- qEI_loop(X0, Y0, mod, q = 10, nrep = 1, fun = ftest, budget = budget,
#'   lower = rep(0, d), upper = rep(1, d))
#' est <- predict(opt$model, opt$model$X0)$mean
#' xbest <- opt$model$X0[which.min(est),,drop=FALSE]
#' par(mfrow = c(1, 2))
#' contour(matrix(Ygrid[,1], ngrid), nlevels = 21, 
#'  main = "True function")
#' points(opt$model$X0, col = 4, pch = 20, cex = opt$model$mult/10)
#' points(rbind(t(x1), t(x2), t(x3)), pch=20, col="red")
#' points(xbest, col = "pink", pch = 15)
#' contour(matrix(sqrt(Ygrid[,2]^2), ngrid), nlevels = 21,
#'  main = "True variance")
#' points(rbind(t(x1), t(x2), t(x3)), pch=20, col="red")
#' points(opt$model$X0, col = 4, pch = 20, cex = opt$model$mult/10)
#' points(xbest, col = "pink", pch = 15)
#' par(mfrow = c(1, 1))
qEI_loop <- function(X0, Y0, model, q, nrep=1, fun, budget, lower, upper, control = NULL){
  d <- ncol(X0)
  if(is.null(control)) control <- list()
  if(is.null(control$nunif)) control$nunif <- 100*d # initial candidate points
  if(is.null(control$ncb)) control$ncb <- 100*d # initial candidate batches
  if(is.null(control$maxit)) control$maxit <- 50 # max iteration for optim call
  if(is.null(control$trackbest)) control$trackbest <- TRUE
  if(is.null(control$maxtime)) control$maxtime <- Inf
  if(is(model, "hetGP")) modfun <- mleHetGP else modfun <- mleHomGP
  
  curtime <- Sys.time() # Inital time
  timings <- 0 # Store times
  nuniq <- NULL # Track number of unique designs over time
  
  if(control$trackbest){
    Xbests <- X0[which.min(Y0),,drop = F]
    Ybests <- min(Y0)
  } else{
    Xbests <- Ybests <- NULL
  }
  
  
  while(nrow(X0) < budget && timings[length(timings)] < control$maxtime){
    if(budget - nrow(X0) < q) q <- budget - nrow(X0) # To stay within budget limits
    Xs <- matrix(runif(control$nunif * d), control$nunif) %*% diag(upper - lower) + matrix(lower, control$nunif, d, byrow = TRUE)
    tmp <- predict(model, model$X0)$mean
    cst <- min(tmp)
    
    # Add few samples around current mean for exploitation
    localpert <- pmin(upper, pmax(lower, X0[which.min(tmp),] + rnorm(d*q, sd = 0.05)))
    localpert <- matrix(localpert, q, byrow = TRUE)
    Xs <- rbind(Xs, localpert)
    
    # Add existing designs to favor replication
    Xs <- unique(rbind(Xs, X0))
    
    EIs <- apply(Xs, 1, crit_EI, cst = cst, model = model)
    
    solEI <- try(optim(Xs[which.max(EIs),], fn = crit_EI, gr = deriv_crit_EI, model = model, cst = cst,
                       lower = lower, upper = upper, method = "L-BFGS-B", control = list(fnscale = -1)))
    if(is(solEI, "try-error")){
      solEI <- list(par = Xs[which.max(EIs),,drop=F], value = max(EIs))
    }
    Xs <- rbind(Xs, solEI$par)
    EIs <- c(EIs, solEI$value)
    if(sum(EIs > 1e-8) < q) EIs <- EIs + 1e-8  # In case the search failed
    
    Xbcands <- array(NA, dim = c(control$ncb, q, d))
    for(i in 1:control$ncb) Xbcands[i,,] <- Xs[sample(1:(nrow(Xs)), q, prob = pmax(0, EIs)),]
    
    # More robust versions
    qEI2 <- function(x, ...){
      tmp <- try(crit_qEI(matrix(x, q), ...), silent = TRUE)
      if(is.nan(tmp)) return(-1)
      if(is(tmp, "try-error")) return(-1)
      return(tmp)
    }
    
    qEI_grid <- apply(Xbcands, 1, qEI2, model = model, cst = cst)
    
    Xq <- Xbcands[which.max(qEI_grid),,]
    
    oGPO <- try(optim(as.vector(Xq), qEI2, gr = NULL,
                      model = model, cst = cst,
                      lower = lower, upper = upper, method = "L-BFGS-B",
                      control = list(fnscale = -1, maxit = control$maxit)))
    
    if(is(oGPO,"try-error")){
      sds <- predict(model, Xs, type = "UK", light.return = TRUE)$sd
      Xnew <- Xs[order(sds, decreasing = TRUE)[1:q],,drop = FALSE]
      print("Optimisation failed, took max variance")
    } else {
      Xnew <- matrix(oGPO$par, q)
    }
    
    # Replicate?
    if(nrep > 1) Xnew <- Xnew[rep(1:nrow(Xnew), times = nrep),,drop = FALSE]
    
    Ynew <- apply(Xnew, 1, fun)
    
    X0 <- rbind(X0, Xnew)
    Y0 <- c(Y0, Ynew)
    
    # model update
    new.model <- try(modfun(X = X0, Z = Y0, covtype = model$covtype,
                            known = model$used_args$known, lower = model$used_args$lower, upper = model$used_args$upper,
                            noiseControl = model$used_args$noiseControl, settings = model$used_args$settings),
                     silent = TRUE)
    if(is(new.model,"try-error")){
      new.model <- try(update(model, Xnew = Xnew, Znew = Ynew, ginit = model$g * 1.01),
                       silent = TRUE)
    }
    
    if(is(new.model,"try-error")) break else model <- new.model
    if(control$trackbest){
      ps <- predict(model, model$X0)$mean
      Xbests <- rbind(Xbests, model$X0[which.min(ps),])
      Ybests <- c(Ybests, ps[which.min(ps)])
    }
    timings <- c(timings, difftime(Sys.time(), curtime, units = "sec"))
    curtime <- Sys.time()
    nuniq <- c(nuniq, nrow(model$X0))
    
  }
  return(list(par = X0, value = Y0, model = model, membest = Xbests, estbest = Ybests, timings = timings, nuniq = nuniq))
}


#' Extract non-dominated points from a set, or with respect to a reference Pareto front
#' @title Generic non-domination computation
#' @param points matrix (one point per row) from which to extract non-dominated points, or,
#' if a reference \code{ref} is provided, non-dominated points with respect to  \code{ref}
#' @param ref matrix (one point per row) of reference (faster if they are already Pareto optimal)
#' @param return.idx if \code{TRUE}, return indices instead of points
#' @return Non-dominated points from \code{points}, unless a \code{ref} is provided, in which case return points from \code{points} non-dominated by \code{ref}.
#' If \code{return.idx} is \code{TRUE}, only returns indices
#' @details Use Kung non-domination sorting
#' @noRd
#' @references
#' Kung, H. T., Luccio, F., & Preparata, F. P. (1975). On finding the maxima of a set of vectors. Journal of the ACM (JACM), 22(4), 469-476.
#' @examples
#' d <- 6
#' n <- 1000
#' n2 <- 1000
#'
#' test <- matrix(runif(d * n), n)
#' ref <- matrix(runif(d * n), n)
#' indPF <- nonDom(ref, return.idx = TRUE)
#' all(nonDom(ref) == ref[indPF,])
#'
#' system.time(res <- nonDom(test, ref[indPF,,drop = FALSE], return.idx = TRUE))
# '
# ' res2 <- rep(NA, n2)
# ' library(emoa)
# ' t0 <- Sys.time()
# ' for(i in 1:n2){
# '   res2[i] <- !is_dominated(t(rbind(test[i,, drop = FALSE], ref[indPF,])))[1]
# ' }
# ' print(Sys.time() - t0)
# '
# ' all(res == which(res2))
# '
# ' all(nonDom(test, ref) == test[res2,])
# '
nonDom <- function(points, ref = NULL, return.idx = FALSE){
  if(is.null(ref)){
    ordrs <- order(points[,1])
    if(return.idx) return(ordrs[nonDomInd_cpp(points[ordrs, , drop = FALSE])])
    return(points[ordrs[nonDomInd_cpp(points[ordrs, , drop = FALSE])],, drop = FALSE])
  }
  
  res <- nonDomSet(points, ref)
  if(return.idx) return(which(res))
  return(points[res, , drop = FALSE])
}


#' Computes variance reduction (s_n^2(x) - s_n+1^2(x))
#' @param model hetGP model
#' @param x new designs matrix
#' @param nr number of replicates at each x (default 1)
#' @examples
#' \dontrun{
#' x <- matrix(c(0.1, 0.2),1)
#' p0 <- predict(mod, x)$sd2
#' res1 <- PBBO:::vared(mod, x)
#'
#' # Compare with two alternatives
#' res2 <- hetGP:::update_predHGP(x = x, Xnew = x, model = mod)
#'
#' modup <- update(mod, Xnew = x, Znew = NA, maxit = 0)
#' res3 <- predict(modup, x)
#'
#' res1
#' p0 - res2$sd2
#' p0 - res3$sd2
#' }
#' @references
#' Gramacy, R. B. Surrogates: Gaussian Process Modeling, Design, and Optimization for the Applied Sciences
#' CRC Press, 2020
#' @noRd
vared <- function(model, x, nr = NULL){
  if(is.null(dim(x))) x <- matrix(x, 1)
  if(is.null(nr)) nr <- rep(1, nrow(x))
  
  kn <- cov_gen(x, model$X0, theta = model$theta, type = model$covtype)
  new_lambda <- predict(object = model, x = x, nugs.only = TRUE)$nugs/model$nu_hat
  vn <- drop(1 - kn %*% tcrossprod(model$Ki, kn)) + new_lambda / nr + model$eps
  gn <- - tcrossprod(model$Ki, kn) / vn
  res <- kn %*% gn %*% t(gn) %*% t(kn) * vn + 2 * kn %*% gn + 1 / vn
  res <- res * model$nu_hat
}

#' Hypervolume Sharpe ratio return and covariance
#' @param A matrix of assets, in R^p
#' @param l,u vector of lower and upper bounds in R^p
#' @return list with the return vector r and the covariance Q.
#' @references A. P. Guerreiro, C. M. Fonseca,
#' Hypervolume Sharpe-Ratio indicator: Formalization and first theoretical results,
#' International Conference on Parallel Problem Solving from Nature, 2016, 814-823.
#' @export
#' @examples
#' ################################################################################
#' ### 2 objectives example
#' ################################################################################
#' A <- matrix(runif(20*2),20)
#' res <- hyperSharperQ(A = A, l = c(0,0), u = c(1,1))
#' plot(A, pch = 20, xlim = c(0,1), ylim = c(0,1))
hyperSharperQ <- function(A, l, u){
  P <- hyperSharperP(A, l, u)
  r <- diag(P)
  Q <- P - tcrossprod(r)
  return(list(r = r, Q = Q))
}

#' Hypervolume Sharpe ratio maximization
#' @param A matrix of assets, in R^p
#' @param l,u vector of lower and upper bounds in R^p
#' @param eps jitter used in the inversion of the covariance matrix for numerical stability
#' @return list with the allocation vector a, corresponding Sharpe ratio value, return vector r and the covariance Q.
#' @references A. P. Guerreiro, C. M. Fonseca,
#' Hypervolume Sharpe-Ratio indicator: Formalization and first theoretical results,
#' International Conference on Parallel Problem Solving from Nature, 2016, 814-823.
#' @importFrom quadprog solve.QP
#' @export
#' @examples
#' ################################################################################
#' ### 2 objectives example
#' ################################################################################
#' set.seed(42)
#' nA <- 20 # Number of assets
#' p <- 2 # Number of objectives
#' A <- matrix(runif(nA * p), nA)
#' sol <- hyperSharpeMax(A = A, l = c(0, 0), u = c(1, 1))
#' plot(A, pch = 20, xlim = c(0, 1), ylim = c(0, 1))
#' points(A[which(sol$par > 1e-6),,drop = FALSE], col = 2)
hyperSharpeMax <- function(A, l, u, eps = sqrt(.Machine$double.eps)){
  tmp <- hyperSharperQ(A, l, u)
  res <- solve.QP(Dmat = tmp$Q + diag(eps, nrow(tmp$Q)), dvec = rep(0, nrow(A)), meq = 1,
                  Amat = cbind(tmp$r, diag(nrow(A))),
                  bvec = c(1, rep(0, nrow(A))))
  return(list(par = pmax(0,res$solution/sum(res$solution)), value = res$value,
              r = tmp$r, Q = tmp$Q))
}

#' Allocate replicates based on portfolio weights, with constraints on the number of possible designs
#' @title Allocation under maximum replication constraints
#' @param w vector of weights
#' @param q scalar batch size
#' @param mr vector of maximum number of evaluation of each element
#' @return allocation of integer number of runs depending on weights
#' @details proceeds by dichotomy
#' @export
#' @examples
#' set.seed(42)
#' n <- 10
#' w <- runif(n)^4
#' w <- w/sum(w)
#' q <- 50
#' mr <- c(rep(2, round(n/2)), rep(q, round(n/2)))
#' al <- allocq(w = w, q = q)
#' al_c <- allocq_c(w = w, q = q, mr = mr)
#' par(mfrow = c(1,2))
#' plot(w, pch = 20)
#' plot(mr, ylim = c(0, q), type = "b", col = 4, lty = 3)
#' segments(x0 = (1:n)-0.02, x1 = (1:n) - 0.02, y0 = rep(0, n), y1 = al)
#' segments(x0 = (1:n)+0.02, x1 = (1:n) + 0.02, y0 = rep(0, n), y1 = al_c, col = "red")
#' legend("topleft", legend = c("without max rep", "with max rep"), col = c(1, 2), lty= 1)
#' par(mfrow = c(1,1))
#' print(sum(al))
#' print(sum(al_c))
allocq_c <- function(w, q, mr){
  if(sum(mr) < q) stop("Not enough evaluations to allocate q \n")
  q0 <- q
  a <- 0
  w <- w + runif(n = length(w), min = 0, max = 1e-5) #to break ties
  w <- w/sum(w)
  b <- (1 / max(w)) * q + 1
  while(sum(round(pmin(w * b, mr))) < q){
    b <- b * 2
  }
  while(sum(round(pmin(q0 * w, mr))) != q){
    if(sum(round(pmin(q0 * w, mr))) < q) a <- (a + q0)/2 else b <- (b + q0)/2
    q0 <- (a + b)/2
  }
  return(round(pmin(q0 * w, mr)))
}

#' Allocate replicates based on portfolio weights
#' @title Allocation
#' @param w weights
#' @param q batch size
#' @return allocation of integer number of runs depending on weights
#' @details proceeds by dichotomy
#' @export
#' @examples
#' set.seed(42)
#' n <- 10
#' w <- runif(n)^4
#' w <- w/sum(w)
#' q <- 5
#' al <- allocq(w = w, q = q)
#' plot(w, pch = 20)
#' segments(x0 = 1:n, x1 = 1:n, y0 = rep(0, n), y1 = al/10)
#' # Asynchronous case
#' q2 <- q + 2
#' al2 <- allocq(w = w, q = q2)
#' plot(w, pch = 20, main = "q = 5, q' = 2")
#' for(i in 1:length(al)){
#'  if(al[i] > 0){
#'   for(j in 1:al[i]) arrows(x0 = i, x1 = i, y0 = (j-1)/10, y1 = j/10, code = 3, 
#'    angle = 90, lwd = 2, length = 0.1)
#'  }
#' }
#' for(i in 1:length(al2)){
#'  if(al2[i] > 0){
#'   for(j in 1:al2[i]) arrows(x0 = i + 0.05, x1 = i + 0.05, y0 = (j-1)/10, y1 = j/10,
#'    code = 3, angle = 90, col = "red", lwd = 2, lty = 2, length = 0.1)
#'  }
#' }
allocq <- function(w, q){
  q0 <- q
  a <- 0
  b <- (1/max(w)) * q + 1
  w <- w + runif(n = length(w), min = 0, max = 1e-5) #to break ties
  w <- w/sum(w)
  while(sum(round(q0*w)) != q){
    if(sum(round(q0*w)) < q) a <- (a + q0)/2 else b <- (b + q0)/2
    q0 <- (a + b)/2
  }
  return(round(q0*w))
}

#' Bayesian optimization loop with portfolio based batch EI
#' @title BO loop with massive batches
#' @param X0 initial design of experiments matrix
#' @param Y0 initial vector of responses at \code{X0}
#' @param model \code{homGP} or \code{hetGP} model
#' @param q batch size
#' @param fun test function to minimize
#' @param budget optimization budget
#' @param lower,upper domain bounds
#' @param control list with parameters 
#' \itemize{
#' \item nunif: number of uniformly sampled candidates points for acquisition function optimization
#' \item pop,gen: population size and number of generations for \code{\link[mco]{nsga2}}
#' \item maxit: passed to \code{\link[stats]{optim}}
#' \item trackbest: if \code{TRUE}, the best estimated solution will be saved at each iteration
#' \item maxtime: alternative stopping criterion
#' \item trace: boolean for printing messages
#' \item minPI: minimum probability of improvement for filtering solutions on the exploration-exploitation Pareto front
#' \item extendper: parameter to set the reference point for hypervolume computation
#' \item tol_dist: minimal distance to an existing design
#' \item tol_bdist: minimal distance to a design in the batch
#' \item max_rep: maximal degree of replication for a given design
#' \item max_brep maximal degree of replication per batch point
#' }
#' @importFrom mco nsga2
#' @importFrom stats cov2cor optim runif pnorm rnorm
#' @importFrom graphics axis filled.contour pairs
#' @return  A list with components:
#' \itemize{
#' \item \code{par}: all points evaluated,
#' \item \code{value}: the matrix of objective values at the points given in \code{par},
#' \item \code{model}: the last kriging models fitted.
#' \item \code{membest}: a matrix of best estimated designs at each iteration.
#' \item \code{estbest}: corresponding predicted mean values.
#' }
#' @references
#' Binois, M., Collier, N., Ozik, J. (2025). A portfolio approach to massively parallel Bayesian optimization. Journal of Artificial Intelligence Research, 82, pp. 137-167.
#' @export
#' @examples
#' d <- 2
#' n <- 10*d
#' N <- 5*n
#' budget <- 150 # Increase for better results
#' 
#' ## Noise field via standard deviation
#' noiseFun <- function(x){
#'   if(is.null(nrow(x)))
#'     x <- matrix(x, nrow = 1)
#'   return(1/5*(3*(2 + 2*sin(x[,1]*pi)*cos(x[,2]*3*pi) + 5*rowSums(x^2))))
#' }
#' 
#' ## Branin redefined in [0,1]^2
#' branin <- function(x){
#'  if(is.null(nrow(x))) x <- matrix(x, nrow = 1)
#'   x1 <- x[,1] * 15 - 5
#'   x2 <- x[,2] * 15
#'   return((x2 - 5/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(x1) + 10)
#'  }
#' 
#' ## data generating function combining mean and noise fields
#' ftest <- function(x){
#'   if(is.null(nrow(x))) x <- matrix(x, nrow = 1)
#'   return(branin(x) + rnorm(nrow(x), mean = 0, sd = noiseFun(x)))
#' }
#'
#'
#' ngrid <- 51
#' Xgrid <- as.matrix(expand.grid(seq(0,1,length.out = ngrid), seq(0,1,length.out = ngrid)))
#' Ygrid <- branin(Xgrid)
#' Ygrid <- cbind(Ygrid, noiseFun(Xgrid)^2)
#' x1 <- c(0.9616520, 0.15); x2 <- c(0.1238946, 0.8166644); x3 <- c(0.5427730, 0.15)
#' fstar <- 0.3978874
#'
#' X0 <- matrix(runif(n*d),n, d)
#' X0 <- X0[sample(1:n, size = N, replace = TRUE),]
#' Y0 <- ftest(X0)
#' 
#' mod <- mleHetGP(X0, Y0, covtype = "Matern5_2", known = list(beta0 = 0))
#' opt <- qHSRI_loop(X0, Y0, mod, q = 10, fun = ftest, budget = budget,
#'   lower = rep(0, d), upper = rep(1, d))
#' est <- predict(opt$model, opt$model$X0)$mean
#' xbest <- opt$model$X0[which.min(est),,drop=FALSE]
#' preds <- predict(opt$model, Xgrid)
#' 
#' par(mfrow = c(2, 2))
#' contour(matrix(preds$mean, ngrid), levels = c(1, seq(0,300,10)), 
#'  main="Mean prediction")
#' points(opt$model$X0, col = 4, pch = 20, cex = opt$model$mult/10)
#' contour(matrix(sqrt(preds$nugs), ngrid), nlevels = 21, 
#'  main="Variance prediction")
#' points(opt$model$X0, col = 4, pch = 20, cex = opt$model$mult/10)
#' contour(matrix(Ygrid[,1], ngrid), levels = c(1, seq(0,300,10)),
#'  main = "True function")
#' points(opt$model$X0, col = 4, pch = 20, cex = opt$model$mult/10)
#' points(rbind(t(x1), t(x2), t(x3)), pch=20, col="red")
#' points(xbest, col = "pink", pch = 15)
#' contour(matrix(sqrt(Ygrid[,2]^2), ngrid), nlevels = 21,
#'  main = "True variance")
#' points(rbind(t(x1), t(x2), t(x3)), pch=20, col="red")
#' points(opt$model$X0, col = 4, pch = 20, cex = opt$model$mult/10)
#' points(xbest, col = "pink", pch = 15)
#' par(mfrow = c(1, 1))
qHSRI_loop <- function(X0, Y0, model, q, fun, budget, lower, upper, control = NULL){
  d <- ncol(X0)
  if(is.null(control)) control <- list()
  if(is.null(control$nunif)) control$nunif <- 100*d # initial candidate points
  if(is.null(control$pop)) control$pop <- max(300, floor(50+2*sqrt(d))) # population size for nsga
  if(is.null(control$gen)) control$gen <- 100 # max iteration for nsga
  if(is.null(control$plot)) control$plot <- FALSE
  if(is.null(control$trace)) control$trace <- TRUE
  if(is.null(control$minPI)) control$minPI <- 1/3 # Minimum probability of improvement considered
  if(is.null(control$extendper)) control$extendper <- 0.2 # To define the refPoint as in GPareto
  if(is.null(control$trackbest)) control$trackbest <- TRUE
  if(is.null(control$tol_dist)) control$tol_dist <- 1e-2 # Minimal distance to an existing design
  if(is.null(control$tol_bdist)) control$tol_bdist <- 5e-2 # Minimal distance to a design in the batch
  if(is.null(control$max_rep)) control$max_rep <- 1e2 # Maximal degree of replication for a given design
  if(is.null(control$max_brep)) control$max_brep <- q # Maximal degree of replication per batch point
  if(is.null(control$maxtime)) control$maxtime <- Inf
  if(is.null(control$Gupta2018)) control$Gupta2018 <- FALSE
  if(control$Gupta2018 && control$minPI > 0) stop("The Gupta et al. method should use minPI = 0.")
  
  curtime <- Sys.time() # Inital time
  timings <- 0 # Store times
  nuniq <- NULL # Track number of unique designs over time
  modelfittimes <- model$time
  eval_times <- NULL
  
  if(is(model, "hetGP")) modfun <- mleHetGP else modfun <- mleHomGP
  
  if(control$trackbest){
    Xbests <- X0[which.min(Y0),,drop = F]
    Ybests <- min(Y0)
  }  else {
    Xbests <- Ybests <- NULL
  }
  
  Xs_old <- NULL # To store candidate points over iterations
  
  while(nrow(X0) < budget && timings[length(timings)] < control$maxtime){
    if(budget - nrow(X0) < q) q <- budget - nrow(X0) # To stay within budget limits
    Xs <- rbind(Xs_old,
                matrix(runif(control$nunif * d), control$nunif) %*% diag(upper - lower) + matrix(lower, control$nunif, d, byrow = TRUE))
    tmp <- predict(model, model$X0)$mean
    cst <- min(tmp)
    
    # Add few samples around current mean for exploitation
    localpert <- pmin(upper, pmax(lower, X0[which.min(tmp),] + rnorm(d*q, sd = 0.05)))
    localpert <- matrix(localpert, q, byrow = TRUE)
    Xs <- rbind(Xs, localpert)
    
    # Add existing designs to favor replication
    Xs <- unique(rbind(Xs, X0))
    
    # Add best EI point to alternatives
    EIs <- apply(Xs, 1, crit_EI, cst = cst, model = model)
    
    solEI <- optim(Xs[which.max(EIs),], fn = crit_EI, gr = deriv_crit_EI, model = model, cst = cst,
                   lower = lower, upper = upper, method = "L-BFGS-B", control = list(fnscale = -1))
    Xs <- rbind(Xs, solEI$par)
    
    # Find full PF
    fn <- function(x, model){
      ps <- predict(model, x)
      return(rbind(ps$mean, -sqrt(ps$sd2)))
    }
    
    solNS <- nsga2(fn = fn, idim = d, odim = 2, lower.bounds = lower, upper.bounds = upper, model = model,
                   popsize = control$pop, generations = control$gen, vectorized = TRUE)
    
    Xs <- rbind(Xs, solNS$par)
    rs <- t(fn(Xs, model = model))
    
    if(control$plot && nrow(X0) %% 25 == 0){
      if(d == 2){
        Xgrid <- as.matrix(expand.grid(seq(0,1,,51), seq(0,1,,51)))
        EIs <- apply(Xgrid, 1, crit_EI, model = model, cst = cst)
        filled.contour(matrix(EIs, 51), plot.axes = {axis(1);axis(2);
          points(solNS$par); points(model$X0, col = "blue", pch = 15)})
        Sys.sleep(1)
      }
      Xsp <- Xs
    }
    
    ids <- nonDom(rs, return.idx = T)
    Xs <- Xs[ids,, drop = FALSE]
    rs <- rs[ids,, drop = FALSE]
    
    # Keep some of the Xs for next round, but not too many
    if(nrow(Xs) > control$pop) Xs_old <- Xs[sample(1:nrow(Xs), size = control$pop),, drop = F] else Xs_old <- Xs
    
    # Filter point based on PI
    pps <- predict(model, Xs)
    Pis <- pnorm((cst-pps$mean)/sqrt(pps$sd2))
    ids2 <- which(Pis > control$minPI)
    
    if(length(ids2) < q && nrow(rs) >=q){
      if(control$trace) print("Not enough diversity (< q) in selected part of the PF, consider increasing pop or decrease minPI")
      Xs <- Xs[order(rs[,1])[1:q],,drop = FALSE]
      rs <- rs[order(rs[,1])[1:q],,drop = FALSE]
    }else{
      Xs <- Xs[ids2,,drop = FALSE]
      rs <- rs[ids2,,drop = FALSE]
    }
    
    # # Now compute variance reduction as additional objective
    vareds <- sqrt(apply(Xs, 1, vared, model = model))
    rs <- cbind(rs, -vareds)
    
    if(is.null(control$l)) l <- apply(rs, 2, min) else l <- control$l
    if(is.null(control$u)) u <- apply(rs, 2, max) else u <- control$u
    ul <- (u - l) # Extend both since there is the option to replicate
    if(is.null(control$u)) u <- u + control$extendper * ul
    if(is.null(control$l)) l <- l - control$extendper * ul
    
    if(control$Gupta2018){
      tmp <- rep(0, nrow(rs))
      tmp[sample(x = 1:nrow(rs), size = q, replace = FALSE)] <- 1/q
      res <- list(par = tmp)
      res$par <- res$par / sum(res$par)
    }else{
      res <- try(hyperSharpeMax(A = (rs - matrix(l, nrow = nrow(rs), ncol = ncol(rs), byrow = T)) %*% diag(1/(u-l)),
                                l = rep(0, length(l)), u = rep(1, length(l))))
    }
    if(is(res, "try-error")){
      print("Problem with Sharpe ratio")
      res <- list(par = c(rep(1/min(q, nrow(rs)), min(q, nrow(rs))),
                          rep(0, max(0, nrow(rs) - q))))
    }
    
    # Specific allocation for the noisy case
    if(control$Gupta2018){
      ids_new <- order(res$par, decreasing = TRUE)[1:(q-1)] # add EI solution
      Xnew <- rbind(Xs[ids_new,,drop = FALSE], solEI$par)
    }else{
      ids_new <- allocq_c(res$par, q - 1, mr = rep(control$max_brep, length(res$par)))
      ids_new[which.min(rs[,1])] <- ids_new[which.min(rs[,1])] + 1
      ids_new <- rep(1:length(ids_new), times = ids_new)
      Xnew <- Xs[ids_new,,drop = FALSE]
    }
    
    
    # Force replication if too close to existing design
    if(control$tol_dist > 0){
      for(i in 1:q) {
        dis <- sqrt(distance_cpp(Xnew[i,,drop = F], model$X0))
        mindis <- which.min(dis)
        if(dis[mindis] < control$tol_dist & model$mult[mindis] < control$max_rep)
          Xnew[i,] <- model$X0[mindis,]
      }
    }
    # Same between batch points
    if(control$tol_bdist > 0){
      for(i in 2:q){
        dis <- sqrt(distance_cpp(Xnew[i,,drop = F], Xnew[1:(i-1),,drop = F]))
        mindis <- which.min(dis)
        if(dis[mindis] < control$tol_bdist & dis[mindis] > 0 & sum(dis == 0) < control$max_brep)
          Xnew[i,] <- Xnew[mindis,]
      }
    }
    ttmp <- Sys.time()
    Ynew <- apply(Xnew, 1, fun)
    eval_times <- c(eval_times, difftime(Sys.time(), ttmp, units = "sec"))
    
    if(control$plot && nrow(X0) %% 25 == 0){
      Ysp <- t(fn(rbind(Xsp, model$X0, Xnew, solEI$par), model))
      Ysp <- cbind(Ysp, - sqrt(apply(rbind(Xsp, model$X0, Xnew, solEI$par), 1, vared, model = model)))
      Ysp <- rbind(Ysp, l, u)
      pps <- predict(model, rbind(Xsp, model$X0, Xnew, solEI$par))
      Pis <- pnorm((cst-pps$mean)/sqrt(pps$sd2))
      cols <- ifelse(Pis > control$minPI, "black", "gray")
      cols[(nrow(Ysp) - 1):nrow(Ysp)] <- "red"
      cols[(nrow(Xsp) + 1):(nrow(Xsp) + nrow(model$X0))] <- "blue"
      cols[nrow(Ysp) - 3] <- "cyan"
      cols[(nrow(Xsp) + nrow(model$X0) + 1):(nrow(Ysp) - 4)] <- "orange"
      pchs <- rep(1, nrow(Ysp))
      pchs[(nrow(Ysp) - 1):nrow(Ysp)] <- 20
      pchs[(nrow(Xsp) + 1):(nrow(Xsp) + nrow(model$X0))] <- 3
      pchs[nrow(Ysp) - 3] <- 17
      pchs[(nrow(Xsp) + nrow(model$X0) + 1):(nrow(Ysp) - 4)] <- 15
      pairs(Ysp, col = cols, pch = pchs)
    }
    
    
    X0 <- rbind(X0, Xnew)
    Y0 <- c(Y0, Ynew)
    
    # model update
    new.model <- try(modfun(X = X0, Z = Y0, covtype = model$covtype,
                            known = model$used_args$known, lower = model$used_args$lower, upper = model$used_args$upper,
                            noiseControl = model$used_args$noiseControl, settings = model$used_args$settings),
                     silent = TRUE)
    if(is(new.model,"try-error")){
      new.model <- try(update(model, Xnew = Xnew, Znew = Ynew, ginit = model$g * 1.01),
                       silent = TRUE)
    }
    modelfittimes <- c(modelfittimes, new.model$time)
    
    if(is(new.model,"try-error")) break else model <- new.model
    
    if(control$trackbest){
      ps <- predict(model, model$X0)$mean
      Xbests <- rbind(Xbests, model$X0[which.min(ps),])
      Ybests <- c(Ybests, ps[which.min(ps)])
    }
    timings <- c(timings, difftime(Sys.time(), curtime, units = "sec"))
    curtime <- Sys.time()
    nuniq <- c(nuniq, nrow(model$X0))
  }
  return(list(par = X0, value = Y0, model = model, membest = Xbests, estbest = Ybests, timings = timings, nuniq = nuniq,
              modelfittimes = modelfittimes, eval_times = eval_times))
  
}
