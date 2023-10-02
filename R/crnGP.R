## Model: noisy observations with unknown homoskedastic noise and common random numbers 
## K = nu^2 * (C + g * I) with C = Cx + Cs
# X0 unique designs matrix
# S0 corresponding seeds
# Z observations vector (all observations)
# theta vector of lengthscale hyperparameters (or one for isotropy)
# g noise variance for the process
# rho crn parameter
# beta0 trend
# stype: type of kronecker structure of the data
# @return loglikelihood value
logLikHomCRN <- function(X0, S0, Z, theta, g, rho, stype, beta0 = NULL, covtype = "Gaussian", eps = sqrt(.Machine$double.eps), env = NULL){
  n <- nrow(X0)
  d <- ncol(X0)
  
  # Temporarily store Cholesky transform of K in Ki
  Cx <- cov_gen(X1 = X0, theta = theta, type = covtype)
  if(is.null(env) || is.null(env$idS)){
    idS <- outer(S0, S0, "==")
    if(!is.null(env)) 
      env$idS <- idS
  }else{
    idS <- env$idS
  }
  Cs <- matrix(rho, n, n)
  Cs[idS] <- 1
  C <- Cx * Cs
  
  if(!is.null(env)){
    env$Cx <- Cx
    env$Cs <- Cs
    env$C <- C
  }
  
  Ki <- chol(C + diag(eps + g, n))
  ldetKi <- - 2 * sum(log(diag(Ki))) # log determinant from Cholesky
  Ki <- chol2inv(Ki)
  
  if(!is.null(env)) env$Ki <- Ki
  
  if(is.null(beta0))
    beta0 <- drop(colSums(Ki) %*% Z / sum(Ki))
  
  psi <- drop(crossprod(Z - beta0, Ki) %*% (Z - beta0))/n
  
  loglik <- -n/2 * log(2*pi) - n/2 * log(psi) + 1/2 * ldetKi - n/2
}


# derivative of log-likelihood for logLikHom with respect to theta with all observations (Gaussian kernel)
## Model: noisy observations with unknown homoskedastic noise
## K = nu^2 * (C + g * I)
# X0 unique designs matrix (last entry is the seed)
# Z observations vector (all observations)
# theta vector of lengthscale hyperparameters (or one for isotropy)
# g noise variance for the process
# rho crn parameter
# beta0 trend
# stype: type of kronecker structure of the data
## @return gradient with respect to theta and g
dlogLikHomCRN <- function(X0, S0, Z, theta, g, rho, stype, beta0 = NULL, covtype = "Gaussian",
                          eps = sqrt(.Machine$double.eps), components = c("theta", "g", "rho"), env = NULL){
  n <- nrow(X0)
  d <- ncol(X0)
  
  
  if(!is.null(env)){
    C <- env$C
    Cx <- env$Cx
    Cs <- env$Cs
    Ki <- env$Ki
  }else{
    Cx <- cov_gen(X1 = X0, theta = theta, type = covtype)
    idS <- outer(S0, S0, "==")
    Cs <- matrix(rho, n, n)
    Cs[idS] <- 1
    C <- Cx * Cs
    Ki <- chol2inv(chol(C + diag(eps + g, n)))
  }
  
  
  if(is.null(beta0))
    beta0 <- drop(colSums(Ki) %*% Z / sum(Ki))
  
  Z <- Z - beta0
  
  KiZ <- Ki %*% Z ## to avoid recomputing  
  
  psi <- drop(crossprod(Z, KiZ))/n
  
  tmp1 <- tmp2 <- tmp3 <- NULL
  
  # First component, derivative with respect to theta
  if("theta" %in% components){
    tmp1 <- rep(NA, length(theta))
    if(length(theta)==1){
      
      # Kfun <- function(X, theta, rho, i, j){
      #   d <- ncol(X)
      #   kx <- cov_gen(matrix(X[,-d,drop = F]), matrix(X[,-d, drop = F]), theta = theta)
      #   tmp <- outer(X[,d], X[,d], "==")
      #   ks <- matrix(rho, nrow(X), nrow(X))
      #   ks[tmp] <- 1
      #   res <- kx * ks
      #   return(res[i, j])
      # }
      # grad(Kfun, X = X0, x = theta, rho = rho, i = 1, j = 2)
      
      dC_dthetak <- partial_cov_gen(X1 = X0, theta = theta, type = covtype, arg = "theta_k") * C
      tmp1 <- 1/2 * crossprod(KiZ, dC_dthetak) %*% KiZ / psi - 1/2 * trace_sym(Ki, dC_dthetak)
    }else{
      for(i in 1:length(theta)){
        dC_dthetak <- partial_cov_gen(X1 = X0[,i, drop = F], theta = theta[i], type = covtype, arg = "theta_k") * C
        tmp1[i] <- 1/2 * crossprod(KiZ, dC_dthetak) %*% KiZ / psi - 1/2 * trace_sym(Ki, dC_dthetak)
      }
    } 
  }
  
  # Second component derivative with respect to g
  if("g" %in% components)
    tmp2 <- 1/2 * sum(KiZ^2) / psi - 1/2 * sum(diag(Ki))
  
  if("rho" %in% components){
    dC_drho <- matrix(1, n, n)
    ids <- outer(S0, S0, "==")
    dC_drho[ids] <- 0
    dC_drho <- dC_drho * Cx
    tmp3 <- 1/2 * crossprod(KiZ, dC_drho) %*% KiZ / psi - 1/2 * trace_sym(Ki, dC_drho)
    
  }
  
  return(c(tmp1, tmp2, tmp3))
}


#' Gaussian process regression when seed (or trajectory) information is provided, based on maximum likelihood estimation of the 
#' hyperparameters. Trajectory handling involves observing all times for any given seed.
#' @title Gaussian process modeling with correlated noise
#' @param X matrix of all designs, one per row. The last column is assumed to contain the integer seed value.
#' @param T0 optional vector of times (same for all \code{X}s)
#' @param Z vector of all observations. If \code{ts} is provided, the \code{Z} is a matrix of size \code{nrow(X) x length(ts)}.
#' @param stype structural assumptions, options include:
#' \itemize{
#' \item \code{none}: no structure, regular matrix inversion is used (only when no time is present);
# ' \item \code{XS}: a Kronecker structure is used between X and S.
# ' It becomes interesting if many seeds have been evaluated at different x locations (and mostly the same seeds).
#' }
#' When time is present, the Kronecker structure is always used (the alternative is to provide times as an extra variable in \code{X})
#' Using the Kronecker structure becomes efficient when the product (nx x ns) x nt becomes large.
#' @param lower,upper optional bounds for the \code{theta} parameter (see \code{\link[hetGP]{cov_gen}} for the exact parameterization).
#' In the multivariate case, it is possible to give vectors for bounds (resp. scalars) for anisotropy (resp. isotropy) 
#' @param known optional list of known parameters, e.g., \code{beta0}, \code{theta}, \code{g} or \code{rho}.
#' @param covtype covariance kernel type, either 'Gaussian', 'Matern5_2' or 'Matern3_2', see \code{\link[hetGP]{cov_gen}}
#' @param noiseControl list with element, 
#' \itemize{
#' \item \code{g_bounds}, vector providing minimal and maximal noise to signal ratio;
#' \item \code{rho_bounds}, vector providing minimal and maximal correlation between seed values, in [0,1];
#' } 
#' @param init optional list specifying starting values for MLE optimization, with elements:
#' \itemize{
#'  \item \code{theta_init} initial value of the theta parameters to be optimized over (default to 10\% of the range determined with \code{lower} and \code{upper})
#'  \item \code{g_init} initial value of the nugget parameter to be optimized over.
#'  \item \code{rho_init} initial value of the seed correlation parameter.
#' }
#' @param maxit maximum number of iteration for L-BFGS-B of \code{\link[stats]{optim}}
#' @param eps jitter used in the inversion of the covariance matrix for numerical stability
#' @param settings list with argument \code{return.Ki}, to include the inverse covariance matrix in the object for further use (e.g., prediction).
#' Arguments \code{factr} (default to 1e9) and \code{pgtol} are available to be passed to \code{control} for L-BFGS-B in \code{\link[stats]{optim}}. 
#' @return a list which is given the S3 class "\code{CRNGP}", with elements:
#' \itemize{
#' \item \code{theta}: maximum likelihood estimate of the lengthscale parameter(s),
#' \item \code{g}: maximum likelihood estimate of the nugget variance,
#' \item \code{rho}: maximum likelihood estimate of the seed correlation parameter,
#' \item \code{trendtype}: either "\code{SK}" if \code{beta0} is given, else "\code{OK}" 
#' \item \code{beta0}: estimated trend unless given in input,
#' \item \code{nu_hat}: plugin estimator of the variance,
#' \item \code{ll}: log-likelihood value,
#' \item \code{X0}, \code{S0}, \code{T0}: values for the spatial, seed and time designs 
#' \item \code{Z}, \code{eps}, \code{covtype}, \code{stype},: values given in input,
#' \item \code{call}: user call of the function
#' \item \code{used_args}: list with arguments provided in the call
#' \item \code{nit_opt}, \code{msg}: \code{counts} and \code{msg} returned by \code{\link[stats]{optim}}
#' \item \code{Ki}: inverse covariance matrix (not scaled by \code{nu_hat}) (if \code{return.Ki} is \code{TRUE} in \code{settings})
#' \item \code{Ct}: if time is used, corresponding covariance matrix.
#' \item \code{time}: time to train the model, in seconds.
#' 
#'}
#' @details
#' The global covariance matrix of the model is parameterized as \code{nu_hat * (Cx + g Id) * Cs = nu_hat * K},
#' with \code{Cx} the spatial correlation matrix between unique designs, depending on the family of kernel used (see \code{\link[hetGP]{cov_gen}} for available choices) and values of lengthscale parameters.
#' \code{Cs} is the correlation matrix between seed values, equal to 1 if the seeds are equal, \code{rho} otherwise.
#' \code{nu_hat} is the plugin estimator of the variance of the process.
#' 
#' Compared to \code{\link[hetGP]{mleHomGP}}, here the replications have a specific identifier, i.e., the seed.
#' @note This function is experimental at this time and could evolve in the future.
#' @seealso \code{\link[hetGP]{predict.CRNGP}} for predictions, \code{\link[hetGP]{simul.CRNGP}} for generating conditional simulation on a Kronecker grid.
#' \code{summary} and \code{plot} functions are available as well. 
#' @references 
#' Xi Chen, Bruce E Ankenman, and Barry L Nelson. The effects of common random numbers on stochastic kriging metamodels. ACM Transactions on Modeling and Computer Simulation (TOMACS), 22(2):1-20, 2012.\cr \cr
#' 
#' Michael Pearce, Matthias Poloczek, and Juergen Branke. Bayesian simulation optimization with common random numbers. In 2019 Winter Simulation Conference (WSC), pages 3492-3503. IEEE, 2019. \cr \cr
#' 
#' A Fadikar, M Binois, N Collier, A Stevens, KB Toh, J Ozik. Trajectory-oriented optimization of stochastic epidemiological models. arXiv preprint arXiv:2305.03926
#' @export
#' @examples
#' ##------------------------------------------------------------
#' ## Example 1: CRN GP modeling on 1d sims
#' ##------------------------------------------------------------
#' #' set.seed(42)
#' nx <- 50
#' ns <- 5
#' x <- matrix(seq(0,1, length.out = nx), nx)
#' s <- matrix(seq(1, ns, length.out = ns))
#' g <- 1e-3
#' theta <- 0.01
#' KX <- cov_gen(x, theta = theta)
#' rho <- 0.3
#' KS <- matrix(rho, ns, ns)
#' diag(KS) <- 1
#' YY <- MASS::mvrnorm(n = 1, mu = rep(0, nx*ns), Sigma = kronecker(KX, KS) + g * diag(nx*ns))
#' YYmat <- matrix(YY, ns, nx)
#' matplot(x, t(YYmat), pch = 1, type = "b", lty = 3)
#' 
#' Xgrid <- as.matrix(expand.grid(s, x))
#' Xgrid <- cbind(Xgrid[,2], Xgrid[,1])
#' ids <- sample(1:nrow(Xgrid), 20)
#' X0 <- Xgrid[ids,]
#' Y0 <-  YY[ids]
#' points(X0[,1], Y0, pch = 20, col = 1 + ((X0[,2] - 1) %% 6))
#' 
#' model <- mleCRNGP(X0, Y0, known = list(theta = 0.01, g = 1e-3, rho = 0.3))
#' 
#' preds <- predict(model, x = Xgrid, xprime = Xgrid)
#' matlines(x, t(matrix(preds$mean, ns, nx)), lty = 1)

#' # prediction on new seed (i.e., average prediction)
#' xs1 <- cbind(x, ns+1)
#' predsm <- predict(model, x = xs1)
#' lines(x, predsm$mean, col = "orange", lwd = 3)
#' lines(x, predsm$mean + 2 * sqrt(predsm$sd2), col = "orange", lwd = 2, lty = 3)
#' lines(x, predsm$mean - 2 * sqrt(predsm$sd2), col = "orange", lwd = 2, lty = 3)
#' 
#' # Conditional realizations
#' sims <- MASS::mvrnorm(n = 1, mu = preds$mean, Sigma = 1/2 * (preds$cov + t(preds$cov)))
#' plot(Xgrid[,1], sims, col = 1 + ((Xgrid[,2] - 1) %% 6))
#' points(X0[,1], Y0, pch = 20, col = 1 + ((X0[,2] - 1) %% 6))
#' \dontrun{
#' ##------------------------------------------------------------
#' ## Example 2: Homoskedastic GP modeling on 2d sims
#' ##------------------------------------------------------------
#' set.seed(2)
#' nx <- 31
#' ns <- 5
#' d <- 2
#' x <- as.matrix(expand.grid(seq(0,1, length.out = nx), seq(0,1, length.out = nx)))
#' s <- matrix(seq(1, ns, length.out = ns))
#' Xgrid <- as.matrix(expand.grid(seq(1, ns, length.out = ns), seq(0,1, length.out = nx), 
#'                                seq(0,1, length.out = nx)))
#' Xgrid <- Xgrid[,c(2, 3, 1)]
#' g <- 1e-3
#' theta <- c(0.02, 0.05)
#' KX <- cov_gen(x, theta = theta)
#' rho <- 0.33
#' KS <- matrix(rho, ns, ns)
#' diag(KS) <- 1
#' YY <- MASS::mvrnorm(n = 1, mu = rep(0, nx*nx*ns), Sigma = kronecker(KX, KS) + g * diag(nx*nx*ns))
#' YYmat <- matrix(YY, ns, nx*nx)
#' filled.contour(matrix(YYmat[1,], nx))
#' filled.contour(matrix(YYmat[2,], nx))
#' 
#' ids <- sample(1:nrow(Xgrid), 80)
#' X0 <- Xgrid[ids,]
#' Y0 <-  YY[ids]
#' 
#' ## Uncomment below for For 3D visualisation
#' # library(rgl)
#' # plot3d(Xgrid[,1], Xgrid[,2], YY, col = 1 + (Xgrid[,3] - 1) %% 6)
#' # points3d(X0[,1], X0[,2], Y0, size = 10, col = 1 + ((X0[,3] - 1) %% 6))
#' 
#' model <- mleCRNGP(X0, Y0, know = list(beta0 = 0))
#' 
#' preds <- predict(model, x = Xgrid, xprime = Xgrid)
#' # surface3d(unique(Xgrid[1:nx^2,1]),unique(Xgrid[,2]), matrix(YY[Xgrid[,3]==1], nx), 
#' #   front = "lines", back = "lines")
#' # aspect3d(1, 1, 1)
#' # surface3d(unique(Xgrid[1:nx^2,1]),unique(Xgrid[,2]), matrix(preds$mean[Xgrid[,3]==1], nx), 
#' #   front = "lines", back = "lines", col = "red")
#' plot(preds$mean, YY)
#' 
#' # prediction on new seed (i.e., average prediction)
#' xs1 <- cbind(x, ns+1)
#' predsm <- predict(model, x = xs1)
#' # surface3d(unique(x[,1]), unique(x[,2]), matrix(predsm$mean, nx), col = "orange", 
#' #   front = "lines", back = "lines")
#' 
#' # Conditional realizations
#' sims <- MASS::mvrnorm(n = 1, mu = preds$mean, Sigma = 1/2 * (preds$cov + t(preds$cov)))
#' # plot3d(X0[,1], X0[,2], Y0, size = 10, col = 1 + ((X0[,3] - 1) %% 6))
#' # surface3d(unique(x[,1]), unique(x[,2]), matrix(sims[Xgrid[,3] == 1], nx), col = 1, 
#' #   front = "lines", back = "lines")
#' # surface3d(unique(x[,1]), unique(x[,2]), matrix(sims[Xgrid[,3] == 2], nx), col = 2, 
#' #   front = "lines", back = "lines")
#'
#' # Faster alternative for conditional realizations 
#' # (note: here the design points are part of the simulation points)
#' Xgrid0 <- unique(Xgrid[, -(d + 1), drop = FALSE])
#' sims2 <- simul(object = model,Xgrid = Xgrid, ids = ids, nsim = 5, check = TRUE) 
#' 
#' ##------------------------------------------------------------
#' ## Example 3: Homoskedastic GP modeling on 1d trajectories (with time)
#' ##------------------------------------------------------------
#' set.seed(42)
#' nx <- 11
#' nt <- 9
#' ns <- 7
#' x <- matrix(sort(seq(0,1, length.out = nx)), nx)
#' s <- matrix(sort(seq(1, ns, length.out = ns)))
#' t <- matrix(sort(seq(0, 1, length.out = nt)), nt)
#' covtype <- "Matern5_2"
#' g <- 1e-3
#' theta <- c(0.3, 0.5)
#' KX <- cov_gen(x, theta = theta[1], type = covtype)
#' KT <- cov_gen(t, theta = theta[2], type = covtype)
#' rho <- 0.3
#' KS <- matrix(rho, ns, ns)
#' diag(KS) <- 1
#' XST <- as.matrix(expand.grid(x, s, t))
#' 
#' Kmc <- kronecker(chol(KT), kronecker(chol(KS), chol(KX)))
#' YY <- t(Kmc) %*% rnorm(nrow(Kmc))
#' 
#' ninit <- 50
#' XS <- as.matrix(expand.grid(x, s))
#' ids <- sort(sample(1:nrow(XS), ninit))
#' XST0 <- cbind(XS[ids[rep(1:ninit, each = nt)],], rep(t[,1], times = ninit))
#' X0 <- XST[which(duplicated(rbind(XST, XST0), fromLast = TRUE)),]
#' Y0 <-  YY[which(duplicated(rbind(XST, XST0), fromLast = TRUE))]
#' 
#' # tmp <- hetGP:::find_reps(X = X0[,-3], Y0)
#' model <- mleCRNGP(X = XS[ids,], T0=t, Z = matrix(Y0, ncol = nt), covtype = covtype)
#' 
#' preds <- predict(model, x = XS, xprime = XS)
#' 
#' # compare with regular CRN GP
#' mref <- mleCRNGP(X = X0[, c(1, 3, 2)], Z = Y0, covtype = covtype)
#' pref <- predict(mref, x = XST[, c(1, 3, 2)], xprime = XST[, c(1, 3, 2)])
#' 
#' print(model$time) # Use Kronecker structure for time
#' print(mref$time)
#' 
#' plot(as.vector(preds$mean), YY)
#' plot(pref$mean, YY) 
#' 
#' }
mleCRNGP <- function(X, Z, T0 = NULL, stype = c("none", "XS"), lower = NULL, upper = NULL, known = NULL,
                     noiseControl = list(g_bounds = c(sqrt(.Machine$double.eps)*10, 1e2),
                                         rho_bounds = c(0.001, 0.9)),
                     init = NULL,
                     covtype = c("Gaussian", "Matern5_2", "Matern3_2"),
                     maxit = 100, eps = sqrt(.Machine$double.eps), settings = list(return.Ki = TRUE, factr = 1e7)){
  
  if(is.null(dim(X))) X <- matrix(X, ncol = 1)
  if(is.null(T0) && nrow(X) != length(Z)) stop("Dimension mismatch between Z and X")
  if(!is.null(T0) && nrow(X) != nrow(Z)) stop("Dimension mismatch between Z and X")
  stype <- match.arg(stype)
  
  # elem <- find_reps(X, Z, return.Zlist = F)
  X0 <- X[,-ncol(X), drop = F] #elem$X0
  # Z <- elem$Z
  S0 <- X[,ncol(X)]
  if(!is.null(T0) && is.null(dim(T0))) T0 <- matrix(T0, ncol = 1)
  
  d <- ncol(X0) # spatial dimension
  if(sum(abs(S0 - as.integer(S0))) > 0) 
    stop("The last column of X0 is assumed to contain integer values seed information.")
  
  covtype <- match.arg(covtype)
  
  if(is.null(lower) || is.null(upper)){
    auto_thetas <- auto_bounds(X = X0, covtype = covtype)
    if(is.null(lower)) lower <- auto_thetas$lower
    if(is.null(upper)) upper <- auto_thetas$upper
    if(is.null(known[["theta"]]) && is.null(init$theta)) init$theta <- sqrt(upper * lower)
  }
  if(!is.null(T0)){
    auto_thetasT <- auto_bounds(X = T0, covtype = covtype)
    lower <- c(lower, auto_thetasT$lower)
    upper <- c(upper, auto_thetasT$upper)
    if(length(init$theta) < length(lower)) init$theta <- sqrt(upper * lower)
  }
  
  if(length(lower) != length(upper)) stop("upper and lower should have the same size")
  
  ## Save time to train model
  tic <- proc.time()[3]
  
  if(is.null(settings$return.Ki)) settings$return.Ki <- TRUE
  if(is.null(noiseControl$g_bounds)) noiseControl$g_bounds <- c(sqrt(.Machine$double.eps), 1e2)
  if(is.null(noiseControl$rho_bounds)) noiseControl$rho_bounds <- c(0, 0.9)
  
  
  g_min <- noiseControl$g_bounds[1]
  g_max <- noiseControl$g_bounds[2]
  
  beta0 <- known$beta0
  
  n <- nrow(X0)
  
  if(is.null(n))
    stop("X0 should be a matrix. \n")
  
  if(is.null(known[["theta"]]) && is.null(init$theta)) 
    init$theta <- 0.9 * lower + 0.1 * upper # useful for mleHetGP
  if(is.null(known$g) && is.null(init$g)){
    init$g <- 0.1
  }
  if(is.null(known$rho) && is.null(init$rho)) init$rho <- 0.1
  
  trendtype <- 'OK'
  if(!is.null(beta0))
    trendtype <- 'SK'
  
  ## General definition of fn and gr
  fn <- function(par, X0, S0, T0, Z, beta0, theta, g, rho, env){
    idx <- 1 # to store the first non used element of par
    
    if(is.null(theta)){
      theta <- par[1:length(init$theta)]
      idx <- idx + length(init$theta)
    }
    
    if(is.null(g)){
      g <- par[idx]
      idx <- idx + 1
    }
    
    if(is.null(rho)){
      rho <- par[idx]
    }
    
    if(is.null(T0)){
      loglik <- logLikHomCRN(X0 = X0, S0 = S0, Z = Z, theta = theta, g = g, rho = rho, 
                             stype = stype, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    }else{
      loglik <- logLikHomCRNT(X0 = X0, S0 = S0, T0 = T0, Z = Z, theta = theta, g = g, rho = rho,
                              beta0 = beta0, covtype = covtype, eps = eps, env = env)
    }
    
    
    # library(numDeriv)
    # grad(logLikHomCRN, x = theta, X0 = X0, Z0 = Z0, Z = Z, mult = mult, g = g, rho = rho, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    # grad(logLikHomCRN, x = g, X0 = X0, Z0 = Z0, Z = Z, mult = mult, theta = theta, rho = rho, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    # grad(logLikHomCRN, x = rho, X0 = X0, Z0 = Z0, Z = Z, mult = mult, theta = theta, g = g, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    # dlogLikHomCRN(X0 = X0, Z0 = Z0, Z = Z, mult = mult, theta = theta, g = g, rho = rho, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    
    # grad(logLikHomCRNT, x = theta, X0 = X0, S0 = S0, T0 = T0, Z = Z, g = g, rho = rho, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    # grad(logLikHomCRNT, x = g, X0 = X0, S0 = S0, T0 = T0, Z = Z, theta = theta, rho = rho, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    # grad(logLikHomCRNT, x = rho, X0 = X0, S0 = S0, T0 = T0, Z = Z, theta = theta, g = g, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    # dlogLikHomCRNT(X0 = X0, S0 = S0, T0 = T0, Z = Z, theta = theta, g = g, rho = rho, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    
    if(!is.null(env) && !is.na(loglik)){
      if(is.null(env$max_loglik) || loglik > env$max_loglik){
        env$max_loglik <- loglik
        env$arg_max <- par
      }
    } 
    
    return(loglik)
  }
  
  gr <- function(par, X0, S0, T0, Z, beta0, theta, g, rho, env){
    idx <- 1
    components <- NULL
    
    if(is.null(theta)){
      theta <- par[1:length(init$theta)]
      idx <- idx + length(init$theta)
      components <- "theta"
    }
    
    if(is.null(g)){
      g <- par[idx]
      idx <- idx + 1
      components <- c(components, "g")
    }
    
    if(is.null(rho)){
      rho <- par[idx]
      components <- c(components, "rho")
    }
    if(is.null(T0)){
      return(dlogLikHomCRN(X0 = X0, S0 = S0, Z = Z, theta = theta, g = g, rho = rho, stype = stype,
                           beta0 = beta0, covtype = covtype, eps = eps,
                           components = components, env = env))
    }else{
      return(dlogLikHomCRNT(X0 = X0, T0 = T0, S0 = S0, Z = Z, theta = theta, g = g, rho = rho, beta0 = beta0, covtype = covtype, eps = eps,
                            components = components, env = env))
    }
    
  }
  
  ## All known
  envtmp <- environment()
  
  if(!is.null(known$g) && !is.null(known[["theta"]]) && !is.null(known$rho)){
    theta_out <- known[["theta"]]
    g_out <- known$g
    rho_out <- known$rho
    if(is.null(T0)){
      out <- list(value = logLikHomCRN(X0 = X0, S0 = S0, Z = Z, theta = theta_out, g = g_out, 
                                       rho = rho_out, stype = stype, beta0 = beta0, covtype = covtype, eps = eps, env = envtmp),
                  message = "All hyperparameters given", counts = 0, time = proc.time()[3] - tic)
    }else{
      out <- list(value = logLikHomCRNT(X0 = X0, S0 = S0, T0 = T0, Z = Z, theta = theta_out, 
                                        g = g_out, rho = rho_out, stype = stype, beta0 = beta0, covtype = covtype, eps = eps, env = envtmp),
                  message = "All hyperparameters given", counts = 0, time = proc.time()[3] - tic)
    }
    
  }else{
    parinit <- lowerOpt <- upperOpt <- NULL
    if(is.null(known[["theta"]])){
      parinit <- init$theta
      lowerOpt <- c(lower)
      upperOpt <- c(upper)
    }
    if(is.null(known$g)){
      parinit <- c(parinit, init$g)
      lowerOpt <- c(lowerOpt, g_min)
      upperOpt <- c(upperOpt, g_max)
    }
    if(is.null(known$rho)){
      parinit <- c(parinit, init$rho)
      lowerOpt <- c(lowerOpt, noiseControl$rho_bounds[1])
      upperOpt <- c(upperOpt, noiseControl$rho_bounds[2])
    }
    
    out <- try(optim(par = parinit, fn = fn, gr = gr, method = "L-BFGS-B", lower = lowerOpt, upper = upperOpt, 
                     theta = known[["theta"]], g = known$g, rho = known$rho,
                     X0 = X0, S0 = S0, T0 = T0, Z = Z, beta0 = beta0,
                     control = list(fnscale = -1, maxit = maxit, factr = settings$factr, pgtol = settings$pgtol), env = envtmp))
    ## Catch errors when at least one likelihood evaluation worked
    if(is(out, "try-error"))
      out <- list(par = envtmp$arg_max, value = envtmp$max_loglik, counts = NA,
                  message = "Optimization stopped due to NAs, use best value so far")
    
    idx <- 1
    if(is.null(known[["theta"]])){
      theta_out <- out$par[1:length(init$theta)]
      idx <- idx + length(theta_out)
    } else  theta_out <- known[["theta"]]
    
    
    if(is.null(known$g)){
      g_out <- out$par[idx]
    } else g_out <- known$g
    
    if(is.null(known$rho)) rho_out <- out$par[length(out$par)] else rho_out <- known$rho
    
  }
  
  Cx <- cov_gen(X1 = X0, theta = theta_out[1:d], type = covtype)
  idS <- envtmp$idS #outer(S0, S0, "==")
  Cs <- matrix(rho_out, n, n)
  Cs[idS] <- 1
  
  if(!is.null(T0)){
    Ct <- cov_gen(X1 = T0, theta = theta_out[d+1], type = covtype)
    SCxs <- svd(Cx*Cs) # C
    SCt <- svd(Ct) # R
    Ki <- tcrossprod(kronecker(SCt$u, SCxs$u) * rep(1/(kronecker(SCt$d, SCxs$d) + g_out), each = length(Z)), kronecker(SCt$u, SCxs$u))
  }else{
    Ct <- NULL
    C <- Cx * Cs
    Ki <- chol2inv(chol(C + diag(eps + g_out, n)))
  }
  
  if(is.null(beta0))
    beta0 <- drop(colSums(Ki) %*% as.vector(Z) / sum(Ki))
  
  nu <- drop(crossprod(as.vector(Z) - beta0, Ki) %*% (as.vector(Z) - beta0)) / length(Z)
  
  
  res <- list(theta = theta_out, g = g_out, rho = rho_out, nu_hat = as.numeric(nu), ll = out$value, nit_opt = out$counts,
              beta0 = beta0, trendtype = trendtype, covtype = covtype, stype = stype, msg = out$message, eps = eps,
              X0 = X0, S0 = S0, T0 = T0, Z = Z, call = match.call(), Ct = Ct, 
              used_args = list(lower = lower, upper = upper, known = known, noiseControl = noiseControl),
              time = proc.time()[3] - tic) # hess = hess)
  
  if(settings$return.Ki) res <- c(res, list(Ki = Ki))
  
  class(res) <- "CRNGP"
  return(res)
}


#' @method summary CRNGP
#' @export
summary.CRNGP <- function(object,...){
  ans <- object
  class(ans) <- "summary.CRNGP"
  ans
}

#' @export
print.summary.CRNGP <- function(x, ...){
  
  if(is.null(x$T0)){
    cat("N = ", length(x$Z), " n = ", length(x$Z), " d = ", ncol(x$X0), "\n")
  }else{
    cat("N = ", length(x$Z), " nxs = ", nrow(x$X0), " nt =  ", length(x$T0), "d = ", ncol(x$X0), "\n")
  }
  cat(x$covtype, " covariance lengthscale values: ", x$theta, "\n")
  
  cat("Homoskedastic nugget value: ", x$g, "\n")
  
  cat("CRN Correlation: ", x$rho, "\n")
  
  cat("Variance/scale hyperparameter: ", x$nu_hat, "\n")
  
  if(x$trendtype == "SK"){
    cat("Given constant trend value: ", x$beta0, "\n")
  }else{
    cat("Estimated constant trend value: ", x$beta0, "\n")
  }
  
  cat("MLE optimization: \n", "Log-likelihood = ", x$ll, "; Nb of evaluations (obj, gradient) by L-BFGS-B: ", x$nit_opt, "; message: ", x$msg, "\n")
}

#' @method print CRNGP
#' @export
print.CRNGP <- function(x, ...){
  print(summary(x))
}



#' Gaussian process predictions using a GP object for correlated noise (of class \code{CRNGP})
#' @param x matrix of designs locations to predict at (one point per row). Last column is for the integer valued seed. 
#' If trajectories are considered, i.e., with time, the prediction will occur at the same times as the training data unless \code{t0} is provided.
#' @param object an object of class \code{CRNGP}; e.g., as returned by \code{\link[hetGP]{mleCRNGP}}
#' @param xprime optional second matrix of predictive locations to obtain the predictive covariance matrix between \code{x} and \code{xprime}
#' @param t0 single column matrix of times to predict at, if trajectories are considered. By default the prediction is at the same times as the training data.
#' @param ... no other argument for this method
#' @return list with elements
#' \itemize{
#' \item \code{mean}: kriging mean;
#' \item \code{sd2}: kriging variance (filtered, e.g. without the nugget value)
#' \item \code{cov}: predictive covariance matrix between \code{x} and \code{xprime}
#' \item \code{nugs}: nugget value at each prediction location, for consistency with \code{\link[hetGP]{mleHomGP}}.
#' }
#' @importFrom MASS ginv
#' @details The full predictive variance corresponds to the sum of \code{sd2} and \code{nugs}. See \code{\link[hetGP]{mleHomGP}} for examples.
#' @method predict CRNGP
#' @export
predict.CRNGP <- function(object, x, xprime = NULL, t0 = NULL, ...){
  if(is.null(dim(x))){
    x <- matrix(x, nrow = 1)
    if(ncol(x) - 1 != ncol(object$X0)) stop("Problem with x format")
  }
  s <- x[, ncol(x), drop = F]
  x <- x[, -ncol(x), drop = F]
  
  if(!is.null(xprime)){
    if(is.null(dim(xprime)))
       xprime <- matrix(xprime, nrow = 1)
    if(ncol(xprime) - 1 != ncol(object$X0)) stop("Problem with xprime format")
    sp <- xprime[, ncol(xprime), drop = F]
    xprime <- xprime[, -ncol(xprime), drop = F]
  }
  
  d <- ncol(x)
  
  if(is.null(object$Ki)){
    warning("GP model does not have a Ki component - prediction is slower.")
    Cx <- cov_gen(X1 = object$X0, theta = object$theta[1:d], type = object$covtype)
    idS <- outer(object$S0, object$S0, "==")
    Cs <- matrix(object$rho, nrow(object$X0), nrow(object$X0))
    Cs[idS] <- 1
    if(is.null(object$T0)){
      C <- Cx * Cs
      Ki <- chol2inv(chol(C + diag(object$eps + object$g, nrow(object$X0))))
    }else{
      Ct <- cov_gen(X1 = object$T0, theta = object$theta[d+1], type = object$covtype)
      SCxs <- svd(Cx*Cs) # C
      SCt <- svd(Ct) # R
      Ki <- tcrossprod(kronecker(SCt$u, SCxs$u) * rep(1/(kronecker(SCt$d, SCxs$d) + object$g), each = length(object$Z)), 
                       kronecker(SCt$u, SCxs$u))
    }
  }
  
  object$Ki <- object$Ki / object$nu_hat
  
  kx <- object$nu_hat * cov_gen(X1 = x, X2 = object$X0, theta = object$theta[1:d], type = object$covtype)
  tmp <- outer(s, object$S0, "==")
  ks <- matrix(object$rho, nrow(x), nrow(object$X0))
  ks[tmp] <- 1
  kx <- kx * ks
  
  if(!is.null(object$T0)){
    if(is.null(t0))
      kt <- object$Ct
    else
      kt <- cov_gen(X1 = t0, X2 = object$T0, theta = object$theta[d+1], type = object$covtype)
    kx <- kronecker(kt, kx)
  }
  
  nugs <- rep(object$nu_hat * object$g, nrow(x))
  mean <- as.vector(object$beta0 + kx %*% (object$Ki %*% (as.vector(object$Z) - object$beta0)))
  
  if(object$trendtype == 'SK'){
    sd2 <- as.vector(object$nu_hat - fast_diag(kx, tcrossprod(object$Ki, kx)))
  }else{
    sd2 <- as.vector(object$nu_hat - fast_diag(kx, tcrossprod(object$Ki, kx)) + (1 - tcrossprod(rowSums(object$Ki), kx))^2/sum(object$Ki))
  }
  
  ## In case of numerical errors, some sd2 values may become negative
  if(any(sd2 < 0)){
    # object$Ki <- ginv(add_diag(cov_gen(X1 = object$X0, theta = object$theta, type = object$covtype), object$g/object$mult + object$eps))/object$nu_hat
    # mean <- as.vector(object$beta0 + kx %*% (object$Ki %*% (object$Z0 - object$beta0)))
    # 
    # if(object$trendtype == 'SK'){
    #   sd2 <- as.vector(object$nu_hat - fast_diag(kx, tcrossprod(object$Ki, kx)))
    # }else{
    #   sd2 <- as.vector(object$nu_hat - fast_diag(kx, tcrossprod(object$Ki, kx)) + (1 - tcrossprod(rowSums(object$Ki), kx))^2/sum(object$Ki))
    # }
    sd2 <- pmax(0, sd2)
    warning("Numerical errors caused some negative predictive variances to be thresholded to zero. Consider using ginv via rebuild.CRNGP")
  }
  
  if(!is.null(xprime)){
    kxprime <- object$nu_hat * cov_gen(X1 = object$X0, X2 = xprime, theta = object$theta[1:d], type = object$covtype)
    tmp <- outer(object$S0, sp, "==")
    ksprime <- matrix(object$rho, nrow(object$X0), nrow(xprime))
    ksprime[tmp] <- 1
    kxprime <- kxprime * ksprime
    if(!is.null(object$T0)){
      kxprime <- kronecker(kt, kxprime)
    }
    
    kxxprime <- object$nu_hat * cov_gen(X1 = x, X2 = xprime, theta = object$theta[1:d], type = object$covtype)
    ksxxprime <- matrix(object$rho, nrow(x), nrow(xprime))
    tmp <- outer(s, sp, "==")
    ksxxprime[tmp] <- 1
    kxxprime <- kxxprime * ksxxprime
    if(!is.null(object$T0)){
      kxxprime <- kronecker(kt, kxxprime)
    }
    
    
    if(object$trendtype == 'SK'){
      if(nrow(x) < nrow(xprime)){
        cov <- kxxprime - kx %*% object$Ki %*% kxprime 
      }else{
        cov <- kxxprime - kx %*% (object$Ki %*% kxprime)
      }
    }else{
      if(nrow(x) < nrow(xprime)){
        cov <- kxxprime - kx %*% object$Ki %*% kxprime + crossprod(1 - tcrossprod(rowSums(object$Ki), kx), 1 - rowSums(object$Ki) %*% kxprime)/sum(object$Ki)
      }else{
        cov <- kxxprime - kx %*% (object$Ki %*%  kxprime) + crossprod(1 - tcrossprod(rowSums(object$Ki), kx), 1 - rowSums(object$Ki) %*% kxprime)/sum(object$Ki)
      }
    }
  }else{
    cov = NULL
  }
  
  if(!is.null(object$T0)){
    mean <- matrix(mean, ncol = length(object$T0))
    sd2 <- matrix(sd2, ncol = length(object$T0))
  }
  
  return(list(mean = mean, sd2 = sd2, nugs = nugs, cov = cov))
}

if(!isGeneric("simul")) {
  setGeneric(name = "simul",
             def = function(object, ...) standardGeneric("simul")
  )
}

#' Conditional simulation for CRNGP
#' @title Conditional simulation for CRNGP
#' @param object \code{CRNGP} model
#' @param Xgrid matrix of (x, seed) locations where the simulation is performed. 
#' Where all design locations are matched with all seed values. In particular, it is assumed that each unique x values is matched with all seeds before going to the next x value.
#' The last column MUST correspond to seeds values. \code{Xgrid} must also contain the evaluated designs (e.g., in model$X0)
#' @param ids vector of indices corresponding to observed values in \code{Xgrid}
#' @param nsim number of simulations to return
#' @param eps jitter used in the Cholesky decomposition of the covariance matrix for numerical stability
#' @param seqseeds is the seed sequence repeated (e.g., 1 2 3 1 2 3), else it is assumed to be ordered (e.g., 1 1 2 2 3 3)
#' @param check if \code{TRUE}, check that Xgrid has the proper structure (slower)
#' @export
#' @return Conditional simulation matrix.
#' 
simul <- function (object, Xgrid, ids, nsim, eps, seqseeds, check) {
  UseMethod("simul", object)
}

#' Fast conditional simulation for a CRNGP model
#' @param object a \code{CRNGP} model obtained with \code{\link[hetGP]{mleCRNGP}}
#' @param Xgrid matrix of (x, seed) locations where the simulation is performed.
#' The last column MUST correspond to seeds values. \code{Xgrid} must also contain the evaluated designs (e.g., in \code{object$X0}).
#' All design locations are matched with all seed values, either by increasing seed values or repeating the seed sequence.
#' @param ids vector of indices corresponding to observed values in \code{Xgrid}
#' @param nsim number of simulations to return
#' @param eps jitter used in the Cholesky decomposition of the covariance matrix for numerical stability
#' @param seqseeds is the seed sequence repeated (e.g., 1 2 3 1 2 3), else it is assumed to be ordered (e.g., 1 1 2 2 3 3)
#' @param check if \code{TRUE}, check that \code{Xgrid} has the proper structure (slower)
#' @returns A matrix of size \code{nrow(Xgrid) x nsim}.
#' @references 
#' 
#' Chiles, J. P., & Delfiner, P. (2012). Geostatistics: modeling spatial uncertainty (Vol. 713). John Wiley & Sons. \cr \cr
#' 
#' Chevalier, C.; Emery, X.; Ginsbourger, D.
#' Fast Update of Conditional Simulation Ensembles
#' Mathematical Geosciences, 2014
#' @examples
#' \dontrun{
#' ##------------------------------------------------------------
#' ## Example: Homoskedastic GP modeling on 2d sims
#' ##------------------------------------------------------------
#' set.seed(2)
#' nx <- 31
#' ns <- 5
#' d <- 2
#' x <- as.matrix(expand.grid(seq(0,1, length.out = nx), seq(0,1, length.out = nx)))
#' s <- matrix(seq(1, ns, length.out = ns))
#' Xgrid <- as.matrix(expand.grid(seq(1, ns, length.out = ns), seq(0,1, length.out = nx), 
#'                                seq(0,1, length.out = nx)))
#' Xgrid <- Xgrid[,c(2, 3, 1)]
#' g <- 1e-6
#' theta <- c(0.2, 0.5)
#' KX <- cov_gen(x, theta = theta)
#' rho <- 0.33
#' KS <- matrix(rho, ns, ns)
#' diag(KS) <- 1
#' 
#' YY <- MASS::mvrnorm(n = 1, mu = rep(0, nx*nx*ns), Sigma = kronecker(KX, KS) + g * diag(nx*nx*ns))
#' YYmat <- matrix(YY, ns, nx*nx)
#' filled.contour(matrix(YYmat[1,], nx))
#' filled.contour(matrix(YYmat[2,], nx))
#' 
#' ids <- sample(1:nrow(Xgrid), 80)
#' X0 <- Xgrid[ids,]
#' Y0 <-  YY[ids]
#' 
#' # For 3d visualization
#' # library(rgl)
#' # plot3d(Xgrid[,1], Xgrid[,2], YY, col = 1 + (Xgrid[,3] - 1) %% 6)
#' # points3d(X0[,1], X0[,2], Y0, size = 10, col = 1 + ((X0[,3] - 1) %% 6))
#' 
#' model <- mleCRNGP(X0, Y0, known = list(g = 1e-6))
#' 
#' preds <- predict(model, x = Xgrid, xprime = Xgrid)
#' # surface3d(unique(Xgrid[1:nx^2,1]),unique(Xgrid[,2]), matrix(YY[Xgrid[,3]==1], nx), 
#' #   front = "lines", back = "lines")
#' # aspect3d(1, 1, 1)
#' # surface3d(unique(Xgrid[1:nx^2,1]),unique(Xgrid[,2]), matrix(preds$mean[Xgrid[,3]==1], nx), 
#' #   front = "lines", back = "lines", col = "red")
#' 
#' # Conditional realizations (classical way)
#' set.seed(2)
#' t0 <- Sys.time()
#' SigmaCond <- 1/2 * (preds$cov + t(preds$cov))
#' sims <- t(chol(SigmaCond + diag(sqrt(.Machine$double.eps), nrow(Xgrid)))) %*% rnorm(nrow(Xgrid))
#' sims <- sims + preds$mean
#' print(difftime(Sys.time(), t0))
#' # sims <- MASS::mvrnorm(n = 1, mu = preds$mean, Sigma = 1/2 * (preds$cov + t(preds$cov)))
#' # plot3d(X0[,1], X0[,2], Y0, size = 10, col = 1 + ((X0[,3] - 1) %% 6))
#' # surface3d(unique(x[,1]), unique(x[,2]), matrix(sims[Xgrid[,3] == 1], nx), col = 1, 
#' #   front = "lines", back = "lines")
#' # surface3d(unique(x[,1]), unique(x[,2]), matrix(sims[Xgrid[,3] == 2], nx), col = 2, 
#' #   front = "lines", back = "lines")
#'
#' # Alternative for conditional realizations 
#' # (note: here the design points are part of the simulation points)
#' set.seed(2)
#' t0 <- Sys.time()
#' condreas <- simul(model, Xgrid, ids = ids)
#' print(difftime(Sys.time(), t0))
#' # plot3d(X0[,1], X0[,2], Y0, size = 10, col = 1 + ((X0[,3] - 1) %% 6))
#' # surface3d(unique(x[,1]), unique(x[,2]), matrix(condreas[Xgrid[,3] == 1], nx), col = 1, 
#' #   front = "lines", back = "lines")
#' # surface3d(unique(x[,1]), unique(x[,2]), matrix(condreas[Xgrid[,3] == 2], nx), col = 2, 
#' #   front = "lines", back = "lines")
#'
#' # Alternative using ordered seeds:
#' Xgrid2 <- as.matrix(expand.grid(seq(0,1, length.out = nx), 
#'   seq(0,1, length.out = nx), seq(1, ns, length.out = ns)))
#' condreas2 <- simul(model, Xgrid2, ids = ids, seqseeds = FALSE)
#' 
#' ## Check that values at X0 are coherent:
#' # condreas[ids,1] - Y0
#' # sims[ids,1] - Y0
#' 
#' ## Check that the empirical mean/covariance is correct
#' condreas2 <- simul(model, Xgrid, ids = ids, nsim = 1000)
#' print(range(rowMeans(condreas2) - preds$mean))
#' print(range(cov(t(condreas2)) - preds$cov))
#' }
#' @method simul CRNGP
#' @export
simul.CRNGP <- function(object, Xgrid, ids = NULL, nsim = 1, eps = sqrt(.Machine$double.eps), 
                        seqseeds = TRUE, check = TRUE){
  d <- ncol(Xgrid)
  
  if(is.null(ids)){
    ids <- find_corres(Xgrid, object$X0)
  }
  
  if(!is.null(object$T0)) stop("Not implemented yet.")
  if(object$trendtype == "OK") warning("The uncertainty on the trend is not taken into account, perhaps fit the model using a fixed known beta0")
  
  Xgrid0 <- unique(Xgrid[, -d, drop = FALSE])
  sgrid <- unique(Xgrid[, d])
  ns <- length(sgrid)
  kx <- object$nu_hat * cov_gen(X1 = Xgrid[, -d, drop = FALSE], X2 = object$X0, theta = object$theta[1:(d-1)], type = object$covtype)
  tmp <- outer(Xgrid[, d], object$S0, "==")
  ks <- matrix(object$rho, nrow(Xgrid), nrow(object$X0))
  ks[tmp] <- 1
  kx <- kx * ks
  kw <- kx %*% object$Ki / object$nu_hat # Kriging weights
  
  kxxprime <- object$nu_hat * cov_gen(X1 = Xgrid0, X2 = Xgrid0, theta = object$theta, type = object$covtype)
  ksxxprime <- matrix(object$rho, ns, ns)
  tmp <- outer(sgrid, sgrid, "==")
  ksxxprime[tmp] <- 1
  
  if(seqseeds){
    ## Check that the covariance matrix is correct
    if(check){
      print("The check of Xgrid is on (Turn off for speed).")
      kxprime <- object$nu_hat * cov_gen(X1 = Xgrid[, -d, drop=F], X2 = Xgrid[, -d, drop = F], theta = object$theta, type = object$covtype)
      tmp <- outer(Xgrid[,d], Xgrid[,d], "==")
      ksprime <- matrix(object$rho, nrow(Xgrid), nrow(Xgrid))
      ksprime[tmp] <- 1
      kxprime <- kxprime * ksprime
      if(max(abs(as.vector(kronecker(kxxprime, ksxxprime) - kxprime))) > eps) 
        stop("Something off with Xgrid, perhaps look at the example.")
    }
    
    cholK <- kronecker(chol(kxxprime + diag(eps, nrow(Xgrid0))), chol(ksxxprime)) # Verif: cK2 <- chol(kronecker(kxxprime, ksxxprime) + diag(eps, nrow(Xgrid)))
  }else{
    ## Check that the covariance matrix is correct
    if(check){
      print("The check of Xgrid is on.")
      kxprime <- object$nu_hat * cov_gen(X1 = Xgrid[, -d, drop=F], X2 = Xgrid[, -d, drop = F], theta = object$theta, type = object$covtype)
      tmp <- outer(Xgrid[,d], Xgrid[,d], "==")
      ksprime <- matrix(object$rho, nrow(Xgrid), nrow(Xgrid))
      ksprime[tmp] <- 1
      kxprime <- kxprime * ksprime
      if(max(abs(as.vector(kronecker(ksxxprime, kxxprime) - kxprime))) > eps) 
        stop("Something off with Xgrid, perhaps look at the example.")
    }
    
    cholK <- kronecker(chol(ksxxprime), chol(kxxprime + diag(eps, nrow(Xgrid0)))) # Verif: cK2 <- chol(kronecker(ksxxprime, kxxprime) + diag(eps, nrow(Xgrid)))
    
  }
  reas <- object$beta0 + t(cholK) %*% matrix(rnorm(nrow(Xgrid0) * ns * nsim), ncol = nsim)
  mean <- as.vector(object$beta0 + kw %*% (object$Z - object$beta0))

  for(i in 1:nsim){
    altmean <- as.vector(object$beta0 + kw %*% (reas[ids, i] - object$beta0)) # maybe it is not ids (?)
    reas[, i] <- reas[, i] - altmean + mean
  }
  return(reas)
}

## Model: noisy observations with unknown homoskedastic noise and common random numbers
## All times are evaluated, leading to a Kronecker Structure
## K = nu^2 * (C + g * I) with C = Cx + Cs
# X0 unique designs matrix (last entry is the seed), here X0 = kronecker(XS, Ts)
# Z0 averaged observations at X0
# Z observations vector (all observations)
# mult number of replicates at each unique design
# theta vector of lengthscale hyperparameters (or one for isotropy)
# g noise variance for the process
# rho crn parameter
# beta0 trend
# @return loglikelihood value
logLikHomCRNT <- function(X0, S0, T0, Z, theta, g, rho, stype, beta0 = NULL, covtype = "Gaussian", eps = sqrt(.Machine$double.eps), env = NULL){
  n <- nrow(X0)
  N <- length(Z)
  d <- ncol(X0) # spatial dims
  
  # Temporarily store Cholesky transform of K in Ki
  Cx <- cov_gen(X1 = X0, theta = theta[1:d], type = covtype)
  if(is.null(env) || is.null(env$idS)){
    idS <- outer(S0, S0, "==")
    if(!is.null(env)) 
      env$idS <- idS
  }else{
    idS <- env$idS
  }
  Cs <- matrix(rho, n, n)
  Cs[idS] <- 1
  Ct <- cov_gen(X1 = T0, theta = theta[d+1], type = covtype)
  # C <- kronecker(Ct, Cx * Cs)
  
  if(!is.null(env)){
    env$Cx <- Cx
    env$Cs <- Cs
    env$Ct <- Ct
    # env$C <- C
  }
  
  # ## Verif
  # C2 <- matrix(NA, N, N)
  # for(i in 1:nrow(X)){
  #   for(j in i:nrow(X)){
  #     if(X[i, 2] == X[j, 2]) tmp <- 1 else tmp <- rho
  #     C2[i, j] <- C2[j, i] <- cov_gen(X[i,c(1,3), drop = F], X[j,c(1,3), drop = F], theta = theta) * tmp
  #   }
  # }
  # Ki2 <- chol(C2 + diag(eps + g, N))
  # ldetKi2 <- - 2 * sum(log(diag(Ki2))) # log determinant from Cholesky
  # Ki2 <- chol2inv(Ki2)
  
  SCxs <- svd(Cx*Cs) # C
  SCt <- svd(Ct) # R
  Ki <- tcrossprod(kronecker(SCt$u, SCxs$u) * rep(1/(kronecker(SCt$d, SCxs$d) + g), each = N), kronecker(SCt$u, SCxs$u))
  
  ldetKi <- -sum(log(kronecker(SCt$d, SCxs$d) + g))
  
  if(!is.null(env)) env$Ki <- Ki
  
  if(is.null(beta0))
    beta0 <- drop(colSums(Ki) %*% as.vector(Z) / sum(Ki))
  
  psi <- drop(crossprod(as.vector(Z) - beta0, Ki) %*% (as.vector(Z) - beta0)) / N
  
  loglik <- -N/2 * log(2*pi) - N/2 * log(psi) + 1/2 * ldetKi - N/2
}


# derivative of log-likelihood for logLikHom with respect to theta with all observations (Gaussian kernel)
## All times are evaluated, leading to a Kronecker Structure
## Model: noisy observations with unknown homoskedastic noise
## K = nu^2 * (C + g * I)
# X0  design matrix (no replicates)
# Z0 averaged observations
# Z observations vector (all observations)
# mult number of replicates per unique design point
# theta vector of lengthscale hyperparameters (or one for isotropy)
# g noise variance for the process
# beta0 trend
## @return gradient with respect to theta and g
dlogLikHomCRNT <- function(X0, T0, S0, Z, theta, g, rho, stype, beta0 = NULL, covtype = "Gaussian",
                           eps = sqrt(.Machine$double.eps), components = c("theta", "g", "rho"), env = NULL){
  n <- nrow(X0)
  N <- length(Z)
  d <- ncol(X0) # spatial dims
  
  
  if(!is.null(env)){
    # C <- env$C
    Cx <- env$Cx
    Cs <- env$Cs
    Ct <- env$Ct
    Ki <- env$Ki
  }else{
    Cx <- cov_gen(X1 = X0, theta = theta[1:d], type = covtype)
    idS <- outer(S0, S0, "==")
    Cs <- matrix(rho, n, n)
    Cs[idS] <- 1
    Ct <- cov_gen(X1 = T0, theta = theta[d+1], type = covtype)
    SCxs <- svd(Cx*Cs) # C
    SCt <- svd(Ct) # R
    Ki <- tcrossprod(kronecker(SCt$u, SCxs$u) * rep(1/(kronecker(SCt$d, SCxs$d) + g), each = N), kronecker(SCt$u, SCxs$u))
  }
  
  # ## Verif
  # C2 <- matrix(NA, N, N)
  # for(i in 1:nrow(X)){
  #   for(j in i:nrow(X)){
  #     if(X[i, 2] == X[j, 2]) tmp <- 1 else tmp <- rho
  #     C2[i, j] <- C2[j, i] <- cov_gen(X[i, c(1, 3), drop = F], X[j, c(1, 3), drop = F], theta = theta) * tmp
  #   }
  # }
  # Ki2 <- chol(C2 + diag(eps + g, N))
  # ldetKi2 <- - 2 * sum(log(diag(Ki2))) # log determinant from Cholesky
  # Ki2 <- chol2inv(Ki2)
  
  if(is.null(beta0))
    beta0 <- drop(colSums(Ki) %*% as.vector(Z) / sum(Ki))
  
  Z <- as.vector(Z) - beta0
  
  KiZ <- Ki %*% Z ## to avoid recomputing  
  
  psi <- drop(crossprod(Z, KiZ))
  
  tmp1 <- tmp2 <- tmp3 <- NULL
  
  # First component, derivative with respect to theta
  if("theta" %in% components){
    tmp1 <- rep(NA, length(theta))
    if(d==1){
      
      # Kfun <- function(X, theta, rho, i, j){
      #   d <- ncol(X)
      #   kx <- cov_gen(matrix(X[,-d,drop = F]), matrix(X[,-d, drop = F]), theta = theta)
      #   tmp <- outer(X[,d], X[,d], "==")
      #   ks <- matrix(rho, nrow(X), nrow(X))
      #   ks[tmp] <- 1
      #   res <- kx * ks
      #   return(res[i, j])
      # }
      # grad(Kfun, X = X0, x = theta, rho = rho, i = 1, j = 2)
      
      dC_dthetak <- partial_cov_gen(X1 = X0, theta = theta[1], type = covtype, arg = "theta_k") * Cx * Cs
      dC_dthetak <- kronecker(Ct, dC_dthetak)
      tmp1 <- N/2 * crossprod(KiZ, dC_dthetak) %*% KiZ / psi - 1/2 * trace_sym(Ki, dC_dthetak)
    }else{
      for(i in 1:(length(theta) - 1)){
        dC_dthetak <- partial_cov_gen(X1 = X0[,i, drop = F], theta = theta[i], type = covtype, arg = "theta_k") * Cx * Cs
        dC_dthetak <- kronecker(Ct, dC_dthetak)
        tmp1[i] <- N/2 * crossprod(KiZ, dC_dthetak) %*% KiZ / psi - 1/2 * trace_sym(Ki, dC_dthetak)
      }
    } 
    # Now for time
    dC_dthetak <- partial_cov_gen(X1 = T0, theta = theta[d+1], type = covtype, arg = "theta_k") * Ct
    dC_dthetak <- kronecker(dC_dthetak, Cx * Cs)
    tmp1 <- c(tmp1, N/2 * crossprod(KiZ, dC_dthetak) %*% KiZ / psi - 1/2 * trace_sym(Ki, dC_dthetak))
  }
  
  # Second component derivative with respect to g
  if("g" %in% components)
    tmp2 <- N/2 * sum(KiZ^2) / psi - 1/2 * sum(diag(Ki))
  
  if("rho" %in% components){
    dC_drho <- matrix(1, n, n)
    ids <- outer(S0, S0, "==")
    dC_drho[ids] <- 0
    dC_drho <- dC_drho * Cx
    dC_drho <- kronecker(Ct, dC_drho)
    tmp3 <- N/2 * crossprod(KiZ, dC_drho) %*% KiZ /psi - 1/2 * trace_sym(Ki, dC_drho)
    
  }
  
  return(c(tmp1, tmp2, tmp3))
}



#' Gaussian process regression under homoskedastic noise based on maximum likelihood estimation of the 
#' hyperparameters. This function is enhanced to deal with replicated observations. [Deprecated: it is included in CRNGP now]
#' @title Gaussian process modeling with homoskedastic noise
#' @param X matrix of all designs, one per row, or list with elements:
#' \itemize{
#'   \item \code{X0} matrix of unique design locations, one point per row
#'   \item \code{Z0} vector of averaged observations, of length \code{nrow(X0)}
#'   \item \code{mult} number of replicates at designs in \code{X0}, of length \code{nrow(X0)}
#' } 
#' @param Z vector of all observations. If using a list with \code{X}, \code{Z} has to be ordered with respect to \code{X0}, and of length \code{sum(mult)}
#' @param lower,upper optional bounds for the \code{theta} parameter (see \code{\link[hetGP]{cov_gen}} for the exact parameterization).
#' In the multivariate case, it is possible to give vectors for bounds (resp. scalars) for anisotropy (resp. isotropy) 
#' @param known optional list of known parameters, e.g., \code{beta0}, \code{theta} or \code{g}
#' @param covtype covariance kernel type, either 'Gaussian', 'Matern5_2' or 'Matern3_2', see \code{\link[hetGP]{cov_gen}}
#' @param noiseControl list with element , 
#' \itemize{
#' \item \code{g_bounds}, vector providing minimal and maximal noise to signal ratio
#' } 
#' @param init optional list specifying starting values for MLE optimization, with elements:
#' \itemize{
#'  \item \code{theta_init} initial value of the theta parameters to be optimized over (default to 10\% of the range determined with \code{lower} and \code{upper})
#'  \item \code{g_init} initial value of the nugget parameter to be optimized over (based on the variance at replicates if there are any, else \code{0.1})
#' }
#' @param maxit maximum number of iteration for L-BFGS-B of \code{\link[stats]{optim}}
#' @param eps jitter used in the inversion of the covariance matrix for numerical stability
#' @param settings list with argument \code{return.Ki}, to include the inverse covariance matrix in the object for further use (e.g., prediction).
#' Arguments \code{factr} (default to 1e9) and \code{pgtol} are available to be passed to \code{control} for L-BFGS-B in \code{\link[stats]{optim}} (for the joint likelihood only). 
#' @return a list which is given the S3 class "\code{CRNGP}", with elements:
#' \itemize{
#' \item \code{theta}: maximum likelihood estimate of the lengthscale parameter(s),
#' \item \code{g}: maximum likelihood estimate of the nugget variance,
#' \item \code{trendtype}: either "\code{SK}" if \code{beta0} is given, else "\code{OK}" 
#' \item \code{beta0}: estimated trend unless given in input,
#' \item \code{nu_hat}: plugin estimator of the variance,
#' \item \code{ll}: log-likelihood value,
#' \item \code{X0}, \code{Z0}, \code{Z}, \code{mult}, \code{eps}, \code{covtype}: values given in input,
#' \item \code{call}: user call of the function
#' \item \code{used_args}: list with arguments provided in the call
#' \item \code{nit_opt}, \code{msg}: \code{counts} and \code{msg} returned by \code{\link[stats]{optim}}
#' \item \code{Ki}: inverse covariance matrix (not scaled by \code{nu_hat}) (if \code{return.Ki} is \code{TRUE} in \code{settings})
#' \item \code{time}: time to train the model, in seconds.
#' 
#'}
#' @details
#' The global covariance matrix of the model is parameterized as \code{nu_hat * (C + g * diag(1/mult)) = nu_hat * K},
#' with \code{C} the correlation matrix between unique designs, depending on the family of kernel used (see \code{\link[hetGP]{cov_gen}} for available choices) and values of lengthscale parameters.
#' \code{nu_hat} is the plugin estimator of the variance of the process.
#' 
#' It is generally recommended to use \code{\link[hetGP]{find_reps}} to pre-process the data, to rescale the inputs to the unit cube and to normalize the outputs.
#' 
#' @seealso \code{\link[hetGP]{predict.CRNGP}} for predictions, 
#' \code{summary} and \code{plot} functions are available as well. 
#' \code{\link[hetGP]{mleHomTP}} provide a Student-t equivalent.
#' @references 
#' M. Binois, Robert B. Gramacy, M. Ludkovski (2018), Practical heteroskedastic Gaussian process modeling for large simulation experiments,
#' Journal of Computational and Graphical Statistics, 27(4), 808--821.\cr 
#' Preprint available on arXiv:1611.05902. \cr \cr
#' @noRd
## ' @importFrom numDeriv hessian
#' @examples
#' ##------------------------------------------------------------
#' ## Example 1: Homoskedastic GP modeling on 1d sims (+ time)
#' ##------------------------------------------------------------
#' set.seed(42)
#' nx <- 50
#' nt <- 30
#' ns <- 20
#' x <- matrix(sort(seq(0,1, length.out = nx)), nx)
#' s <- matrix(sort(seq(1, ns, length.out = ns)))
#' t <- matrix(sort(seq(0, 1, length.out = nt)), nt)
#' g <- 1e-3
#' theta <- c(0.01, 0.02)
#' KX <- cov_gen(x, theta = theta[1])
#' KT <- cov_gen(t, theta = theta[2])
#' rho <- 0.3
#' KS <- matrix(rho, ns, ns)
#' diag(KS) <- 1
#' XST <- as.matrix(expand.grid(x, s, t))
#' KM <- kronecker(KT, kronecker(KS, KX))
#' # KM2 <- matrix(NA, nrow(XST), nrow(XST))
#' # for(i in 1:nrow(XST)){
#' #  for(j in i:nrow(XST)){
#' #   if(XST[i, 2] == XST[j, 2]) tmp <- 1 else tmp <- rho
#' #   KM2[i, j] <- cov_gen(XST[i,c(1,3), drop = F], XST[j,c(1,3), drop = F], theta = theta) * tmp
#' #   KM2[j, i] <- KM2[i, j]
#' #   }
#' # }
#' # range(abs(KM - KM2))
#' 
#' ## YY <- MASS::mvrnorm(n = 1, mu = rep(0, nx*nt*ns), Sigma = KM + g * diag(nx*nt*ns))
#' 
#' Kmc <- kronecker(chol(KT), kronecker(chol(KS), chol(KX)))
#' YY <- t(Kmc) %*% rnorm(nrow(Kmc))
#' 
#' ninit <- 40
#' XS <- as.matrix(expand.grid(x, s))
#' ids <- sort(sample(1:nrow(XS), ninit))
#' XS0 <- cbind(XS[ids[rep(1:ninit, each = nt)],], rep(t[,1], times = ninit))
#' X0 <- XST[which(duplicated(rbind(XST, XS0), fromLast = TRUE)),]
#' Y0 <-  YY[which(duplicated(rbind(XST, XS0), fromLast = TRUE))]
#' 
#' model <- mleCRNGPT(X0, Y0)
#' 
#' preds <- predict(model, x = XST, xprime = XST)
#' 
#' # Conditional realizations
#' sims <- MASS::mvrnorm(n = 1, mu = preds$mean, Sigma = 1/2 * (preds$cov + t(preds$cov)))
#'
#' # compare with regular CRN GP
#' mref <- mleCRNGP(X = X0[, c(1, 3, 2)], Z = Y0)#, known = list(theta = model$theta, g = model$g)) 
#' pref <- predict(mref, x = XST[, c(1, 3, 2)], xprime = XST[, c(1, 3, 2)])
#' 
#' range(preds$mean - pref$mean)
#' range(preds$cov - pref$cov)
#' range(preds$sd2 - pref$sd2)
#' 
#' # Higher dimensional example
#' d <- 3
#' nx <- 50
#' nt <- 20
#' ns <- 20
#' x <- matrix(matrix(runif(d * nx)), nx)
#' s <- matrix(sort(seq(1, ns, length.out = ns)))
#' t <- matrix(sort(seq(0, 1, length.out = nt)), nt)
#' 
#' function(X){
#'   x <- X[1:3]
#'   s <- X[4]
#'   t <- X[5]
#'   return(sin(sum(x) * sin(s) * sin(t)))
#' }
#' 
#' 
#' 
#' 
#' 
mleCRNGPT <- function(X, Z, lower = NULL, upper = NULL, known = NULL,
                      noiseControl = list(g_bounds = c(sqrt(.Machine$double.eps)*10, 1e2),
                                          rho_bounds = c(0.001, 0.9)),
                      init = NULL,
                      covtype = c("Gaussian", "Matern5_2", "Matern3_2"),
                      maxit = 100, eps = sqrt(.Machine$double.eps), settings = list(return.Ki = TRUE, factr = 1e7)){
  
  # if(typeof(X)=="list"){
  #   X0 <- X$X0
  #   Z0 <- X$Z0
  #   mult <- X$mult
  #   if(sum(mult) != length(Z)) stop("Length(Z) should be equal to sum(mult)")
  #   if(is.null(dim(X0))) X0 <- matrix(X0, ncol = 1)
  #   if(length(Z0) != nrow(X0)) stop("Dimension mismatch between Z0 and X0")
  # }else{
  #   if(is.null(dim(X))) X <- matrix(X, ncol = 1)
  #   if(nrow(X) != length(Z)) stop("Dimension mismatch between Z and X")
  #   elem <- find_reps(X, Z, return.Zlist = F)
  #   X0 <- elem$X0
  #   Z0 <- elem$Z0
  #   Z <- elem$Z
  #   mult <- elem$mult
  # }
  
  T0 <- unique(X[, ncol(X), drop = F])
  elem <- find_reps(X = X[,-ncol(X), drop = F], Z = Z)
  S0 <- elem$X0[, ncol(elem$X0), drop = F]
  X0 <- elem$X0[, - ncol(elem$X0), drop = F]
  d <- ncol(X0) # spatial dim
  
  # Some checks
  if(any(elem$mult != elem$mult[1])) warning("Some X do not have the same number of time observations.")
  
  covtype <- match.arg(covtype)
  
  if(is.null(lower) || is.null(upper)){
    auto_thetas1 <- auto_bounds(X = X0, covtype = covtype)
    auto_thetas2 <- auto_bounds(X = T0, covtype = covtype)
    if(is.null(lower)) lower <- c(auto_thetas1$lower, auto_thetas2$lower)
    if(is.null(upper)) upper <- c(auto_thetas1$upper, auto_thetas2$upper)
    if(is.null(known[["theta"]]) && is.null(init$theta)) init$theta <- sqrt(upper * lower)
  }
  if(length(lower) != length(upper)) stop("upper and lower should have the same size")
  
  ## Save time to train model
  tic <- proc.time()[3]
  
  if(is.null(settings$return.Ki)) settings$return.Ki <- TRUE
  if(is.null(noiseControl$g_bounds)) noiseControl$g_bounds <- c(sqrt(.Machine$double.eps), 1e2)
  if(is.null(noiseControl$rho_bounds)) noiseControl$rho_bounds <- c(0, 0.9)
  
  
  g_min <- noiseControl$g_bounds[1]
  g_max <- noiseControl$g_bounds[2]
  
  beta0 <- known$beta0
  
  N <- length(Z)
  n <- nrow(X0)
  
  if(is.null(n))
    stop("X0 should be a matrix. \n")
  
  if(is.null(known[["theta"]]) && is.null(init$theta)) init$theta <- 0.9 * lower + 0.1 * upper # useful for mleHetGP
  if(is.null(known$g) && is.null(init$g)){
    # if(any(mult > 2)) init$g <- mean((fast_tUY2(mult, (Z - rep(Z0, times = mult))^2)/mult)[which(mult > 2)])/var(Z0) else 
    init$g <- 0.1
  }
  if(is.null(known$rho)) init$rho <- 0.1
  
  trendtype <- 'OK'
  if(!is.null(beta0))
    trendtype <- 'SK'
  
  ## General definition of fn and gr
  fn <- function(par, X0, T0, S0, Z, beta0, theta, g, rho, env){
    idx <- 1 # to store the first non used element of par
    
    if(is.null(theta)){
      theta <- par[1:length(init$theta)]
      idx <- idx + length(init$theta)
    }
    
    if(is.null(g)){
      g <- par[idx]
      idx <- idx + 1
    }
    
    if(is.null(rho)){
      rho <- par[idx]
    }
    
    loglik <- logLikHomCRNT(X0 = X0, S0 = S0, T0 = T0, Z = Z, theta = theta, g = g, rho = rho, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    
    # library(numDeriv)
    # grad(logLikHomCRNT, x = theta, X0 = X0, S0 = S0, T0 = T0, Z = Z, g = g, rho = rho, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    # grad(logLikHomCRNT, x = g, X0 = X0, S0 = S0, T0 = T0, Z = Z, theta = theta, rho = rho, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    # grad(logLikHomCRNT, x = rho, X0 = X0, S0 = S0, T0 = T0, Z = Z, theta = theta, g = g, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    # dlogLikHomCRNT(X0 = X0, S0 = S0, T0 = T0, Z = Z, theta = theta, g = g, rho = rho, beta0 = beta0, covtype = covtype, eps = eps, env = env)
    
    if(!is.null(env) && !is.na(loglik)){
      if(is.null(env$max_loglik) || loglik > env$max_loglik){
        env$max_loglik <- loglik
        env$arg_max <- par
      }
    } 
    
    return(loglik)
  }
  
  gr <- function(par, X0, T0, S0, Z, beta0, theta, g, rho, env){
    idx <- 1
    components <- NULL
    
    if(is.null(theta)){
      theta <- par[1:length(init$theta)]
      idx <- idx + length(init$theta)
      components <- "theta"
    }
    
    if(is.null(g)){
      g <- par[idx]
      idx <- idx + 1
      components <- c(components, "g")
    }
    
    if(is.null(rho)){
      rho <- par[idx]
      components <- c(components, "rho")
    }
    return(dlogLikHomCRNT(X0 = X0, T0 = T0, S0 = S0, Z = Z, theta = theta, g = g, rho = rho, beta0 = beta0, covtype = covtype, eps = eps,
                          components = components, env = env))
  }
  
  ## All known
  envtmp <- environment()
  if(!is.null(known$g) && !is.null(known[["theta"]]) && !is.null(known$rho)){
    theta_out <- known[["theta"]]
    g_out <- known$g
    rho_out <- known$rho
    out <- list(value = logLikHomCRNT(X0 = X0, T0 = T0, S0 = S0, Z = Z, theta = theta_out, g = g_out, rho = rho_out, beta0 = beta0, covtype = covtype, eps = eps),
                message = "All hyperparameters given", counts = 0, time = proc.time()[3] - tic)
  }else{
    parinit <- lowerOpt <- upperOpt <- NULL
    if(is.null(known[["theta"]])){
      parinit <- init$theta
      lowerOpt <- c(lower)
      upperOpt <- c(upper)
    }
    if(is.null(known$g)){
      parinit <- c(parinit, init$g)
      lowerOpt <- c(lowerOpt, g_min)
      upperOpt <- c(upperOpt, g_max)
    }
    if(is.null(known$rho)){
      parinit <- c(parinit, init$rho)
      lowerOpt <- c(lowerOpt, noiseControl$rho_bounds[1])
      upperOpt <- c(upperOpt, noiseControl$rho_bounds[2])
    }
    
    out <- try(optim(par = parinit, fn = fn, gr = gr, method = "L-BFGS-B", lower = lowerOpt, upper = upperOpt, 
                     theta = known[["theta"]], g = known$g, rho = known$rho,
                     X0 = X0, T0 = T0, S0 = S0, Z = Z, beta0 = beta0,
                     control = list(fnscale = -1, maxit = maxit, factr = settings$factr, pgtol = settings$pgtol), env = envtmp))
    ## Catch errors when at least one likelihood evaluation worked
    if(is(out, "try-error"))
      out <- list(par = envtmp$arg_max, value = envtmp$max_loglik, counts = NA,
                  message = "Optimization stopped due to NAs, use best value so far")
    
    idx <- 1
    if(is.null(known[["theta"]])){
      theta_out <- out$par[1:length(init$theta)]
      idx <- idx + length(theta_out)
    } else  theta_out <- known[["theta"]]
    
    
    if(is.null(known$g)){
      g_out <- out$par[idx]
    } else g_out <- known$g
    
    if(is.null(known$rho)) rho_out <- out$par[length(out$par)] else rho_out <- known$rho
    
  }
  
  Cx <- cov_gen(X1 = X0, theta = theta_out[1:d], type = covtype)
  tmp <- outer(S0, S0, "==")
  Cs <- matrix(rho_out, n, n)
  Cs[tmp] <- 1
  Ct <- cov_gen(X1 = T0, theta = theta_out[d+1], type = covtype)
  SCxs <- svd(Cx*Cs) # C
  SCt <- svd(Ct) # R
  Ki <- tcrossprod(kronecker(SCt$u, SCxs$u) * rep(1/diag(kronecker(diag(SCt$d), diag(SCxs$d)) + g_out * diag(N)), each = N), kronecker(SCt$u, SCxs$u))  
  # Cx <- cov_gen(X1 = X0[,-d, drop = F], theta = theta_out, type = covtype)
  # tmp <- outer(X0[,d], X0[,d], "==")
  # Cs <- matrix(rho_out, n, n)
  # Cs[tmp] <- 1
  # C <- Cx * Cs
  # Ki <- chol2inv(chol(C + diag(eps + g_out / mult)))
  
  if(is.null(beta0))
    beta0 <- drop(colSums(Ki) %*% Z / sum(Ki))
  
  psi <- drop(crossprod(Z - beta0, Ki) %*% (Z - beta0))
  
  nu <- psi / N
  
  # # Get hessian of our cost function/
  # if (is.null(known[["theta"]])) {
  #   # Jacobian is more precise numerically but doesn't seem to work for some reason
  #   #hess <- jacobian(func = gr, x = out$par, theta = known[["theta"]], g = known$g,
  #   #             X0 = X0, Z0 = Z0, Z = Z, mult = mult, beta0 = beta0, env = envtmp)
  #   fwrap <- function(par, ...) fn(c(par, out$par[length(out$par)]), ...)
  #   hess <- hessian(func = fwrap, x = out$par[1:(length(out$par)-1)], theta = known[["theta"]], g = known$g,
  #                   X0 = X0, Z0 = Z0, Z = Z, mult = mult, beta0 = beta0, env = NULL)
  # } else {
  #   hess <- NULL
  # }
  
  res <- list(theta = theta_out, g = g_out, rho = rho_out, nu_hat = as.numeric(nu), ll = out$value, nit_opt = out$counts,
              beta0 = beta0, trendtype = trendtype, covtype = covtype, msg = out$message, eps = eps,
              X0 = X0, S0 = S0, T0 = T0, X = X, Z = Z, call = match.call(),
              used_args = list(lower = lower, upper = upper, known = known, noiseControl = noiseControl),
              time = proc.time()[3] - tic) # hess = hess)
  
  if(settings$return.Ki) res <- c(res, list(Ki = Ki))
  
  class(res) <- "CRNGPT"
  return(res)
}

#' @method summary CRNGPT
#' @noRd
summary.CRNGPT <- function(object,...){
  ans <- object
  class(ans) <- "summary.CRNGPT"
  ans
}

#' @noRd
print.summary.CRNGPT <- function(x, ...){
  
  cat("N = ", length(x$Z), " nxs = ", nrow(x$X0), " nt =  ", length(x$T0), "d = ", ncol(x$X0), "\n")
  cat(x$covtype, " covariance lengthscale values (x then t): ", x$theta, "\n")
  
  cat("Homoskedastic nugget value: ", x$g, "\n")
  
  cat("CRN Correlation: ", x$rho, "\n")
  
  cat("Variance/scale hyperparameter: ", x$nu_hat, "\n")
  
  if(x$trendtype == "SK"){
    cat("Given constant trend value: ", x$beta0, "\n")
  }else{
    cat("Estimated constant trend value: ", x$beta0, "\n")
  }
  
  cat("MLE optimization: \n", "Log-likelihood = ", x$ll, "; Nb of evaluations (obj, gradient) by L-BFGS-B: ", x$nit_opt, "; message: ", x$msg, "\n")
}

#' @method print CRNGPT
#' @noRd
print.CRNGPT <- function(x, ...){
  print(summary(x))
}


#' Gaussian process predictions using a homoskedastic noise GP object (of class \code{CRNGPT}) [Deprecated: it is included in CRNGP now]
#' @param x matrix of designs locations to predict at (one point per row)
#' @param object an object of class \code{CRNGPT}; e.g., as returned by \code{\link[hetGP]{mleCRNGPT}}
#' @param xprime optional second matrix of predictive locations to obtain the predictive covariance matrix between \code{x} and \code{xprime}
#' @param ... no other argument for this method
#' @return list with elements
#' \itemize{
#' \item \code{mean}: kriging mean;
#' \item \code{sd2}: kriging variance (filtered, e.g. without the nugget value)
#' \item \code{cov}: predictive covariance matrix between \code{x} and \code{xprime}
#' \item \code{nugs}: nugget value at each prediction location, for consistency with \code{\link[hetGP]{mleHomGP}}.
#' }
#' @importFrom MASS ginv
#' @details The full predictive variance corresponds to the sum of \code{sd2} and \code{nugs}. See \code{\link[hetGP]{mleHomGP}} for examples.
#' @method predict CRNGPT
#' @noRd
predict.CRNGPT <- function(object, x, xprime = NULL, ...){
  if(is.null(dim(x))){
    x <- matrix(x, nrow = 1)
    if(ncol(x) != ncol(object$X0)) stop("x is not a matrix")
  }
  
  if(!is.null(xprime) && is.null(dim(xprime))){
    xprime <- matrix(xprime, nrow = 1)
    if(ncol(xprime) != ncol(object$X0)) stop("xprime is not a matrix")
  }
  
  t0 <- unique(x[, ncol(x), drop = F])
  elem <- find_reps(X = x[,-ncol(x), drop = F], Z = rep(0, nrow(x)))
  s0 <- elem$X0[, ncol(elem$X0), drop = F]
  x0 <- elem$X0[, - ncol(elem$X0), drop = F]
  
  
  d <- ncol(object$X0) # spatial dim
  n <- nrow(object$X0)
  N <- length(object$Z)
  
  if(is.null(object$Ki)){
    Cx <- cov_gen(X1 = object$X0, theta = object$theta[1:d], type = object$covtype)
    tmp <- outer(object$S0, object$S0, "==")
    Cs <- matrix(object$rho, n, n)
    Cs[tmp] <- 1
    Ct <- cov_gen(X1 = object$T0, theta = object$theta[d+1], type = object$covtype)
    SCxs <- svd(Cx*Cs) # C
    SCt <- svd(Ct) # R
    Ki <- tcrossprod(kronecker(SCt$u, SCxs$u) * rep(1/(kronecker(SCt$d, SCxs$d) + object$g), each = N), kronecker(SCt$u, SCxs$u))  
  }
  
  object$Ki <- object$Ki/object$nu_hat
  
  kx <- object$nu_hat * cov_gen(X1 = x0, X2 = object$X0, theta = object$theta[1:d], type = object$covtype)
  tmp <- outer(s0, object$S0, "==")
  ks <- matrix(object$rho, nrow(x0), nrow(object$X0))
  ks[tmp] <- 1
  kt <- cov_gen(X1 = t0, X2 = object$T0, theta = object$theta[d+1], type = object$covtype)
  
  kx <- kronecker(kt, kx * ks)
  
  nugs <- rep(object$nu_hat * object$g, nrow(x0))
  mean <- as.vector(object$beta0 + kx %*% (object$Ki %*% (object$Z - object$beta0)))
  
  if(object$trendtype == 'SK'){
    sd2 <- as.vector(object$nu_hat - fast_diag(kx, tcrossprod(object$Ki, kx)))
  }else{
    sd2 <- as.vector(object$nu_hat - fast_diag(kx, tcrossprod(object$Ki, kx)) + (1 - tcrossprod(rowSums(object$Ki), kx))^2/sum(object$Ki))
  }
  
  ## In case of numerical errors, some sd2 values may become negative
  if(any(sd2 < 0)){
    sd2 <- pmax(0, sd2)
    warning("Numerical errors caused some negative predictive variances to be thresholded to zero. Consider using ginv via rebuild.CRNGP")
  }
  
  if(!is.null(xprime)){
    t0p <- unique(xprime[, ncol(xprime), drop = F])
    elem <- find_reps(X = xprime[,-ncol(xprime), drop = F], Z = rep(0, nrow(xprime)))
    s0p <- elem$X0[, ncol(elem$X0), drop = F]
    x0p <- elem$X0[, -ncol(elem$X0), drop = F]
    
    kxprime <- object$nu_hat * cov_gen(X1 = object$X0, X2 = x0p, theta = object$theta[1:d], type = object$covtype)
    tmp <- outer(object$S0, s0p, "==")
    ksprime <- matrix(object$rho, nrow(object$X0), nrow(x0p))
    ksprime[tmp] <- 1
    ktprime <- cov_gen(X1 = object$T0, X2 = t0p, theta = object$theta[d + 1], type = object$covtype)
    
    kxprime <- kronecker(ktprime, kxprime * ksprime)
    
    kxxprime <- object$nu_hat * cov_gen(X1 = x0, X2 = x0p, theta = object$theta[1:d], type = object$covtype)
    ksxxprime <- matrix(object$rho, nrow(x0), nrow(x0p))
    tmp <- outer(s0, s0p, "==")
    ksxxprime[tmp] <- 1
    ktxxprime <- cov_gen(t0, t0p, theta = object$theta[d + 1], type = object$covtype)
    kxxprime <- kronecker(ktxxprime, kxxprime * ksxxprime)
    
    if(object$trendtype == 'SK'){
      if(nrow(x0) < nrow(x0p)){
        cov <-  kxxprime - kx %*% object$Ki %*% kxprime 
      }else{
        cov <- kxxprime - kx %*% (object$Ki %*% kxprime)
      }
    }else{
      if(nrow(x0) < nrow(x0p)){
        cov <- kxxprime - kx %*% object$Ki %*% kxprime + crossprod(1 - tcrossprod(rowSums(object$Ki), kx), 1 - rowSums(object$Ki) %*% kxprime)/sum(object$Ki)
      }else{
        cov <- kxxprime - kx %*% (object$Ki %*% kxprime) + crossprod(1 - tcrossprod(rowSums(object$Ki), kx), 1 - rowSums(object$Ki) %*% kxprime)/sum(object$Ki)
      }
    }
  }else{
    cov = NULL
  }
  
  return(list(mean = mean, sd2 = sd2, nugs = nugs, cov = cov))
}


