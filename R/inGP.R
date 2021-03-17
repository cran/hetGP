# Functions to predict at noisy inputs

#' Gaussian process prediction prediction at a noisy input \code{x}, with centered Gaussian noise of variance \code{sigma_x}. 
#' Several options are available, with different efficiency/accuracy tradeoffs.
#' @param x design considered
#' @param model GP
#' @param sigma_x input variance
#' @param type available options include
#' \itemize{
#' \item \code{simple} relying on a corrective term, see (McHutchon2011);
#' \item \code{taylor} based on a Taylor expansion, see, e.g., (Girard2003);
#' \item \code{exact} for exact moments (only for the \code{Gaussian} covariance).
#' }
#' @note Beta version.
#' @references 
#' A. McHutchon and C.E. Rasmussen (2011), 
#' Gaussian process training with input noise, 
#' Advances in Neural Information Processing Systems, 1341-1349. \cr
#' 
#' A. Girard, C.E. Rasmussen, J.Q. Candela and R. Murray-Smith (2003), 
#' Gaussian process priors with uncertain inputs application to multiple-step ahead time series forecasting,
#' Advances in Neural Information Processing Systems, 545-552.
#' 
#' @export
#' @examples
#' ################################################################################
#' ### Illustration of prediction with input noise
#' ################################################################################
#' 
#' ## noise std deviation function defined in [0,1]
#' noiseFun <- function(x, coef = 1.1, scale = 0.25){
#'   if(is.null(nrow(x))) x <- matrix(x, nrow = 1)
#'   return(scale*(coef + sin(x * 2 * pi)))
#' }
#' 
#' ## data generating function combining mean and noise fields
#' ftest <- function(x, scale = 0.25){
#' if(is.null(nrow(x))) x <- matrix(x, ncol = 1)
#'   return(f1d(x) + rnorm(nrow(x), mean = 0, sd = noiseFun(x, scale = scale)))
#' }
#' 
#' ntest <- 101; xgrid <- seq(0,1, length.out = ntest); Xgrid <- matrix(xgrid, ncol = 1)
#' set.seed(42)
#' Xpred <- Xgrid[rep(1:ntest, each = 100),,drop = FALSE]
#' Zpred <- matrix(ftest(Xpred), byrow = TRUE, nrow = ntest)
#' n <- 10
#' N <- 20
#' X <- matrix(seq(0, 1, length.out = n))
#' if(N > n) X <- rbind(X, X[sample(1:n, N-n, replace = TRUE),,drop = FALSE])
#' X <- X[order(X[,1]),,drop = FALSE]
#' 
#' Z <- apply(X, 1, ftest)
#' par(mfrow = c(1, 2))
#' plot(X, Z, ylim = c(-10,15), xlim = c(-0.1,1.1))
#' lines(xgrid, f1d(xgrid))
#' lines(xgrid, drop(f1d(xgrid)) + 2*noiseFun(xgrid), lty = 3)
#' lines(xgrid, drop(f1d(xgrid)) - 2*noiseFun(xgrid), lty = 3)
#' model <- mleHomGP(X, Z, known = list(beta0 = 0))
#' preds <- predict(model, Xgrid)
#' lines(xgrid, preds$mean, col = "red", lwd = 2)
#' lines(xgrid, preds$mean - 2*sqrt(preds$sd2), col = "blue")
#' lines(xgrid, preds$mean + 2*sqrt(preds$sd2), col = "blue")
#' lines(xgrid, preds$mean - 2*sqrt(preds$sd2 + preds$nugs), col = "blue", lty = 2)
#' lines(xgrid, preds$mean + 2*sqrt(preds$sd2 + preds$nugs), col = "blue", lty = 2)
#' 
#' sigmax <- 0.1
#' X1 <- matrix(0.5)
#' 
#' lines(xgrid, dnorm(xgrid, X1, sigmax) - 10, col = "darkgreen")
#' 
#' # MC experiment
#' nmc <- 1000
#' XX <- matrix(rnorm(nmc, X1, sigmax))
#' pxx <- predict(model, XX)
#' YXX <- rnorm(nmc, mean = pxx$mean, sd = sqrt(pxx$sd2 + pxx$nugs))
#' points(XX, YXX, pch = '.')
#' 
#' hh <- hist(YXX, breaks = 51, plot = FALSE)
#' dd <- density(YXX)
#' plot(hh$density, hh$mids, ylim = c(-10, 15))
#' lines(dd$y, dd$x)
#' 
#' # GP predictions
#' pin1 <- pred_noisy_input(X1, model, sigmax^2, type = "exact")
#' pin2 <- pred_noisy_input(X1, model, sigmax^2, type = "taylor")
#' pin3 <- pred_noisy_input(X1, model, sigmax^2, type = "simple")
#' ygrid <- seq(-10, 15,, ntest)
#' lines(dnorm(ygrid, pin1$mean, sqrt(pin1$sd2)), ygrid, lty = 2, col = "orange")
#' lines(dnorm(ygrid, pin2$mean, sqrt(pin2$sd2)), ygrid, lty = 2, col = "violet")
#' lines(dnorm(ygrid, pin3$mean, sqrt(pin3$sd2)), ygrid, lty = 2, col = "grey")
#' abline(h = mean(YXX), col = "red") # empirical mean
#' 
#' par(mfrow = c(1, 1))
pred_noisy_input <- function(x, model, sigma_x, type = c("simple", "taylor", "exact")){
  if(type == "exact" & model$covtype == "Gaussian") 
    return(pred_input_var_exact(x = x, model = model, sigma_x = sigma_x))
  if(type == "taylor")
    return(pred_input_var_taylor(x = x, model = model, sigma_x = sigma_x))
  return(pred_input_var_simple(x = x, model = model, sigma_x = sigma_x))
}

#' Exact version for Gaussian kernel
#' @param x design considered
#' @param model GP
#' @param sigma_x input variance (scalar or vector)
#' @noRd
pred_input_var_exact <- function(x, model, sigma_x){
  d <- ncol(model$X0)
  n <- nrow(model$X0)
  # Sigma_x <- diag(sigma_x, d)
  W <- model$theta/2 # implicitly a diagonal matrix
  beta <- (model$Z0 - model$beta0) %*% model$Ki
  
  ## Mean part
  detWI <- sqrt(prod(sigma_x/model$theta * 2 + 1))
  
  q <- cov_gen(model$X0, x, theta = 2 * sigma_x + model$theta, type = "Gaussian")
  q <- q/detWI
  
  mean <- model$beta0 + drop(beta %*% q)
  
  ## Variance part
  detWI2 <- sqrt(prod(2 * sigma_x / W + 1))
  kx <- cov_gen(model$X0, x, theta = model$theta, type = "Gaussian")
  Q <- matrix(NA, n, n)
  Q2 <- Q3 <-  matrix(NA, n, n)
  
  for(i in 1:n){
    for(j in i:n){
      zx <- (model$X0[i,] + model$X0[j,])/2 - x
      Q[i,j] <- Q[j,i] <- exp(sum(zx^2/(W + 0.5 * W^2/sigma_x)))
    }
  }
  Q <- tcrossprod(kx)/detWI2 * Q
 
  if(model$trendtype != "SK") warning("Ordinary kriging variance not taken into account.")
  
  sd2 <- model$nu_hat * (1 + fast_trace(crossprod(beta)/model$nu_hat - model$Ki, Q)) - mean^2
  return(list(mean = mean, sd2 = sd2))
}

#' Version Girard 2003 for a new point
#' @param x design considered
#' @param model GP
#' @param sigma_x input variance (scalar or vector)
#' @noRd
pred_input_var_taylor <- function(x, model, sigma_x){
  preds <- predict(model, x)
  
  # for the variance
  p_gr <- predict_gr(model, x)$mean
  p_He <- predict_Hess(model, x)$sd2
  sd2 <- preds$sd2 + 0.5 * sum(diag(p_He)*sigma_x) + tcrossprod(p_gr * sigma_x, p_gr)
  
  return(list(mean = preds$mean, sd2 = sd2))
}


#' Version McHutchon simple
#' Simpler version developed in Section Formalization
#' @param x design considered
#' @param model GP
#' @param sigma_x input variance (scalar or vector)
#' @noRd
pred_input_var_simple <- function(x, model, sigma_x){
  preds <- predict(model, x)
  
  # For the variance
  p_gr <- predict_gr(model, x)$mean
  sd2 <- preds$sd2 + tcrossprod(p_gr * sigma_x, p_gr)
  return(list(mean = preds$mean, sd2 = sd2))
}


#' Hessian of the predictive mean and variance of hetGP models
#' @param object GP/TP model
#' @param x matrix of new designs, one point per row (size n x d)
#' @noRd
predict_Hess <- function(object, x){
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  kvec <-  cov_gen(X1 = object$X0, X2 = x, theta = object$theta, type = object$covtype)
  ZKi <-  crossprod(object$Z0 - object$beta0, object$Ki)
  d <- ncol(x)
  
  # dvec_num <- matrix(NA, nrow = length(kvec), ncol = ncol(x))
  # d2vec_num <- matrix(NA, nrow = length(kvec), ncol = ncol(x)^2)
  # for (i in 1:length(kvec)){
  #   dvec_num[i,] <- grad(cov_gen, object$X0[i,,drop = F], X2 = x, theta = object$theta, type = object$covtype)
  #   d2vec_num[i,] <- as.vector(hessian(cov_gen, object$X0[i,,drop = F], X2 = x, theta = object$theta, type = object$covtype))
  # }
  
  d2m <- d2s2 <- matrix(NA, ncol(x), ncol(x))
  
  for(i in 1:d){
    for(j in i:d){
      if(i == j){
        if(object$covtype == "Gaussian"){
          d2kvec <- 2 * (2 * (object$X0[,i] - x[i])^2 - object$theta[i]) / object$theta[i]^2 * kvec
          dkvec <- -2 * (object$X0[,i] - x[i]) / object$theta[i] * kvec
        }
        if(object$covtype == "Matern5_2"){
          tmp <- (object$X0[,i] - x[i])/object$theta[i]
          ids <- which(tmp < 0)
          tmp <- abs(tmp)
          dkvec <- (- 5/3 * tmp - 5 * sqrt(5.)/3. * tmp^2) / (1. + sqrt(5.) * tmp + 5./3. * tmp^2) * kvec / object$theta[i]
          dkvec[ids] <- - dkvec[ids]
          d2kvec <- (-5/3 - 5/3 * sqrt(5.) * tmp + 25/3. * tmp^2) / (1. + sqrt(5.) * tmp + 5./3. * tmp^2) * kvec / object$theta[i]^2
        }
        if(object$covtype == "Matern3_2"){
          tmp <- (object$X0[,i] - x[i])/object$theta[i]
          ids <- which(tmp < 0)
          tmp <- abs(tmp)
          dkvec <-  (- 3. * tmp) / (1. + sqrt(3.) * tmp) * kvec / object$theta[i]
          dkvec[ids] <- - dkvec[ids]
          d2kvec <- (3 * sqrt(3) * tmp - 3) / (1. + sqrt(3.) * tmp) * kvec / object$theta[i]^2
        }

        d2m[i, i] <- ZKi %*% d2kvec
        if(class(object) %in% c("hetGP", "homGP") && object$trendtype == "OK") 
          tmp <- -(colSums(object$Ki) %*% dkvec)^2 / sum(object$Ki) + drop(1 - colSums(object$Ki) %*% kvec)/(sum(object$Ki)) * colSums(object$Ki) %*% d2kvec else tmp <- 0
        d2s2[i, i] <- -  2 * (crossprod(dkvec, object$Ki) %*% dkvec + crossprod(kvec, object$Ki) %*% d2kvec + tmp)
      }else{
        if(object$covtype == "Gaussian"){
          d2kvec <- 4 * (object$X0[,i] - x[i]) / object$theta[i] * (object$X0[,j] - x[j]) / object$theta[j] * kvec
          dkveci <- -2 * (object$X0[,i] - x[i]) / object$theta[i] * kvec
          dkvecj <- -2 * (object$X0[,j] - x[j]) / object$theta[j] * kvec
        }
        if(object$covtype == "Matern5_2"){
          tmpi <- (object$X0[,i] - x[i]) / object$theta[i]
          tmpj <- (object$X0[,j] - x[j]) / object$theta[j]
          idsi <- which(tmpi < 0)
          idsj <- which(tmpj < 0)
          tmpi <- abs(tmpi)
          tmpj <- abs(tmpj)
          dkveci <-  (- 5/3 * tmpi - 5 * sqrt(5.)/3. * tmpi^2) / (1. + sqrt(5.) * tmpi + 5./3. * tmpi^2) / object$theta[i]
          dkvecj <-  (- 5/3 * tmpj - 5 * sqrt(5.)/3. * tmpj^2) / (1. + sqrt(5.) * tmpj + 5./3. * tmpj^2) / object$theta[j]
          dkveci[idsi] <- - dkveci[idsi]
          dkvecj[idsj] <- - dkvecj[idsj]
          d2kvec <- dkveci * dkvecj * kvec
          dkveci <- dkveci * kvec
          dkvecj <- dkvecj * kvec
        }
        if(object$covtype == "Matern3_2"){
          tmpi <- (object$X0[,i] - x[i]) / object$theta[i]
          tmpj <- (object$X0[,j] - x[j]) / object$theta[j]
          idsi <- which(tmpi < 0)
          idsj <- which(tmpj < 0)
          tmpi <- abs(tmpi)
          tmpj <- abs(tmpj)
          dkveci <-  (- 3. * tmpi)/(1. + sqrt(3.) * tmpi) / object$theta[i]
          dkvecj <-  (- 3. * tmpj)/(1. + sqrt(3.) * tmpj) / object$theta[j]
          dkveci[idsi] <- - dkveci[idsi]
          dkvecj[idsj] <- - dkvecj[idsj]
          d2kvec <- dkveci * dkvecj * kvec
          dkveci <- dkveci * kvec
          dkvecj <- dkvecj * kvec
        }

        d2m[i, j] <- d2m[j, i] <- ZKi %*% d2kvec
        if(class(object) %in% c("hetGP", "homGP") && object$trendtype == "OK") 
          tmp <- -(colSums(object$Ki) %*% dkveci) * (colSums(object$Ki) %*% dkvecj) / sum(object$Ki) + drop(1 - colSums(object$Ki) %*% kvec) / (sum(object$Ki)) * colSums(object$Ki) %*% d2kvec else tmp <- 0
        d2s2[i,j] <- d2s2[j,i] <- -  2 * (crossprod(dkveci, object$Ki) %*% dkvecj + crossprod(kvec, object$Ki) %*% d2kvec + tmp)
      }
    }
  }
  if(class(object) %in% c("hetGP", "homGP")){
    return(list(mean = d2m, sd2 = object$nu_hat * d2s2))
  }else{
    return(list(mean = object$sigma2 * d2m, sd2 =  (object$nu + object$psi - 2) / (object$nu + length(object$Z) - 2) * object$sigma2^2 * d2s2)) 
  }
}


