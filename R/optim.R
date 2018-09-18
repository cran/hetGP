##' Computes EI for minimization
##' @title Expected Improvement criterion
##' @param x matrix of new designs, one point per row (size n x d)
##' @param model \code{homGP} or \code{hetGP} model, or their TP equivalents, including inverse matrices
##' @param cst optional plugin value used in the EI, see details
##' @param preds optional predictions at \code{x} to avoid recomputing if already done
##' @note This is a beta version at this point.
##' @details 
##' \code{cst} is classically the observed minimum in the determininistic case. 
##' In the noisy case, the min of the predictive mean works fine.
##'  
##' @references
##' Mockus, J.; Tiesis, V. & Zilinskas, A. (1978).
##' The application of Bayesian methods for seeking the extremum Towards Global Optimization, Amsterdam: Elsevier, 2, 2.\cr 
##' 
##' A. Shah, A. Wilson, Z. Ghahramani (2014), Student-t processes as alternatives to Gaussian processes, Artificial Intelligence and Statistics, 877--885.
## ' # FIXME: reference for use of mean as plugin    
##' @importFrom stats dt qt
##' @export
##' @examples 
##' ## Optimization example
##' set.seed(42)
##' 
##' 
##' ## Noise field via standard deviation
##' noiseFun <- function(x, coef = 1.1, scale = 1){
##' if(is.null(nrow(x)))
##'  x <- matrix(x, nrow = 1)
##'    return(scale*(coef + cos(x * 2 * pi)))
##' }
##' 
##' ## Test function defined in [0,1]
##' ftest <- function(x){
##' if(is.null(nrow(x)))
##' x <- matrix(x, ncol = 1)
##' return(f1d(x) + rnorm(nrow(x), mean = 0, sd = noiseFun(x)))
##' }
##' 
##' n_init <- 10 # number of unique designs
##' N_init <- 100 # total number of points
##' X <- seq(0, 1, length.out = 10)
##' X <- matrix(X[sample(1:n_init, N_init, replace = TRUE)], ncol = 1)
##' Z <- ftest(X)
##' 
##' ## Predictive grid
##' ngrid <- 51
##' xgrid <- seq(0,1, length.out = ngrid)
##' Xgrid <- matrix(xgrid, ncol = 1)
##' 
##' model <- mleHetGP(X = X, Z = Z, lower = 0.001, upper = 1)
##' 
##' EIgrid <- crit_EI(Xgrid, model)
##' preds <- predict(x = Xgrid, model)
##' 
##' par(mar = c(3,3,2,3)+0.1)
##' plot(xgrid, f1d(xgrid), type = 'l', lwd = 1, col = "blue", lty = 3,
##' xlab = '', ylab = '', ylim = c(-8,16))
##' points(X, Z)
##' lines(Xgrid, preds$mean, col = 'red', lwd = 2)
##' lines(Xgrid, qnorm(0.05, preds$mean, sqrt(preds$sd2)), col = 2, lty = 2)
##' lines(Xgrid, qnorm(0.95, preds$mean, sqrt(preds$sd2)), col = 2, lty = 2)
##' lines(Xgrid, qnorm(0.05, preds$mean, sqrt(preds$sd2 + preds$nugs)), col = 3, lty = 2)
##' lines(Xgrid, qnorm(0.95, preds$mean, sqrt(preds$sd2 + preds$nugs)), col = 3, lty = 2)
##' par(new = TRUE)
##' plot(NA, NA, xlim = c(0, 1), ylim = c(0,max(EIgrid)), axes = FALSE, ylab = "", xlab = "")
##' lines(xgrid, EIgrid, lwd = 2, col = 'cyan')
##' axis(side = 4)
##' mtext(side = 4, line = 2, expression(EI(x)), cex = 0.8)
##' mtext(side = 2, line = 2, expression(f(x)), cex = 0.8)
crit_EI <- function(x, model, cst = NULL, preds = NULL){
  if(is.null(cst)) cst <- min(predict(model, x = model$X)$mean)
  if(is.null(dim(x))) x <- matrix(x, nrow = 1)
  if(is.null(preds)) preds <- predict(model, x = x)
  
  if(class(model) %in% c("homTP", "hetTP")){
    gamma <- (cst - preds$mean)/sqrt(preds$sd2)
    res <- (cst - preds$mean) * pt(gamma, df = model$nu + length(model$Z))
    res <- res + sqrt(preds$sd2) * (1 + (gamma^2 - 1)/(model$nu + length(model$Z) - 1)) * dt(x = gamma, df = model$nu + length(model$Z))
    return(res)
  }
  
  res <- (cst - preds$mean) * pnorm(cst, mean = preds$mean, sd = sqrt(preds$sd2))
  res <- res + sqrt(preds$sd2) * dnorm(cst, mean = preds$mean, sd = sqrt(preds$sd2))
  return(res) 
}



