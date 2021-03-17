library(hetGP)
library(numDeriv)

context("Input noise")

test_that("Hessian 1D",{
  for(covt in c("Gaussian", "Matern3_2", "Matern5_2")){
    # 1D test
    set.seed(32)
    ## motorcycle data
    library(MASS)
    X <- matrix(mcycle$times, ncol = 1)/60
    Z <- mcycle$accel
    ## Model fitting
    mhet <- mleHetGP(X = X, Z = Z, lower = 0.05, upper = NULL, covtype = covt)
    mhom <- mleHomGP(X = X, Z = Z, lower = 0.05, upper = NULL, covtype = covt, known = list(beta0 = 0))
    
    x <- matrix(runif(1))
    
    pmf <- function(x, model) predict(model, x)$mean
    pmv <- function(x, model) predict(model, x)$sd2
 
    # pmk <- function(x, model) cov_gen(matrix(x, nrow = 1), model$X0[50,,drop = F], model$theta, type = model$covtype)
    # pmk1 <- function(x, model) cov_gen(matrix(x, nrow = 1), theta = model$theta)
    # grad(pmk, x, model = mhet)
    # hetGP:::partial_cov_gen(X1 = x[1,,drop = F], X2 = mhet$X0[50,,drop = F], theta = mhet$theta, i1 = 1, i2 = 1, arg = "X_i_j", type = mhet$covtype)*cov_gen(x, mhet$X0[50,,drop = F], mhet$theta, type = covt)
    
    # Gradient of prediction is also tested in optim tests
    expect_equal(c(mean = grad(func = pmf, x = x, model = mhet), sd2 = grad(func = pmv, x = x, model = mhet)),
                 unlist(hetGP:::predict_gr(mhet, x)), tol = 1e-4)
    expect_equal(c(mean = grad(func = pmf, x = x, model = mhom), sd2 = grad(func = pmv, x = x, model = mhom)),
                 unlist(hetGP:::predict_gr(mhom, x)), tol = 1e-4)
    
    # Now Hessian
    expect_equal(c(mean = hessian(func = pmf, x = x, model = mhet), sd2 = hessian(func = pmv, x = x, model = mhet)),
                 unlist(hetGP:::predict_Hess(mhet, x)), tol = 1e-2)
    
    expect_equal(c(mean = hessian(func = pmf, x = x, model = mhom), sd2 = hessian(func = pmv, x = x, model = mhom)),
                 unlist(hetGP:::predict_Hess(mhom, x)), tol = 1e-2)
  }
})


test_that("Hessian 2D",{
  for(covt in c("Gaussian", "Matern3_2", "Matern5_2")){
    # 1D test
    set.seed(32)
    nvar <- 2
    
    ## Branin redefined in [0,1]^2
    branin <- function(x){
      if(is.null(nrow(x)))
        x <- matrix(x, nrow = 1)
      x1 <- x[,1] * 15 - 5
      x2 <- x[,2] * 15
      (x2 - 5/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(x1) + 10
    }
    
    ## data generating function combining mean and noise fields
    ftest <- function(x){
      return(branin(x) + rnorm(nrow(x), mean = 0, sd = 2*rowSums(x^2)))
    }
    
    ## Grid of predictive locations
    ngrid <- 51
    xgrid <- matrix(seq(0, 1, length.out = ngrid), ncol = 1) 
    Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
    
    ## Unique (randomly chosen) design locations
    n <- 15
    Xu <- matrix(runif(n * 2), n)
    ## Select replication sites randomly
    X <- Xu[sample(1:n, 3*n, replace = TRUE),]
    ## obtain training data response at design locations X
    Z <- ftest(X)
    ## Formatting of data for model creation (find replicated observations) 
    prdata <- find_reps(X, Z, rescale = FALSE, normalize = FALSE)
    
    ## Model fitting
    mhet <- mleHetGP(X = list(X0 = prdata$X0, Z0 = prdata$Z0, mult = prdata$mult), Z = prdata$Z,
                      lower = rep(0.1, nvar), upper = rep(0.5, nvar), known = list(beta0 = 0),
                      covtype = covt)
    mhom <- mleHomGP(X = list(X0 = prdata$X0, Z0 = prdata$Z0, mult = prdata$mult), Z = prdata$Z,
                     lower = rep(0.1, nvar), upper = rep(0.5, nvar),
                     covtype = covt)
    
    x <- matrix(runif(nvar), 1)
    
    pmf <- function(x, model) predict(model, x)$mean
    pmv <- function(x, model) predict(model, x)$sd2
    
    # pmk <- function(x, model) cov_gen(matrix(x, nrow = 1), model$X0[50,,drop = F], model$theta, type = model$covtype)
    # pmk1 <- function(x, model) cov_gen(matrix(x, nrow = 1), theta = model$theta)
    # grad(pmk, x, model = mhet)
    # hetGP:::partial_cov_gen(X1 = x[1,,drop = F], X2 = mhet$X0[50,,drop = F], theta = mhet$theta, i1 = 1, i2 = 1, arg = "X_i_j", type = mhet$covtype)*cov_gen(x, mhet$X0[50,,drop = F], mhet$theta, type = covt)
    
    # Gradient of prediction is also tested in optim tests
    expect_equal(c(mean = grad(func = pmf, x = x, model = mhet), sd2 = grad(func = pmv, x = x, model = mhet)),
                 unlist(hetGP:::predict_gr(mhet, x)), tol = 1e-4)
    expect_equal(c(mean = grad(func = pmf, x = x, model = mhom), sd2 = grad(func = pmv, x = x, model = mhom)),
                 unlist(hetGP:::predict_gr(mhom, x)), tol = 1e-4)
    
    # Now Hessian
    expect_equal(c(mean = hessian(func = pmf, x = x, model = mhet), sd2 = hessian(func = pmv, x = x, model = mhet)),
                 unlist(hetGP:::predict_Hess(mhet, x)), tol = 1e-2)
    
    expect_equal(c(mean = hessian(func = pmf, x = x, model = mhom), sd2 = hessian(func = pmv, x = x, model = mhom)),
                 unlist(hetGP:::predict_Hess(mhom, x)), tol = 1e-2)
  }
})



