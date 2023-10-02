library(hetGP)
context("Covariance derivatives")


test_that("Deriv",{
  library(hetGP)
  d <- 3
  
  for(type in c('Gaussian', 'Matern5_2', 'Matern3_2')){
    X <- matrix(runif(4*d), 4, d)
    theta <- runif(d) + 0.2
    
    # Covariance matrix case
    K <- cov_gen(X, theta= theta, type = type)
    
    expect_equal(hetGP:::partial_cov_gen(X1 = X[,1,drop = F], theta = theta[1], type = type, arg = "theta_k")*K,
                 1e6 * (cov_gen(X1 = X, theta = theta + c(1e-6, 0, 0), type = type) - cov_gen(X1 = X, theta = theta, type = type)), tol = 1e-4)
    
    expect_equal(hetGP:::partial_cov_gen(X1 = X[,2,drop = F], theta = theta[2], type = type, arg = "theta_k")*K,
                 1e6 * (cov_gen(X1 = X, theta = theta + c(0, 1e-6, 0), type = type) - cov_gen(X1 = X, theta = theta, type = type)), tol = 1e-4)
    
    # expect_equal(hetGP:::partial_cov_gen(X1 = X[,3,drop = F], theta = theta[3], type = type, arg = "theta_k")*K,
    #              1e6 * (cov_gen(X1 = X, theta = theta + c(0, 0, 1e-6), type = type) - cov_gen(X1 = X, theta = theta, type = type)), tolerance = 1e-4)
    
    # i <- sample(1:4, 1)
    # j <- sample(1:3, 1)
    # eMat <- matrix(0, nrow(X), ncol(X))
    # eMat[i, j] <- 1e-6
    # expect_equal(drop(hetGP:::partial_cov_gen(X1 = X[i,j,drop = F], theta = theta, type = type, arg = "X_i_j", i1 = i, i2 = j))*K,
    #              1e6 * (cov_gen(X1 = X + eMat, theta = theta, type = type) - cov_gen(X1 = X, theta = theta, type = type)), tolerance = 1e-4)
    
    # cross matrix case
    Y <- matrix(runif(6*d), 6, d)
    K <- cov_gen(X, Y, theta= theta, type = type)
    
    expect_equal(hetGP:::partial_cov_gen(X1 = X[,1,drop = F], X2 = Y[,1,drop = F], theta = theta[1], type = type, arg = "theta_k")*K,
                 1e6 * (cov_gen(X1 = X, X2 = Y, theta = theta + c(1e-6, 0, 0), type = type) - cov_gen(X1 = X, X2 = Y, theta = theta, type = type)), tol = 1e-4)
    
    expect_equal(hetGP:::partial_cov_gen(X1 = X[,2,drop = F], X2 = Y[,2,drop = F], theta = theta[2], type = type, arg = "theta_k")*K,
                 1e6 * (cov_gen(X1 = X, X2 = Y, theta = theta + c(0, 1e-6, 0), type = type) - cov_gen(X1 = X, X2 = Y, theta = theta, type = type)), tol = 1e-4)
    
    expect_equal(hetGP:::partial_cov_gen(X1 = X[,3,drop = F], X2 = Y[,3,drop = F], theta = theta[3], type = type, arg = "theta_k")*K,
                 1e6 * (cov_gen(X1 = X, X2 = Y, theta = theta + c(0, 0, 1e-6), type = type) - cov_gen(X1 = X, X2 = Y, theta = theta, type = type)), tol = 1e-4)
    
    # expect_equal(drop(hetGP:::partial_cov_gen(X1 = X[i,j,drop = F], X2 = Y, theta = theta, type = type, arg = "X_i_j", i1 = i, i2 = j))*K,
    #              1e6 * (cov_gen(X1 = X + eMat, X2 = Y, theta = theta, type = type) - cov_gen(X1 = X, X2 = Y, theta = theta, type = type)), tolerance = 1e-4)
    
  }
  
})
