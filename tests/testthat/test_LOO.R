library(hetGP)
context("LOO")


test_that("LOO",{
  
  set.seed(32)
  ## motorcycle data
  library(MASS)
  X <- matrix(mcycle$times, ncol = 1)
  Z <- mcycle$accel
  nvar <- 1
  
  ## Start with GP models
  for(modelfun in c("mleHomGP", "mleHetGP")){
    for(trend in c(NA, 0)){
      
      if(is.na(trend)) trend <- NULL
      ## Model fitting
      model <- match.fun(modelfun)(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(10, nvar),
                                   covtype = "Matern5_2", known = list(beta0 = trend))
      LOO_p <- LOO_preds(model)
      
      # model minus observation(s) at x_i
      d_mot <- find_reps(X, Z)
      
      LOO_ref <- matrix(NA, nrow(d_mot$X0), 2)
      for(i in 1:nrow(d_mot$X0)){
        model_i <- match.fun(modelfun)(X = list(X0 = d_mot$X0[-i,, drop = FALSE], Z0 = d_mot$Z0[-i],
                                                mult = d_mot$mult[-i]), Z = unlist(d_mot$Zlist[-i]),
                                       lower = rep(0.1, nvar), upper = rep(50, nvar), covtype = "Matern5_2",
                                       known = list(theta = model$theta, k_theta_g = model$k_theta_g, g = model$g,
                                                    Delta = model$Delta[-i], beta0 = trend))
        model_i$nu_hat <- model$nu_hat
        
        # For hetGP, need to use the same Lambdas to get the same results 
        if(modelfun == "mleHetGP"){
          model_i$Lambda <- model$Lambda[-i]
          model_i <- strip(model_i)
          model_i <- rebuild(model_i)
        }   
        
        p_i <- predict(model_i, d_mot$X0[i,,drop = FALSE])
        LOO_ref[i,] <- c(p_i$mean, p_i$sd2)
      }
      # Compare results
      expect_equal(LOO_ref[,1], as.numeric(LOO_p$mean), tol = 1e-6)
      expect_equal(LOO_ref[,2], as.numeric(LOO_p$sd2), tol = 1e-5)
    }
  }
  
  ## Then TP versions
  for(modelfun in c("mleHomTP", "mleHetTP")){
    
    ## Model fitting
    model <- match.fun(modelfun)(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(10, nvar),
                                 covtype = "Matern5_2", known = list(beta0 = 0))
    LOO_p <- LOO_preds(model)
    
    # model minus observation(s) at x_i
    d_mot <- find_reps(X, Z)
    
    LOO_ref <- matrix(NA, nrow(d_mot$X0), 2)
    for(i in 1:nrow(d_mot$X0)){
      model_i <- match.fun(modelfun)(X = list(X0 = d_mot$X0[-i,, drop = FALSE], Z0 = d_mot$Z0[-i],
                                              mult = d_mot$mult[-i]), Z = unlist(d_mot$Zlist[-i]),
                                     lower = rep(0.1, nvar), upper = rep(50, nvar), covtype = "Matern5_2",
                                     known = list(theta = model$theta, k_theta_g = model$k_theta_g, g = model$g,
                                                  sigma2 = model$sigma2, nu = model$nu,
                                                  Delta = model$Delta[-i], beta0 = trend))
      
      model_i$psi <- model$psi # psi is taken as fixed
      # For hetTP, need to use the same Lambdas and psi to get the same results 
      if(modelfun == "mleHetTP"){
        model_i$Lambda <- model$Lambda[-i]
        model_i <- strip(model_i)
        model_i <- rebuild(model_i)
      }   
      
      p_i <- predict(model_i, d_mot$X0[i,,drop = FALSE])
      LOO_ref[i,] <- c(p_i$mean, p_i$sd2)
    }
    # Compare results
    expect_equal(LOO_ref[,1], as.numeric(LOO_p$mean), tol = 1e-6)
    expect_equal(LOO_ref[,2], as.numeric(LOO_p$sd2), tol = 1e-6)
  }
})


test_that("LOO error",{
  set.seed(32)
  ## motorcycle data
  library(MASS)
  X <- matrix(mcycle$times, ncol = 1)
  Z <- mcycle$accel
  nvar <- 1
  
  # X <- model$X0
  # Z <- model$Z0
  
  for(trend in c(0)){
    # Only for simple kriging for now
    if(is.na(trend)) trend <- NULL
    
    model <- mleHomGP(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(10, nvar),
                      covtype = "Matern5_2", known = list(beta0 = trend))
    
    if(max(model$mult) == 1){
      LOO_p <- LOO_preds(model)
      LOO_p$sd2 <- LOO_p$sd2/model$nu_hat + model$g
      err_p <- mean((rep(LOO_p$mean, model$mult) - model$Z)^2/LOO_p$sd2 + log(LOO_p$sd2))
    }else{
      err_p <- NULL
      for(i in 1:nrow(X)){
        mtmp <- mleHomGP(X = X[-i,], Z = Z[-i], known = list(theta = model$theta, g = model$g, beta0 = trend),
                         covtype = "Matern5_2")
        mtmp$nu_hat <- model$nu_hat
        p_i <- predict(mtmp, X[i,,drop = F])
        p_i$sd2 <- p_i$sd2/model$nu_hat + model$g
        err_p <- c(err_p, (p_i$mean - Z[i])^2/p_i$sd2 + log(p_i$sd2))
      }
      err_p <- mean(err_p)
    }
    # err_alt <- drop(t(Z) %*% model$Ki %*% diag(diag(model$Ki)^(-2)) %*% diag(1/(diag(model$Ki)^(-1) - model$g)) %*% model$Ki %*% Z / model$nu_hat)/length(Z)
    err_loofun <- hetGP:::LOO_hom(model$X0, model$Z0, model$Z, model$mult, model$theta, model$g, beta0 = trend, covtype = "Matern5_2")
    expect_equal(err_p, err_loofun, tolerance = 1e-6)
    # expect_equal(err_alt, err_loofun, tolerance = 1e-6)
  }
})

test_that("LOO gradient",{
  library(numDeriv)
  set.seed(32)
  
  d <- 2
  
  ## Select replication sites randomly
  n <- 4
  Xu <- matrix(runif(d * n), n)
  mult <- c(1, 2, 3, 2)
  X <- Xu[sample(1:n, 2*n, replace = TRUE),]
  Z <- runif(sum(mult))

  # funlooN <- function(X0, Z, theta, mult, g, covtype = "Gaussian"){
  #   n <- nrow(X0)
  #   XN <- X0[rep(1:n, times = mult),]
  #   SN <- rep(g, nrow(XN)) # S0[rep(1:n, times = mult)]
  #   KN <- hetGP::cov_gen(XN, theta = theta, type = covtype)
  #   KNi <- solve(KN + diag(SN))
  #   
  #   return(t(Z) %*% KNi %*% diag(diag(KNi)^(-2)) %*% KNi %*% Z)
  # }
  
  for(trend in c(0)){
    # Only for simple kriging for now
    if(is.na(trend)) trend <- NULL
    
    model <- mleHomGP(X = X, Z = Z, lower = rep(0.1, d), upper = rep(10, d),
                      covtype = "Matern5_2", known = list(beta0 = trend, g = 1.3e-1))
    
    err_loofun <- hetGP:::LOO_hom(X0 = model$X0, Z0 = model$Z0, Z = model$Z, mult = model$mult, theta = model$theta, g = model$g, beta0 = trend, covtype = "Matern5_2")
    err_dloo <- hetGP:::dLOOHom(X0 = model$X0, Z0 = model$Z0, Z = model$Z, mult = model$mult, theta = model$theta, g = model$g, beta0 = trend, covtype = "Matern5_2")
    grads <- c(grad(hetGP:::LOO_hom, X0 = model$X0, Z0 = model$Z0, Z = model$Z, mult = model$mult, x = model$theta, g = model$g, beta0 = trend, covtype = "Matern5_2"),
               grad(hetGP:::LOO_hom, X0 = model$X0, Z0 = model$Z0, Z = model$Z, mult = model$mult, theta = model$theta, x = model$g, beta0 = trend, covtype = "Matern5_2"))
    expect_equal(err_dloo, grads)
    
  }
})

