library(hetGP)
context("CRN")


test_that("Lik and dlik",{
  library(hetGP)
  
  ##------------------------------------------------------------
  ## Example 1: Homoskedastic GP modeling on 2d sims
  ##------------------------------------------------------------
  
  set.seed(2)
  nx <- 21
  ns <- 5
  d <- 2
  x <- as.matrix(expand.grid(seq(0,1, length.out = nx), seq(0,1, length.out = nx)))
  s <- matrix(seq(1, ns, length.out = ns))
  Xgrid <- as.matrix(expand.grid(seq(1, ns, length.out = ns), seq(0,1, length.out = nx), 
                                 seq(0,1, length.out = nx)))
  Xgrid <- Xgrid[,c(2, 3, 1)]
  
  covtype <- "Gaussian"
  g <- 1e-3
  theta <- c(0.01, 0.05)
  KX <- cov_gen(x, theta = theta, type = covtype)
  rho <- 0.33
  KS <- matrix(rho, ns, ns)
  diag(KS) <- 1
  Km <- kronecker(KX, KS) + g * diag(nx*nx*ns)
  YY <- MASS::mvrnorm(n = 1, mu = rep(0, nx*nx*ns), Sigma = Km)
  YYmat <- matrix(YY, ns, nx*nx)
  filled.contour(matrix(YYmat[1,], nx))
  filled.contour(matrix(YYmat[2,], nx))
  
  n0 <- 40
  ids <- sample(1:nrow(Xgrid), n0)
  X0 <- Xgrid[ids,]
  Y0 <-  YY[ids]
  
  # Check likelihood
  modelref <- mleCRNGP(X0, Y0, known = list(theta = theta, g = g, rho = rho, beta0 = 0), covtype = covtype)
  expect_equal(modelref$ll, mvtnorm::dmvnorm(x = Y0, mean = rep(0, n0), sigma = modelref$nu_hat*chol2inv(chol(modelref$Ki)), log = TRUE))
  
  # Check likelihood derivatives
  l0 <- hetGP:::logLikHomCRN(X0 = X0[,1:d, drop = F], S0 = X0[,d+1,drop = F], Z = Y0, theta = theta, g = g, rho = rho, beta0 = 0, covtype = covtype)
  
  expect_equal(hetGP:::dlogLikHomCRN(X0 = X0[,1:d, drop = F], S0 = X0[,d+1,drop = F], Y0, theta, g, rho, beta0 = 0, covtype = covtype),
               c(1e6 * (hetGP:::logLikHomCRN(X0 = X0[,1:d, drop = F], S0 = X0[,d+1,drop = F], Y0, theta + c(1e-6,0), g, rho, beta0 = 0, covtype = covtype) - l0),
                 1e6 * (hetGP:::logLikHomCRN(X0 = X0[,1:d, drop = F], S0 = X0[,d+1,drop = F], Y0, theta + c(0,1e-6), g, rho, beta0 = 0, covtype = covtype) - l0),
                 1e6 * (hetGP:::logLikHomCRN(X0 = X0[,1:d, drop = F], S0 = X0[,d+1,drop = F], Y0, theta, g + 1e-6, rho, beta0 = 0, covtype = covtype) - l0),
                 1e6 * (hetGP:::logLikHomCRN(X0 = X0[,1:d, drop = F], S0 = X0[,d+1,drop = F], Y0, theta, g, rho + 1e-6, beta0 = 0, covtype = covtype) - l0)), tolerance = 1e-4)
  
  ## Check consistency of the various formulations (knonecker or no kronecker)
  
  
})

test_that("time treatment",{
  library(hetGP)
  
  set.seed(42)
  nx <- 7
  nt <- 5
  ns <- 6
  x <- matrix(sort(seq(0,1, length.out = nx)), nx)
  s <- matrix(sort(seq(1, ns, length.out = ns)))
  t <- matrix(sort(seq(0, 1, length.out = nt)), nt)
  g <- 1e-3
  theta <- c(0.1, 0.2)
  covtype <- "Matern5_2"
  KX <- cov_gen(x, theta = theta[1], type = covtype)
  KT <- cov_gen(t, theta = theta[2], type = covtype)
  rho <- 0.3
  KS <- matrix(rho, ns, ns)
  diag(KS) <- 1
  XST <- as.matrix(expand.grid(x, s, t))
  KM <- kronecker(KT, kronecker(KS, KX))
  
  KM2 <- matrix(NA, nrow(XST), nrow(XST))
  for(i in 1:nrow(XST)){
    for(j in i:nrow(XST)){
      if(XST[i, 2] == XST[j, 2]) tmp <- 1 else tmp <- rho
      KM2[i, j] <- cov_gen(XST[i,c(1,3), drop = F], XST[j,c(1,3), drop = F], theta = theta, type = covtype) * tmp
      KM2[j, i] <- KM2[i, j]
    }
  }
  expect_equal(KM, KM2)
  
  ## YY <- MASS::mvrnorm(n = 1, mu = rep(0, nx*nt*ns), Sigma = KM + g * diag(nx*nt*ns))
  
  Kmc <- kronecker(chol(KT), kronecker(chol(KS), chol(KX)))
  YY <- t(Kmc) %*% rnorm(nrow(Kmc))
  
  ninit <- 40
  XS <- as.matrix(expand.grid(x, s))
  ids <- sort(sample(1:nrow(XS), ninit))
  XST0 <- cbind(XS[ids[rep(1:ninit, each = nt)],], rep(t[,1], times = ninit))
  X0 <- XST[which(duplicated(rbind(XST, XST0), fromLast = TRUE)),]
  Y0 <-  YY[which(duplicated(rbind(XST, XST0), fromLast = TRUE))]
  
  tmp <- hetGP:::find_reps(X = X0[,-3], Y0)
  
  expect_equal(tmp$Zlist[[13]], matrix(Y0, ncol = nt)[13,])
  
  model <- mleCRNGP(X = XS[ids,], T0=t, Z = matrix(Y0, ncol = nt), 
                    known = list(theta = theta, g = g), covtype = covtype)
  model_old <- hetGP:::mleCRNGPT(X0, Y0, known = list(theta = theta, g = g), covtype = covtype)
  
  expect_equal(model$ll, model_old$ll)
  
  preds <- predict(model, x = XS, xprime = XS)
  pold <- hetGP:::predict.CRNGPT(model_old, x = XST, xprime = XST)
  
  # compare with regular CRN GP
  mref <- mleCRNGP(X = X0[, c(1, 3, 2)], Z = Y0, known = list(theta = theta, g = g), covtype = covtype)#, known = list(theta = model$theta, g = model$g))
  pref <- predict(mref, x = XST[, c(1, 3, 2)], xprime = XST[, c(1, 3, 2)])
  
  expect_equal(as.vector(preds$mean), pref$mean, tol = 1e-6)
  expect_equal(preds$cov, pref$cov, tol = 1e-6)
  expect_equal(as.vector(preds$sd2), pref$sd2, tol = 1e-6)
  
  expect_equal(as.vector(preds$mean), pold$mean, tol = 1e-6)
  expect_equal(preds$cov, pold$cov, tol = 1e-6)
  expect_equal(as.vector(preds$sd2), pold$sd2, tol = 1e-6)
  
  # Check likelihood derivatives
  l0 <- hetGP:::logLikHomCRNT(X0 = XS[ids,1,drop = F], T0 = t, S0 = XS[ids,2,drop = F], Z = Y0, theta = theta, g = g, rho = rho, beta0 = 0, covtype = covtype)
  
  expect_equal(hetGP:::dlogLikHomCRNT(X0 = XS[ids,1,drop = F], T0 = t, S0 = XS[ids,2,drop = F], Y0, theta, g, rho, beta0 = 0, covtype = covtype),
               c(1e6 * (hetGP:::logLikHomCRNT(X0 = XS[ids,1,drop = F], T0 = t, S0 = XS[ids,2,drop = F], Y0, theta + c(1e-6,0), g, rho, beta0 = 0, covtype = covtype) - l0),
                 1e6 * (hetGP:::logLikHomCRNT(X0 = XS[ids,1,drop = F], T0 = t, S0 = XS[ids,2,drop = F], Y0, theta + c(0,1e-6), g, rho, beta0 = 0, covtype = covtype) - l0),
                 1e6 * (hetGP:::logLikHomCRNT(X0 = XS[ids,1,drop = F], T0 = t, S0 = XS[ids,2,drop = F], Y0, theta, g + 1e-6, rho, beta0 = 0, covtype = covtype) - l0),
                 1e6 * (hetGP:::logLikHomCRNT(X0 = XS[ids,1,drop = F], T0 = t, S0 = XS[ids,2,drop = F], Y0, theta, g, rho + 1e-6, beta0 = 0, covtype = covtype) - l0)), tolerance = 5e-4)
  
  
})


test_that("simul", {
  # # To be completed
  # nx <- 21
  # ns <- 5
  # d <- 2
  # x <- as.matrix(expand.grid(seq(0,1, length.out = nx), seq(0,1, length.out = nx)))
  # s <- matrix(seq(1, ns, length.out = ns))
  # Xgrid <- as.matrix(expand.grid(seq(1, ns, length.out = ns), seq(0,1, length.out = nx),
  #                                seq(0,1, length.out = nx)))
  # Xgrid <- Xgrid[,c(2, 3, 1)]
  # g <- 1e-3
  # theta <- c(0.02, 0.05)
  # KX <- cov_gen(x, theta = theta)
  # rho <- 0.33
  # KS <- matrix(rho, ns, ns)
  # diag(KS) <- 1
  # 
  # YY <- MASS::mvrnorm(n = 1, mu = rep(0, nx*nx*ns), Sigma = kronecker(KX, KS) + g * diag(nx*nx*ns))
  # YYmat <- matrix(YY, ns, nx*nx)
  # filled.contour(matrix(YYmat[1,], nx))
  # filled.contour(matrix(YYmat[2,], nx))
  # 
  # ids <- sample(1:nrow(Xgrid), 80)
  # X0 <- Xgrid[ids,]
  # Y0 <-  YY[ids]
  # 
  # model <- mleCRNGP(X0, Y0, know = list(beta0 = 0)) #, known = list(theta = 0.01, g = 1e-3, rho = 0.3))
  # 
  # 
  # kx <- model$nu_hat * cov_gen(X1 = Xgrid[, -(d+1), drop = FALSE],
  #                              X2 = model$X0, theta = model$theta, type = model$covtype)
  # tmp <- outer(Xgrid[, d + 1], model$S0, "==")
  # ks <- matrix(model$rho, nrow(Xgrid), nrow(model$X0))
  # ks[tmp] <- 1
  # kx <- kx * ks
  # kw <- kx %*% model$Ki/model$nu_hat # Kriging weights
  # 
  # kxxprime <- model$nu_hat * cov_gen(X1 = Xgrid[,-(d+1)], theta = model$theta, type = model$covtype)
  # ksxxprime <- matrix(model$rho, ns, ns)
  # tmp <- outer(unique(Xgrid[, d + 1]), unique(Xgrid[, d + 1]), "==")
  # ksxxprime[tmp] <- 1
  # 
  # cholK <- kronecker(chol(kxxprime + diag(1e-8, nrow(Xgrid))), chol(ksxxprime))
  # cK2 <- chol(kronecker(kxxprime, ksxxprime) + diag(1e-8, nrow(Xgrid) * ns))
  # expect_equal(cholK, cK2)
  # 
  # uncondreas <- t(cholK) %*% rnorm(nrow(Xgrid0)*ns)
  # mean <- as.vector(model$beta0 + kw %*% (model$Z - model$beta0))
  # altmean <- as.vector(model$beta0 + kw %*% (uncondreas[ids] - model$beta0)) # maybe it is not ids (?)
  # resi <- uncondreas - altmean
  # condreas <- mean + resi
  # # plot3d(X0[,1], X0[,2], Y0, size = 10, col = 1 + ((X0[,3] - 1) %% 6))
  # # surface3d(unique(x[,1]), unique(x[,2]), matrix(condreas[Xgrid[,3] == 1], nx), col = 1,
  # #   front = "lines", back = "lines")
})

