## This is the main functions for spatial envelope
# This functions are prepared by Z. Naji and H. Moradi Rekabdarkolaee

env.fun<-function (X, Y, u,s,cov.form="matern",init.theta=c(1,0.5),init.error=0.0005) 
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  if(cov.form!="matern"){
    if(length(init.theta)>1)
      warning("init.theta should be an scalar for exponential covariance. We set smoothness = 0.5.")
    init.theta<-init.theta[1]
  }
  if (a[1] != nrow(X)) 
    stop("X and Y should have the same number of observations.")
  if (u > r || u < 0) 
    stop("u must be an integer between 0 and r.")
  if (sum(duplicated(cbind(X, Y), MARGIN = 2)) > 0) 
    stop("Some responses also appear in the predictors, or there maybe duplicated columns in X or Y.")
  Ybar<-apply(Y,2,mean)
  Xbar<-apply(X,2,mean)
  one<-rep(1,n)
  H<-Y-kronecker(t(Ybar),one)
  G<-X-kronecker(t(Xbar),one)
  init.value<-starting.value(X, Y, s, u, cov.form,init.theta)
  
  beta0<-init.value$beta
  alpha0<-init.value$alpha
  Omega00<-init.value$Omega0
  Omega10<-init.value$Omega
  Gamma10<-init.value$Gamma1
  Gamma00<-init.value$Gamma0
  
  init.theta=init.value$theta
  Dis<-as.matrix(dist(s))
  e<-init.error+1
  
  maxiter <- 100
  iter <- 0
  ###########
  while(e>init.error&& iter < maxiter){
    iter = iter +1
    if(cov.form=="matern"){
      rho<-Matern(Dis , range = init.theta[1], smoothness =init.theta[2])
    }else{
      rho<-Matern(Dis , range =init.theta)   
    }
    Bmle<-ginv(t(G)%*%ginv(rho)%*%G)%*%t(G)%*%ginv(rho)%*%H
    
    sigY<-t(H)%*%ginv(rho)%*%H
    U<-t(H)%*%ginv(rho)%*%G%*%ginv(t(G)%*%ginv(rho)%*%G)%*%t(G)%*%ginv(rho)%*%H
    M<-sigY-U
    tmp <- envMU(M, U, u)
    
    Gammahat <- tmp$Gammahat
    Gamma0hat <- tmp$Gamma0hat
    objfun <- tmp$objfun
    covMatrix <- NULL
    asySE <- NULL
    ratio <- NULL
    if (u == 0) {
      etahat <- NULL
      Omegahat <- NULL
      Omega0hat <- sigY
      alphahat  <- colMeans(Y)
      betahat <- matrix(0, r, p)
      Sigmahat <- sigY
      Theta<-Theta.fun(X, Y,u, s,Gammahat,
                       Gamma0hat,Omegahat,Omega0hat,cov.form)
      thetahat<-optim(init.theta, fn=Theta, lower = c(0.1,0.1), upper = c(7,5), method = "L-BFGS-B")$par
      
    }
    else if (u == r) {
      etahat <- Bmle
      Omegahat <- M
      Omega0hat <- NULL
      alphahat  <- colMeans(Y) - Bmle %*% colMeans(X)
      betahat <- Bmle
      Sigmahat <- M
      Theta<-Theta.fun(X, Y,u, s,Gammahat,
                       Gamma0hat,Omegahat,Omega0hat,cov.form)
      thetahat<-optim(init.theta, fn=Theta, lower = c(0.1,0.1), upper = c(7,5), method = "L-BFGS-B")$par
      
    }
    else {
      Om <- crossprod(Gammahat, M) %*% Gammahat
      Om0 <- crossprod(Gamma0hat, sigY) %*% Gamma0hat
      Sigma1 <- Gammahat %*% tcrossprod(Om, Gammahat)
      Sigmahat <- Sigma1 + Gamma0hat %*% tcrossprod(Om0, 
                                                    Gamma0hat)
      Omegahat<-Om
      Omega0hat<-Om0
      
      etahat <- crossprod(Gammahat, t(Bmle))
      betahat <- Gammahat %*% etahat
      
      alphahat  <- colMeans(Y) - betahat %*% colMeans(X)
      
      Theta<-Theta.fun(X, Y,u, s,Gammahat,
                       Gamma0hat,Omegahat,Omega0hat,cov.form)
      thetahat<-optim(init.theta, fn=Theta, lower = c(0.1,0.1), upper = c(7,5), method = "L-BFGS-B")$par
      
    }
    e<-mean((betahat-beta0)^2)+mean((alphahat-alpha0)^2)+
      mean((Gammahat-Gamma10)^2)+mean((Gamma0hat-Gamma00)^2)
    
    init.theta<-thetahat
    beta0<-betahat
    alpha0<-alphahat
    Omega00<-Omega0hat
    Omega10<-Omegahat
    Gamma10<-Gammahat
    Gamma00<-Gamma0hat
    
    
    
  }
  
  
  return(list(Gamma1 = Gammahat, Gamma0 = Gamma0hat, alpha = alphahat, 
              beta = betahat, Sigma = Sigmahat, eta = etahat, Omega = Omegahat, 
              Omega0 = Omega0hat, theta=thetahat))
}
