## This Function Generates the M-fold Cross validation for 
## Multivariate Regression, Envelope, and Spatial Envelope
## 
## The inputs for this function are
## X: denotes the matrix of covariates
## Y: denotes the matrix of response
## u: denotes the structural dimension of the envelope
## s: denotes the location
## m: denotes the numebr of fold for cross-validation
## nprem: denotes the number of permutations for data
## Cov.form: denotes the form of spatial covariance function
## init.theta: denotes the initial values for spatial covariance dunction
## init.error: denotes the initial error value

cv.spenv<-function (X, Y, u,s, m, nperm,cov.form="matern",init.theta=c(1,0.5),init.error=0.0005) 
{
      n<-nrow(Y)
      r<-ncol(Y)
      D<-diag(r)
      Dis<-as.matrix(dist(s))
      theta=init.theta
      if(cov.form=="matern"){
            rho<-Matern(Dis , range = theta[1], smoothness =theta[2])
      }else{
            rho<-Matern(Dis , range =theta)
      }
      p <- ncol(X)
      prederr <- rep(0, nperm)
      prederr.sp<- rep(0, nperm)
      prederr.ols<-rep(0, nperm)
      prederr.sp.CN<-rep(0, nperm)
      for (i in 1:nperm) {
            id <- sample(n, n)
            Xn <- as.matrix(X[id, ])
            Yn <- Y[id, ]
            sn<-s[id, ]
            Disn<-as.matrix(dist(sn))
            for (j in 1:m) {
                  id.test <- (floor((j - 1) * n/m) + 1):ceiling(j * 
                                                                      n/m)
                  id.train <- setdiff(1:n, id.test)
                  X.train <- Xn[id.train, ]
                  Y.train <- Yn[id.train, ]
                  s.train<-sn[id.train,]
                  n.train<-length(id.train)
                  one.train<-rep(1,n.train)
                  X.test <- Xn[id.test, ]
                  Y.test <- Yn[id.test, ]
                  n.test <- length(id.test)
                  one.test<-rep(1,n.test)
                  
                  ############## SPATIAL ENV
                  fit.sp <-env.fun(X.train, Y.train, u,s.train,cov.form,init.theta,init.error)
                  alpha<-fit.sp$alpha
                  beta<-fit.sp$beta
                  Gamma1<-fit.sp$Gamma1
                  Gamma0<-fit.sp$Gamma0
                  Omega<-fit.sp$Omega
                  Omega0<-fit.sp$Omega0
                  theta<-fit.sp$theta
                  M1<-kronecker(alpha,one.test)+kronecker(D,X.test)%*%as.vector(t(beta))
                  M2<-kronecker(alpha,one.train)+kronecker(D,X.train)%*%as.vector(t(beta))
                  V<-Gamma1 %*% tcrossprod(Omega, Gamma1)+ Gamma0 %*% tcrossprod(Omega0, 
                                                                                              Gamma0)
                  
                  if(cov.form=="matern"){
                        rho<-Matern(Disn , range = theta[1], smoothness =theta[2])
                  }else{
                        rho<-Matern(Disn , range =theta)
                  }
                  
                  R1<-rho[id.test,]
                  R11<-R1[,id.test]
                  R12<-R1[,-id.test]
                  R2<-rho[-id.test,]
                  R22<-R2[,-id.test]
                  R21<-R2[,id.test]   
                  Sig11<-kronecker(V,R11)
                  Sig12<-kronecker(V,R12)
                  Sig22<-kronecker(V,R22)
                              
                  resi.sp <-as.vector(Y.test)-(M1+Sig12%*%solve(Sig22)%*%(as.vector(Y.train)-M2))
                  prederr.sp.CN[i] <- prederr.sp.CN[i] + sum(resi.sp^2)
                  ############## ENV
                  fit <- env(X.train, Y.train, u, asy = F)
                  betahat <- fit$beta
                  muhat <- fit$mu
                  resi.env <- as.matrix(Y.test - matrix(1, n.test, 1) %*% 
                                          t(muhat) - as.matrix(X.test) %*% t(betahat))
                  prederr[i] <- prederr[i] + sum(resi.env^2)
                  
                  ################ OLS
                  Xols<-cbind(1,X.train)
                  Bols<-solve(t(Xols)%*%Xols)%*%t(Xols)%*% Y.train
                  res.ols<-Y.test-cbind(1,X.test) %*%Bols
                  prederr.ols[i] <- prederr.ols[i] + sum(res.ols^2)
            }
      }
      return(list(cv.env=mean(prederr/(r*n)),cv.ols=mean(prederr.ols/(r*n)),cv.sp.env.CN=mean(prederr.sp.CN/(r*n))))
}

   