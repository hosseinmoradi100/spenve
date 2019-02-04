## This is an example for running spatial envelope 


### Required Packages
library(rstiefel)
library(fields)
library(MASS)
library(Renvlp)
library(matrixcalc)
library(pracma)
######################

GHn<-function(n)
{
  k<-n*(n+1)/2
  Gn<-matrix(0,nrow=n^2,ncol=k)
  z1<-rep(0,k)
  for(j in 1:n)
    for(i in j:n)
    {
      Ik<-z1
      Ik[(j-1)*(n-j/2)+i]<-1
      Gn[(j-1)*n+i,]<-Ik
      Gn[(i-1)*n+j,]<-Ik
    }
  
  Hn<-solve(t(Gn)%*%Gn)%*%t(Gn)
  return(list(Gn=Gn,Hn=Hn))
}


# Set the working directory to the folder that you put the R files
setwd()
source("Spatial-ENV1.R") # Functions to generat starting values
source("Spatial-envelope1.R") # Functions that estimates the spatial envelope model
source("CR.R") # Leave one out Cross validation function
source("k_fold_cr.R") # M-fold Cross-validation (More general version of the CR.R)

################################################

n <- 100 # sample size

arange <- 0.5 # spatial range parameter
a1 <- 0.2 # portion of immaterial part
a2 <- 0.2 # spatial nugget 





r<- 5 # dimension of response
p <- 6 # dimension of covariate
u <- 2 # dimension of the envelope

U <- randortho(r)  # Matrix of (\Gamma_1,\Gamma_0)

g0 <- U[,(u+1):r] # \Gamma_0
g1 <- U[,1:u] # \Gamma_1
eta<-matrix(rnorm(2*p),p,2) #eta
beta <- g1%*%t(eta) # regression coefficients
D1 <- abs(outer(1:u,1:u,"-")) 
D2 <- abs(outer(1:(r-u),1:(r-u),"-"))
omega1 <- (-0.9)^D1 # Omega1
omega0 <- (-0.5)^D2 # Omega0
sig1 <- g1%*%omega1%*%t(g1) # Variance of the material part
sig0 <- a1*g0%*%omega0%*%t(g0) # Variance of the immaterial part


Sig <- sig0+sig1 # Multivariate variance 


###### saving the results 

### time

tols <- 0
tenv <- 0
tspenv <- 0

###MSE

mseols <- 0
mseenv <- 0
msespenv <- 0


### CV (Leave one out)
cvols <- 0
cvenv <- 0
cvspenv <- 0

### CV (M-fold)

m=10 ## number of fold for cross-validation

cv.env.m<-cv.sp.env.m<-cv.sp.env.CN.m<- 0


### beta error

betaerrols <- 0
betaerrenv <- 0
betaerrspenv <- 0


# Distance between two spaces (r^2) [Ye and Weiss, 2003]

dbtenv <- 0
dbtspenv <- 0

## Assymptotic Variance

asympvarenv <- list()
asympvarspenv <- list()


### Data generation


one <- rep(1,n)  # intercept

for(i in 1:100){
  s<-loc<- matrix( runif(n*2), ncol=2) # Location

  d <- as.matrix(rdist(loc,loc)) # Distance matrix

  rhoran <- a2*Matern(d, range = arange, smoothness =0.5 ) # Matern covariance function

  SigYran <- kronecker(Sig,rhoran) # covariance of the response variable

  
  e <- matrix(mvrnorm(1,c(rep(0,n*r)),SigYran),n) # Error

  X<- matrix(rnorm(p*n),n,p)  # Covariates

  Y <- X%*%t(beta) + e

  
  CR.result<-CR.fun(X, Y, u,s,cov.form="matern",init.theta=c(5,0.5),init.error=0.0005)
  CR.result.m<-cv.spenv(X, Y, u,s, m, nperm=10,cov.form="matern",init.theta=c(1,0.5),init.error=0.0005)
  
  cv.env.m[i]<-CR.result.m$cv.env
  cv.sp.env.m[i]<-CR.result.m$cv.sp.env
  cv.sp.env.CN.m[i]<-CR.result.m$cv.sp.env.CN
  
  
  # OLS
  
  tols[i] <- system.time(betahatols <- solve(t(X)%*%X)%*%t(X)%*%Y)[3]
  Yhatols <- X%*%betahatols
  Ytesthatols <- Xtest%*%betahatols
  
  
  mseols[i] <- mean((Y - Yhatols)^2)
  betaerrols[i] <- mean((beta-t(betahatols))^2)
  cvols[i]<-CR.result$cr.ols
  
  
  
  #envelope 
  
  tenv[i] <- system.time(envest<-env(X, Y, u, asy = TRUE, init = NULL))[3]
  betaenv <- envest$beta
  alphaenv <- envest$mu
  
  Yhatenv <- matrix(kronecker(alphaenv,one),n) + X%*%t(betaenv)

  mseenv[i] <- mean((Y-Yhatenv)^2)
  betaerrenv[i] <- mean((beta - betaenv)^2)
  dbtenv[i] <- 0.5*sum(diag(t(envest$Gamma)%*%g1%*%t(g1)%*%envest$Gamma))
  asympvarenv[[i]] <- envest$asySE
  cvenv[i]<-CR.result$cr.env

  
  #spatial envelope
  tspenv[i] <- system.time(envsp <- env.fun(X, Y, u, s=loc, cov.form="matern", init.theta=c(1,0.5), init.error=0.0005))[3]
  betaspenv <- envsp$beta
  alphaspenv <- envsp$alpha
  V <- envsp$Sigma
  theta <- envsp$theta
  Gamma1<-envsp$Gamma1
  
  
  Yhatspenv <- matrix(kronecker(alphaspenv,one),n)+X%*%t(betaspenv)

  msespenv[i] <- mean((Y-Yhatspenv)^2)

  betaerrspenv[i] <- mean((beta-betaspenv)^2)
  dbtspenv[i] <- 0.5*sum(diag(t(envsp$Gamma1)%*%g1%*%t(g1)%*%envsp$Gamma1))
  
  cvspenv[i]<-CR.result$CR.sp

  
  r<-dim(Y)[2]
  B<-GHn(r)
  Cr<-B$Hn
  Er<-B$Gn
  XX<-kronecker(diag(r),X)
  
  Gamma1<-envsp$Gamma1
  Gamma0<-envsp$Gamma0
  alpha<-envsp$alpha
  Omega<-envsp$Omega
  Omega0<-envsp$Omega0
  theta<-envsp$theta
  beta<-envsp$beta
  eta<-envsp$eta
  V<-envsp$Sigma
  
  Dis<-as.matrix(dist(loc))
  rho<-Matern(Dis , range = theta[1], smoothness =theta[2])
  A1<-solve(t(X)%*%solve(rho)%*%X)
  A2<-Gamma1%*%Omega%*%t(Gamma1)
  A3<-kronecker(t(eta),Gamma0)
  Si21<-kronecker(t(eta),Gamma0)
  Si22<-2*Cr%*%(kronecker(Gamma1%*%Omega,Gamma0)-kronecker(Gamma1,Gamma0%*%Omega0))
  J11<-t(XX)%*%kronecker(solve(V),solve(rho))%*%XX
  J22<-0.5*t(Er)%*%kronecker(solve(V),solve(V))%*%Er
  a1<-dim(J11)
  a2<-dim(J22)
  b1<-a1[1]+a2[1]
  b2<-a1[2]+a2[2]
  J<-matrix(0,nrow=b1,ncol=b2)
  J[1:a1[1],1:a1[2]]<-J11
  J[(a1[1]+1):b1,(a1[2]+1):b2]<-J22
  Si<-rbind(Si21,Si22)
  
  
  asympvarspenv[[i]] <- matrix(diag(kronecker(A1,A2)+A3%*%chol2inv(chol(t(Si)%*%J%*%Si))%*%t(A3)),r)
  
  
  
  
  print(i)
}


### For OLS

mean(tols)
mean(mseols)
mean(mspeols)
mean(cvols)
mean(betaerrols)
mean(cv.ols.m)


sd(tols)
sd(mseols)
sd(mspeols)
sd(cvols)
sd(betaerrols)
sd(cv.ols.m)

### For ENV

mean(tenv)
mean(mseenv)
mean(mspeenv)
mean(dbtenv)
mean(cvenv)
mean(betaerrenv)
mean(cv.env.m)

sd(tenv)
sd(mseenv)
sd(mspeenv)
sd(dbtenv)
sd(cvenv)
sd(betaerrenv)
sd(cv.env.m)

### For Spatial ENV

mean(tspenv)
mean(msespenv)
mean(dbtspenv)
mean(cvspenv)
mean(betaerrspenv)
mean(cv.sp.env.CN.m)


sd(tspenv)
sd(msespenv)
sd(dbtspenv)
sd(cvspenv)
sd(betaerrspenv)
sd(cv.sp.env.CN.m)

comptime <- cbind(tols,tenv, tspenv)
compmse <- cbind(mseols,mseenv,msespenv)
compdbt <- cbind(dbtenv,dbtspenv)
compcv <- cbind(cvols, cvenv, cvspenv)
compbeta <- cbind(betaerrols,betaerrenv,betaerrspenv)
compcvmfold <- cbind(cv.ols.m,cv.env.m,cv.sp.env.CN.m)

