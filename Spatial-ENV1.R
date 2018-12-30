## These functions gives estimate for spatial parameters and other starting values
## Important packages
library(Renvlp)
library(fields)
######################
Theta.fun<-function(X, Y,u, s,Gamma1,
                    Gamma0,Omega,Omega0, cov.form="matern",init.theta=c(1,0.5)){
  n<-nrow(Y)
  r<-ncol(Y)
  Dis<-as.matrix(dist(s))
  if(u==0){
    invV1<-matrix(0,r,r)
    V0<-Gamma0%*%Omega0 %*% t(Gamma0)
    invV0<-ginv(V0)
  }else if(u==r){
    invV0<-matrix(0,r,r)
    V1<-Gamma1%*%Omega %*% t(Gamma1)
    invV1<-ginv(V1)
  }else{
    V1<-Gamma1%*%Omega %*% t(Gamma1)
    invV1<-ginv(V1)
    V0<-Gamma0%*%Omega0 %*% t(Gamma0)
    invV0<-ginv(V0)
  }
  
  Ybar<-apply(Y,2,mean)
  Xbar<-apply(X,2,mean)
  one<-rep(1,n)
  H<-Y-kronecker(t(Ybar),one)
  G<-X-kronecker(t(Xbar),one)
  
  function(par){
    
    
    
    
    if(sum(par<0)>0){
      fun<-100000000
    }else{
      if(cov.form=="matern"){
        rho<-Matern(Dis , range = par[1], smoothness = par[2])
      }else{
        rho<-Matern(Dis , range = par, smoothness = 0.5)   
      }
      tmp.M<-eigen(rho)
      #print(sum(tmp.M$values<0)>0)
      if(sum(tmp.M$values<0)>0){fun<-100000000}else{
        rho2<-sweep(tmp.M$vectors, MARGIN = 2, -1/sqrt((tmp.M$values)), 
                    "*") %*% t(tmp.M$vectors)
        Hx<-rho2%*%G
        Q<-diag(n)-Hx%*%solve(t(Hx)%*%Hx)%*%t(Hx)
        A1<-Q%*%rho2%*%H
        A2<-rho2%*%H
        fun<-0.5*r*det(rho)+(1/2)*sum(diag(A1%*%invV1%*%t(A1)+A2%*%invV0%*%t(A2)))
      }}
    return(fun)
  }
}
######################
starting.value<-function(X, Y, s, u, cov.form="matern",init.theta=c(1,0.5))
{
  Hat<-env(X, Y, u,asy=F)
  
  Ga<-Hat$Gamma
  Ga0<-Hat$Gamma0
  alpha<-Hat$mu
  beta<-Hat$beta
  eta<-Hat$eta
  Om<-Hat$Omega
  Om0<-Hat$Omega0 
  
  Omega<-Om
  Omega0<-Om0
  Gamma1<-Ga
  Gamma0<-Ga0
  
  Theta<-Theta.fun(X, Y,u, s,Gamma1,
                   Gamma0,Omega,Omega0,cov.form)
  theta<-optim(init.theta, fn=Theta, lower = c(0.1,0.1), upper = c(7,7), method = "L-BFGS-B")$par
  
  return(list(Gamma1 = Gamma1, Gamma0 = Gamma0, alpha = alpha, 
              beta = beta, eta = eta, Omega = Omega, 
              Omega0 = Omega0,theta=theta))
}
