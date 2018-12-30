# This function Provide the 
#Asymptotic variance for sqrt(n)*vec(B)

Asy_var<-function(X, Y, u,s,cov.form,init.theta,init.error)
{
  r<-dim(Y)[2]
  B<-GHn(r)
  Cr<-B$Hn
  Er<-B$Gn
  XX<-kronecker(diag(r),X)
  
  result<-env.fun(X, Y, u,s,cov.form,init.theta,init.error) 
  Gamma1<-result$Gamma1
  Gamma0<-result$Gamma0
  alpha<-result$alpha
  Omega<-result$Omega
  Omega0<-result$Omega0
  theta<-result$theta
  beta<-result$beta
  eta<-result$eta
  V<-result$Sigma
  
  Dis<-as.matrix(dist(s))
  if(cov.form=="matern"){
    rho<-Matern(Dis , range = theta[1], smoothness =theta[2])
  }else{
    rho<-Matern(Dis , range =theta)   
  }
  
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
  
  asy_var<-kronecker(A1,A2)+A3%*%chol2inv(t(Si)%*%J%*%Si)%*%t(A3)
  matrix(sqrt(diag(AV)), nrow = r)
  return(asy_var)
  
}

############
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

