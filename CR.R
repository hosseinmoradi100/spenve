## This Function Generates the leave one out Cross validation for 
## Multivariate Regression, Envelope, and Spatial Envelope
## 
## The inputs for this function are
## X: denotes the matrix of covariates
## Y: denotes the matrix of response
## u: denotes the structural dimension of the envelope
## s: denotes the location
## Cov.form: denotes the form of spatial covariance function
## init.theta: denotes the initial values for spatial covariance dunction
## init.error: denotes the initial error value

CR.fun<-function(X, Y, u,s,cov.form="matern",init.theta=c(1,0.5),init.error=0.0005)
{
      n<-nrow(Y)
      r<-ncol(Y)
      one<-rep(1,n-1)
      D<-diag(r)
      Dis<-as.matrix(dist(s))
      theta=init.theta
      if(cov.form=="matern"){
            rho<-Matern(Dis , range = theta[1], smoothness =theta[2])
      }else{
            rho<-Matern(Dis , range =theta)
      }
      
      Yihat<-matrix(0,nrow=n,ncol=r)
      Yhat<-matrix(0,nrow=n,ncol=r)
      Yhat.env<-matrix(0,nrow=n,ncol=r)
      Yhat.ols<-matrix(0,nrow=n,ncol=r)
      for(i in 1:n)
      {
            print(i)
            Yi<-Y[-i,]
            Xi<-X[-i,]
            si<-s[-i,]
            result<-env.fun(Xi, Yi, u,si,cov.form,init.theta,init.error) 
            Gamma1<-result$Gamma1
            Gamma0<-result$Gamma0
            alpha<-result$alpha
            Omega<-result$Omega
            Omega0<-result$Omega0
            theta<-result$theta
            beta<-result$beta
             if(cov.form=="matern"){
                           rho<-Matern(Dis , range = theta[1], smoothness =theta[2])
                    }else{
                          rho<-Matern(Dis , range =theta)   
                    }
            M1<-alpha+beta%*%X[i,]
            M2<-kronecker(alpha,one)+kronecker(D,X[-i,])%*%as.vector(t(beta))
            
            
            V<-Gamma1 %*% tcrossprod(Omega, Gamma1)+ Gamma0 %*% tcrossprod(Omega0, 
                                                          Gamma0)
            
            
            R1<-rho[i,]
            R11<-R1[i]
            R12<-R1[-i]
            R2<-rho[-i,]
            R22<-R2[,-i]
            R21<-R2[,i]
            
            
            Sig11<-V*R11
            Sig12<-kronecker(V,R12)
            Sig22<-kronecker(V,R22)
            
            Yhat[i,]<-M1+t(Sig12)%*%solve(Sig22)%*%(as.vector(Yi)-M2)
            
            fit <- env(Xi, Yi, u, asy = F)
            betahat <- fit$beta
            muhat <- fit$mu
            Yhat.env[i,] <-muhat + betahat%*%X[i,]
            
            X0<-cbind(one,Xi)
            Bols<-solve(t(X0)%*%X0)%*%t(X0)%*%Yi
            
            Yhat.ols[i,]<-c(1,X[i,])%*%Bols
            
      }
      res.spenv<-Y-Yhat
      res.env<-Y-Yhat.env
      res.ols<-Y-Yhat.ols
      cr.spenv<-mean(res.spenv^2)
      cr.env<-mean(res.env^2)
      cr.ols<-mean(res.ols^2)
      
      return(list(CR.sp=cr.spenv,cr.env=cr.env,cr.ols=cr.ols))
}

