#' Generate MCMC sample from the full conditional distribution of beta
#'
#' @param y The realizations of the true underlying process at monitored
#' stations and sampling points within grid cells.
#' @param x The covariate matrix at monitored stations and sampling points within
#' grid cells.
#' @param prior.mean Prior mean of the coefficient vector beta.
#' @param prior.var.solve Inverse of the prior variance matrix of the coefficient vector beta.
#' @param sigma.solve Inverse of the spatial covariance matrix of true underlying process.
#'
#' @return MCMC samples of beta
#' @import MASS
#' @export
#' @examples
updatebeta<-function(y,x,prior.mean,prior.var.solve,sigma.solve)
  {  B<-solve(t(x)%*%sigma.solve%*%x+prior.var.solve)
     b<-crossprod(x,sigma.solve)%*%y + prior.var.solve%*%prior.mean
     beta1<-mvrnorm(1,mu=B%*%b,B)
     return(beta1)
}




#' Generate MCMC sample from the full conditional distribution of theta, including sill parameter and range parameter.
#'
#' @param diff The residuals of the true underlying process.
#' @param theta The values of the covariance parameters theta from the previous MCMC iteration.
#' @param Distance The Euclidean distance matrix between all the locations including monitored stations and sampling points within grid cells.
#' @param n Number of monitored stations and sampling points within grid cells.
#' @param cov.model One of the three possible choices for the covariance function: “exponential”, “Gaussian” and “matern”.
#'
#' @return MCMC samples of theta
#' @export
#' @import MASS
#' @examples
updatetheta<-function(diff,theta,Distance,n,cov.model)
  {
   if(cov.model=="exponential")
     {  cov.fun<-function(d,range)
         { return(exp(-d/range)) }
     }
    if(cov.model=="Gaussian")
     { cov.fun<-function(d,range)
         { return(exp(-d^2/range)) }
     }
    if(cov.model=="Matern")
     { cov.fun<-function(d,range,smooth)
         { return(matern(d,range,smooth)) }
     }
  if(cov.model %in%c("exponential","Gaussian") )
  {
  sill0<-theta[1]
  range0<-theta[2]
  range1<-rlnorm(1,mean=log(range0),sdlog=0.25)

  sigma<-cov.fun(Distance,range0)
  sigma.chol<-chol(sigma)
  sigma.det<-prod(diag(sigma.chol)^2)
  v0<-backsolve(sigma.chol,diff)

  sigma1<-cov.fun(Distance,range1)
  sigma1.chol<-chol(sigma1)
  sigma1.det<-prod(diag(sigma1.chol)^2)
  v1<-backsolve(sigma1.chol,diff,upper.tri=F)

  sill.new<-1/rgamma(1,shape=2+n/2,scale=1/(2+0.5*(t(v0)%*%v0)) )
  temp<-range0/range1  #range0/range1
  temp<-temp*dgamma(range1,shape=3.5,scale=2)/dgamma(range0,shape=3.5,scale=2)

  temp1<-t(v0)%*%v0-t(v1)%*%v1
  temp1<-temp1/(2*sill.new)
  temp<-temp*exp(temp1)
  if(sigma1.det>0)
   {
     temp<-temp*sqrt(sigma.det)/sqrt(sigma1.det)
     u<-runif(1)
    if (  u<=temp )
        { range.new<-range1
          rate<-1
        }
    if(  u>temp )
       { range.new<-range0 }
    }
  if(sigma1.det<=0)
    { range.new<-range0 }
   return(c(sill.new,range.new))
  }
  if(cov.model=="matern")
    {
  sill0<-theta[1]
  range0<-theta[2]
  smooth0<-theta[3]
  range1<-rlnorm(1,mean=log(range0),sdlog=0.25)
  smooth1<-rlnorm(1,mean=log(smooth0),sdlog=0.25)

  sigma<-cov.fun(Distance,range0,smooth0)
  sigma.chol<-chol(sigma)
  sigma.det<-prod(diag(sigma.chol)^2)
  v0<-backsolve(sigma.chol,diff)

  sigma1<-cov.fun(Distance,range1,smooth0)
  sigma1.chol<-chol(sigma1)
  sigma1.det<-prod(diag(sigma1.chol)^2)
  v1<-backsolve(sigma1.chol,diff)

  sill.new<-1/rgamma(1,shape=2+n/2,scale=1/(2+0.5*(t(v0)%*%v0)) )
  temp<-range0/range1  #range0/range1
  temp<-temp*dgamma(range1,shape=3.5,scale=2)/dgamma(range0,shape=3.5,scale=2)
  temp1<-t(v0)%*%v0-t(v1)%*%v1
  temp1<-temp1/(2*sill.new)
  temp<-temp*exp(temp1)
  if(sigma1.det>0)
  {
  temp<-temp*sqrt(sigma.det)/sqrt(sigma1.det)
  u<-runif(1)
  if (  u<=temp )
       { range.new<-range1
         rate<-1
       }
  if(  u>temp )
       { range.new<-range0 }
   }
  if(sigma1.det<=0) {range.new<-range0}
  sigma<-cov.fun(Distance,range.new,smooth0)
  sigma.chol<-chol(sigma)
  sigma.det<-prod(diag(sigma.chol)^2)
  v0<-backsolve(sigma.chol,diff)
  sigma1<-cov.fun(Distance,range.new,smooth1)
  sigma1.chol<-chol(sigma1)
  sigma1.det<-prod(diag(sigma1.chol)^2)
  v1<-backsolve(sigma1.chol,diff)
  temp<-smooth0/smooth1  #range0/range1
  temp<-temp*dgamma(smooth1,shape=0.25,scale=2)/dgamma(range0,shape=0.25,scale=2)
  temp1<-t(v0)%*%v0-t(v1)%*%v1
  temp1<-temp1/(2*sill.new)
  temp<-temp*exp(temp1)
  if(sigma1.det>0)
  {
  temp<-temp*sqrt(sigma.det)/sqrt(sigma1.det)
  u<-runif(1)
  if (  u<=temp )
       { smooth.new<-smooth1
         rate<-1
       }
  if(  u>temp )
       { smooth.new<-smooth0 }
   }
  if(sigma1.det<=0) { smooth.new<-smooth0 }
    return(c(sill.new,range.new,smooth.new))
    }
}

#' Generate MCMC sample from the full conditional distribution of sigma_delta^2 and sigma_e^2
#'
#' @param zhat The measurements vector.
#' @param Zbtilde The modeling outputs vector.
#' @param nm Number of sampling points in each grid cell of the modeling output;
#' @param n  Number of monitored stations and sampling points within grid cells;
#' @param m2 Number of grid cells.
#' @param y The realizations of the true underlying process at monitored
#' stations and sampling points within grid cells.
#' @param a Additive bias parameter a
#' @param b Multiplicative bias parameter b
#' @param A2 Matrix A2,
#' @export
#' @return MCMC sample of sigma
#' @import MASS
#' @examples
updatesigma<-function(zhat,Zbtilde,nm,n,m2,y,a,b,A2)
  { sca<-1/(1+0.5*t(zhat-y[(nm*m2+1):n])%*%(zhat-y[(nm*m2+1):n]))
    sigmaesquare<-1/rgamma(1,3+(n-nm*m2)/2,scale=sca)
    Zb<-A2%*%y
    uv<-a+b*Zb
    sca<-1/(1+0.5*t(Zbtilde-uv)%*%(Zbtilde-uv) )
    sigmadsquare<-1/rgamma(1,3+m2/2,scale=sca)
    return( c(sigmadsquare,sigmaesquare) )
}

#' Generate MCMC sample from the full conditional distribution of a and b
#'
#' @param Zbtilde the modeling outputs vector.
#' @param sigmad the output error variance parameter sigma_delta^2
#' @param y the realizations of the true underlying process at monitored
#' stations and sampling points within grid cells.
#' @param ab0 mean of normal prior of vector including a and b
#' @param fb.solve inverse of variance matrix of normal prior of vector including a and b
#' @param nm number of sampling points in each grid cell of the modeling output;
#' @param m2 number of grid cells.
#' @param A2 matrix A2,

#' @return MCMC samples of a,b
#' @import MASS
#' @export
#' @examples
updateab<-function(Zbtilde,sigmad,y,ab0,fb.solve,nm,m2,A2)
  { ### ab0: prior mean of a and b; fb: prior variance of a and b
    Zb<-A2%*%y
    Xb<-cbind(1,Zb)
    B<-solve( crossprod(Xb,Xb)/sigmad+fb.solve )
    b<-crossprod(Xb,Zbtilde)/sigmad+crossprod(fb.solve,ab0)
    abtemp<-mvrnorm(1,mu=B%*%b,B)
    a<-abtemp[1]
    s<-abtemp[2]
    return( c(a,s) )
  }


#' Generate MCMC sample from the full conditional distribution of true underlying process Z
#' @param X The covariate matrix at monitored stations and sampling points within grid cells.
#' @param zhat The measurements vector.
#' @param Zbtilde The modeling outputs vector.
#' @param A1 Matrix A1
#' @param A2 Matrix A2
#' @param beta Regression coefficient
#' @param a Additive bias parameter a
#' @param b Nultiplicative bias parameter b
#' @param nm Number of sampling points in each grid cell of the modeling output;
#' @param m2 Number of grid cells.
#' @param sigmad The output error variance parameter sigma_delta^2
#' @param sigmae Measurement error variance parameter sigma_e^2
#' @param sigma3.solve Inverse of the spatial covariance matrix of true underlying process
#'
#' @return MCMC sample of true underlying process Z
#' @export
#' @import MASS
#' @examples
updatezs<-function(X,zhat,Zbtilde,A1,A2,beta,a,b,nm,m2,sigmad,sigmae,sigma3.solve)
  { ## zhelicuole
    Zbtilde1<-(Zbtilde-a) ## output after removing the additive bias
    Mu<-X%*%beta
    temp<-crossprod(A1,A1)/sigmae+b*b*crossprod(A2,A2)/sigmad+sigma3.solve
    temp.chol<-chol(temp)
    sigmastar<-chol2inv(temp.chol)
    ustar<-sigmastar%*%(crossprod(A1,zhat)/sigmae+b*crossprod(A2,Zbtilde1)/sigmad+sigma3.solve%*%Mu)
    Zs<-mvrnorm(1, ustar, sigmastar)
    return(Zs)
  }
