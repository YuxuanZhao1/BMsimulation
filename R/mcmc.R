#' Run BM model and return a list of outputs
#'
#' @param m1 Number of monitoring locations
#' @param m2 Number of grid cells
#' @param nm Number of points with in each of grid cells
#' @param sam.sloc The coordinates of sampling points and the monitoring stations
#' @param sam.sloc1 The coordinates of the unmonitored stations where the measurements for prediction
#' @param ton Number of MCMC iterations in the Gibbs sampling
#' @param burnin The “burn-in” period of the Gibbs sampling.
#' @param zhat The measurements vector.
#' @param Zbtilde The modeling outputs vector
#' @param degree Degree of the polynomial function, 0 the mean is assumed constant across space,
#' 1 the mean is assumed to be a first order polynomial on the coordinates,
#' 2 the mean is assumed to be a second order polynomial on the coordinates.
#' @param cov.model one of the three possible choices for the covariance function: “exponential”, “Gaussian” and “matern”.
#'
#' @return A list of estimated model parameters
#' @export
#' @import MASS
#' @import geoR
#'
#' @examples
melding<-function(m1,m2,nm,sam.sloc,sam.sloc1,ton,burnin,zhat,Zbtilde,degree,cov.model){
   q1<-function(vec) {return(quantile(vec,0.05))} ## calulate the 5% and 95 quantiles
   q2<-function(vec) {return(quantile(vec,0.95))} ##
    print(degree)
    n<-dim(sam.sloc)[1]
    n.pred<-dim(sam.sloc1)[1]
    ntotal<-m1+m2*nm+n.pred
    distance<-as.matrix( dist(sam.sloc) )
    distance.vec<-as.vector(distance)
    sam.sloc2<-rbind(sam.sloc1,sam.sloc)
    distance1<-as.matrix( dist(sam.sloc2) )
    distance1.vec<-as.vector(distance1)
    if( (degree!=1)&(degree!=2) ) { print("degree has to be 1 or 2!"); return(NULL); }
    if(degree==1){
         X<-cbind(rep(1,n),sam.sloc[,1],sam.sloc[,2])
         X1<-cbind(rep(1,n.pred),sam.sloc1[,1],sam.sloc1[,2])
      }
    if(degree==2){
        X<-cbind(rep(1,n),sam.sloc[,1],sam.sloc[,2],sam.sloc[,1]^2,sam.sloc[,2]^2,sam.sloc[,1]*sam.sloc[,2])
        X1<-cbind(rep(1,n.pred),sam.sloc1[,1],sam.sloc1[,2],sam.sloc1[,1]^2,sam.sloc1[,2]^2,sam.sloc1[,1]*sam.sloc1[,2])
      }
    nbeta<-dim(X)[2]
    cov.model.list<-c("exponential","Gaussian","Matern")
    if((cov.model%in%cov.model.list)=="FALSE" ){
    print("No such covariance function available. Please choose from
         exponential, Gaussian or Matern ")
        }
    if(cov.model=="exponential"){
      ntheta<-2
       cov.fun<-function(d,range)
         { return(exp(-d/range)) }
     }
    if(cov.model=="Gaussian")
     { ntheta<-2
       cov.fun<-function(d,range)
         { return(exp(-d^2/range)) }
     }
    if(cov.model=="Matern")
     { ntheta<-3

       cov.fun<-function(d,range,smooth)
         { return(matern(d,range,smooth)) }
     }

    nab<-2
    nerror<-2
    predict<-rep(0,n.pred)
    predictsd<-predict
    predicttemp<-matrix(rep(0,ton*n.pred),ncol=ton)    ## E(Z)
    predicttemp1<-matrix(rep(0,ton*n.pred),ncol=ton)   ## E(Z^2)
    theta.est<-rep(0,ntheta)
    theta.est.sd<-rep(0,ntheta)
    ab.est<-rep(0,nab)
    ab.est.sd<-rep(0,nab)
    error.est<-rep(0,nerror)
    error.est.sd<-rep(0,nerror)
    a<-rep(1,ton)
    s<-rep(1,ton)
    sigmad<-rep(1,ton)  ## sigma delta square
    sigmae<-rep(1,ton)  ## sigma e square
    thetaouter<-rbind(rep(1,ton),rep(1,ton))
    beta<- matrix(rep(0,ton*nbeta),ncol=ton)   ##  beta matrix
    rate<-rep(0,ton)      ## m-h acceptance rate
    Zsmat<-matrix(rep(1,ton*n),ncol=ton,byrow=T)
    ###### prior of (a,b) is N(ab0,f1)#########
     ab0<-c(0,1)
     fb<-diag(rep(100,2))
     fb.solve<-solve(fb)
    ### the prior of beta is N(betabar,f) #######
     betabar<-rep(1,nbeta)
     F.solve<-diag(rep(0.1,nbeta))
    ############ set up the matrices A1 and A2 ##################
    temp1<-matrix(rep(0,nm*m2*m1),ncol=nm*m2)
    temp2<-diag( rep(1,m1) )
    A1<-cbind(temp1,temp2)
    A2<-matrix(rep(0,m2*(m1+m2*nm) ),nrow=m2)
    for(i in 1:m2){ A2[i,(1+(i-1)*nm):(i*nm)]<-1/nm  }

    for(k in 2:ton)
        {  #print(k)
           #######################  update theta ##############################
           diff<-Zsmat[,k-1]-X%*%beta[,k-1]
           thetaouter[,k]<-updatetheta(diff,thetaouter[,k-1],distance,n,cov.model)
           ############### update underlying true process Z #####################
           if(cov.model %in%c("exponential","Gaussian") )
              { sigma3<-thetaouter[1,k]*cov.fun(distance,thetaouter[2,k]) }
           if(cov.model=="matern" )
              { sigma3<-thetaouter[1,k]*cov.fun(distance,thetaouter[2,k],thetaouter[3,k]) }
           sigma3.chol<-chol(sigma3)
           Sigma3.solve<-chol2inv(sigma3.chol)
           Zsmat[,k]<-updatezs(X,zhat,Zbtilde,A1,A2,beta[,k-1],a[k-1],s[k-1],nm,m2,sigmad[k-1],sigmae[k-1],Sigma3.solve)
           #######################  update beta  ##############################
           beta[,k]<-updatebeta(Zsmat[,k],X,betabar,F.solve,Sigma3.solve)
           ######## update sigma e square and sigma delta square ######
           sigmade<-updatesigma(zhat,Zbtilde,nm,n,m2,Zsmat[,k],a[k-1],s[k-1],A2)
           sigmad[k]<-sigmade[1]
           sigmae[k]<-sigmade[2]
           ######## update biases a and b ##############################
           ab<-updateab(Zbtilde,sigmad[k],Zsmat[,k],ab0,fb.solve,nm,m2,A2)
           a[k]<-ab[1]
           s[k]<-ab[2]
           ##########  prediction in Gibbs Sampling loop #############
           if(cov.model %in%c("exponential","Gaussian") )
             { Sigma<-thetaouter[1,k]*cov.fun(distance1,thetaouter[2,k]) }
           if(cov.model=="matern" )
             { Sigma<-thetaouter[1,k]*cov.fun(distance1,thetaouter[2,k],thetaouter[3,k]) }
           Sigma12<-Sigma[1:n.pred,(n.pred+1):ntotal]
           Sigma11<-Sigma[1:n.pred,1:n.pred]
           predict.mean<-X1%*%beta[,k]+Sigma12%*%Sigma3.solve%*%(Zsmat[,k]-X%*%beta[,k])
           predict.var<-Sigma11-Sigma12%*%Sigma3.solve%*%t(Sigma12)+diag(rep(sigmae[k],n.pred))
           predicttemp[,k]<-mvrnorm(1,mu=predict.mean,predict.var)
        }
   index<-seq(from=burnin,to=ton,by=5)
   theta.est<- apply(thetaouter[,index],1,mean)
   theta.est.sd<- apply(thetaouter[,index],1,sd)
   predict<-apply(predicttemp[,index],1,mean)
   predict.q1<-apply(predicttemp[,index],1,q1)
   predict.q2<-apply(predicttemp[,index],1,q2)
   ab.est<-c(mean(a[index]),mean(s[index]))
   ab.est.sd<-c(sd(a[index]),sd(s[index]))
   beta.est<- apply(beta[,index],1,mean)
   beta.est.sd<-apply(beta[,index],1,sd)
   sigmae.est<-mean(sigmae[index])
   sigmae.est.sd<-sd(sigmae[index])
   sigmad.est<-mean(sigmad[index])
   sigmad.est.sd<-sd(sigmad[index])
   result<-list(beta.est,beta.est.sd,theta.est,theta.est.sd,predict,predict.q1,predict.q2,ab.est,ab.est.sd,sigmae.est,
                sigmae.est.sd,sigmad.est,sigmad.est.sd, a[1:ton], s[1:ton],
                beta[,1:ton], thetaouter[,1:ton])
   names(result)<-c("beta.est","beta.est.sd","theta.est","theta.est.sd","prediction","pred.q1",
                    "pred.q2","ab.est","ab.est.sd",
                    "sigmae.est","sigmae.est.sd","sigmad.est","sigmad.est.sd",
                    "a.trace","b.trace","beta.trace", "theta.trace")
   return(result)
  }
