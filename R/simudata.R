#' Title
#'
#' @param m1
#' @param m2
#' @param n.pred Number of unmonitored locations for prediction
#' @param ran
#' @param nm
#' @param seed
#'
#' @return
#' @export
#' @import MASS
#' @import geoR
#'
#' @examples
simudata_data_position <- function(m1, m2, n.pred , ran , nm , seed = 1234){
  set.seed(seed)
  sloc <- cbind( runif(m1,min=-ran,max=ran),runif(m1,min=-ran,max=ran) )
  sam.sloc1 <- cbind( runif(n.pred,min=-ran,max=ran),runif(n.pred,min=-ran,max=ran))
  sam <- cbind( runif(m2,min=-ran,max=ran),runif(m2,min=-ran,max=ran) )
  sam.sloc <- rbind(sam,sloc)
  return(list(sloc = sloc,
              sam = sam,
              sam.sloc = sam.sloc,
              sam.sloc1 = sam.sloc1))
}

#' Title
#'
#' @param m1
#' @param m2
#' @param n.pred
#' @param ran
#' @param nm
#' @param nhours
#' @param seed
#'
#' @return
#' @export
#' @import MASS
#' @import geoR
#'
#' @examples
simulate_data <- function(m1, m2, n.pred , ran , nm , nhours, seed = 1234){
  sloc <- simudata_data_position(m1, m2, n.pred , ran , nm , seed)$sloc
  sam <- simudata_data_position(m1, m2, n.pred , ran , nm , seed)$sam
  sam.sloc <- simudata_data_position(m1, m2, n.pred , ran , nm , seed)$sam.sloc
  sam.sloc1 <- simudata_data_position(m1, m2, n.pred , ran , nm , seed)$sam.sloc1
  Theta<-c(1.5,5,0.5)
  Beta<-c(2.5,2.9,3.2)
  true<-matrix(rep(0,n.pred*nhours),ncol=nhours)
  airs<-matrix(rep(0,m1*nhours),ncol=nhours) # measured stations mat
  simu<-matrix(rep(0,m2*nhours),ncol=nhours) # grid cells mat
  true<-as.matrix(true)
  airs<-as.matrix(airs)
  simu<-as.matrix(simu)
  n<-length(sam.sloc[,1]) ########## number of stations + points within cells
  n.total<-n+n.pred
  n.pre<-length(sam.sloc1[,1]) ############ mark1
  X<-cbind(rep(1,n),sam.sloc[,1],sam.sloc[,2])
  X1<-cbind(rep(1,n.pred),sam.sloc1[,1],sam.sloc1[,2])
  sam.sloc2<-rbind(sam.sloc,sam.sloc1)
  X2<-rbind(X,X1)   # covariates for the available stations+sampling points within grid cells+stations to be predicted
  distance1<-as.matrix(dist(sam.sloc2))
  set.seed(seed)
  #############################################################################
  for(kk in 1:nhours)
  { #print(kk)
  #################  simulate the data  ##############################
    SIGMA<-Theta[1]*matern(distance1,Theta[2],Theta[3])
    ZS.underly<-mvrnorm(1,mu=X2%*%Beta,SIGMA)
    airs[,kk]<-ZS.underly[(nm*m2+1):n]+rnorm((n-nm*m2),0,0.5)
    temp1<-5+2.5*ZS.underly[1:(nm*m2)]+rnorm(nm*m2,0,0.5)
    temp<-matrix(temp1,ncol=nm,byrow=T)
    simu[,kk]<-apply(temp,1,mean)
    true[,kk]<- ZS.underly[(n+1):n.total]+rnorm(n.pred,0,0.5)
  }
  return(list(simu = simu,
         true = true,
         airs = airs,
         sam.sloc = sam.sloc,
         sam.sloc1 = sam.sloc1))
}
