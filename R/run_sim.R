#' Title
#'
#' @param m1
#' @param m2
#' @param n.pred
#' @param ran
#' @param nm
#' @param nhours
#' @param ton
#' @param burnin
#' @param degree
#' @param nbeta
#' @param ntheta
#' @param nab
#' @param nerror
#' @param cov.model
#' @param seed
#' @param num_for_cores
#'
#' @return
#' @import MASS
#' @import geoR
#' @import doParallel
#' @import foreach
#' @import parallel
#'
#' @export
#'
#'
#'
#' @examples
run_func <- function(m1=20, m2=2, n.pred = 100, ran = 5, nm = 1, nhours = 50,
                     ton = 1000, burnin = ton/10, degree = 1, nbeta= 3,
                     ntheta= 2, nab=2, nerror=2, cov.model = "exponential", seed,
                     num_for_cores = 6){

  cl <- makeCluster(num_for_cores)
  # clusterExport(cl = cl, list("melding", "q1", "q2", "simulate_data","simudata_data_position",
  #                             "updateab", "updatebeta", "updatesigma", "updatetheta",
  #                             "updatezs"),
  #               envir=environment())
  # clusterEvalQ(cl, library(MASS))


  registerDoParallel(cl)

  predictmat<-matrix(rep(0,nhours*n.pred),ncol=nhours)
  predictq1mat<-matrix(rep(0,nhours*n.pred),ncol=nhours)
  predictq2mat<-matrix(rep(0,nhours*n.pred),ncol=nhours)
  betamat<-matrix(rep(0,nbeta*nhours),ncol=nhours)
  betasdmat<-betamat
  thetamat<-matrix(rep(0,ntheta*nhours),ncol=nhours)
  ktheta<-thetamat
  thetasdmat<-thetamat
  abmat<-matrix(rep(0,nab*nhours),ncol=nhours)
  errormat<-matrix(rep(0,nerror*nhours),ncol=nhours)
  absdmat<-abmat
  errorsdmat<-matrix(rep(0,nerror*nhours),ncol=nhours)
  error<-rep(0,nhours)
  ind<-1:n.pred
  out<-rep(0,nhours)
  simulated_data <- simulate_data(m1, m2, n.pred, ran, nm, nhours, seed)
  allresult <- foreach (i = 1:nhours, .errorhandling = "pass") %dopar%
  {
    result<-melding(m1,m2,nm,simulated_data$sam.sloc,simulated_data$sam.sloc1,
                    ton,burnin,simulated_data$airs[,i],simulated_data$simu[,i],
                    degree,cov.model)

    result
  }

  stopCluster(cl)
  return(list(allresult,
              simulated_data))

}
