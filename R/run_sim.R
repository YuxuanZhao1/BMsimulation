#' Run the simulation studies
#'
#' @param m1 Number of monitoring locations
#' @param m2 Number of grid cells
#' @param n.pred Number of unmonitored locations for prediction
#' @param ran Range for the monitoring/unmonitored locations and sampling points within grid cells, Lower bound is -ran and upper bound is ran
#' @param nm Number of points with in each of grid cells
#' @param nhours Number of independent data sets during simulation
#' @param ton Number of MCMC iterations in the Gibbs sampling
#' @param burnin The “burn-in” period of the Gibbs sampling.
#' @param degree Degree of the polynomial function, 0 the mean is assumed constant across space,
#' 1 the mean is assumed to be a first order polynomial on the coordinates,
#' 2 the mean is assumed to be a second order polynomial on the coordinates.
#' @param nbeta Length of regression coefficient vector
#' @param ntheta Length of covariance parameter vector
#' @param nab Length of calibration parameter vector
#' @param nerror Length of error vector
#' @param cov.model one of the three possible choices for the covariance function: “exponential”, “Gaussian” and “matern”.
#' @param seed Number of random seed
#' @param num_for_cores Number of cores for parallel computing
#'
#' @return A lists: the first element will consist of output from melding function and the second list will consist of information about location.

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
#' @examples run_func(m1=20, m2=2, n.pred = 100, ran = 5, nm = 1, nhours = 50,ton = 1000, burnin = ton/10, degree = 1, nbeta= 3,ntheta= 2, nab=2, nerror=2, cov.model = "exponential", seed,num_for_cores = 6)
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
