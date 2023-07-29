#' Make the traceplot for all of the estimated model parameters
#'
#' @param result The first element from the output of the melding function

#' @param num_data_set The number of independent during simulation, take the number from 1 to nhours
#'
#' @return Traceplot of all model parameters
#' @export
#' @import ggplot2
#' @import latex2exp
#' @import gridExtra
#' @examples traceplot_plt(compute_param_error(), 1)
traceplot_plt = function(result, num_data_set){
  len = length(result[[1]][[num_data_set]]$a.trace)
  data_trace = data.frame(MCMC.interactions = seq(1, len, 1),
                          a = result[[1]][[num_data_set]]$a.trace,
                          b = result[[1]][[num_data_set]]$b.trace,
                          beta_0 =  result[[1]][[num_data_set]]$beta.trace[1,],
                          beta_1 =  result[[1]][[num_data_set]]$beta.trace[2,],
                          beta_2 =  result[[1]][[num_data_set]]$beta.trace[3,],
                          sigma = result[[1]][[num_data_set]]$theta.trace[1,],
                          rho = result[[1]][[num_data_set]]$theta.trace[2,])
  p1 = ggplot(data_trace, aes(x = MCMC.interactions))+
    geom_line(aes(x = MCMC.interactions,a))+theme_bw()+ylab("a")+xlab("MCMC iterations")

  p2 = ggplot(data_trace, aes(x = MCMC.interactions))+
    geom_line(aes(x = MCMC.interactions,b))+theme_bw()+ylab("b")+xlab("MCMC iterations")

  p3 = ggplot(data_trace, aes(x = MCMC.interactions))+
    geom_line(aes(x = MCMC.interactions,beta_0))+theme_bw()+ylab(TeX("$\\beta_0$"))+xlab("MCMC iterations")

  p4 = ggplot(data_trace, aes(x = MCMC.interactions))+
    geom_line(aes(x = MCMC.interactions,beta_1))+theme_bw()+xlab("MCMC iterations")+ylab(TeX("$\\beta_1$"))

  p5= ggplot(data_trace, aes(x = MCMC.interactions))+
    geom_line(aes(x = MCMC.interactions,beta_2))+theme_bw()+xlab("MCMC iterations")+ylab(TeX("$\\beta_2$"))

  p6= ggplot(data_trace, aes(x = MCMC.interactions))+
    geom_line(aes(x = MCMC.interactions,sigma))+theme_bw()+xlab("MCMC iterations")+ylab(TeX("$\\sigma$"))
  p7= ggplot(data_trace, aes(x = MCMC.interactions))+
    geom_line(aes(x = MCMC.interactions,rho))+theme_bw()+xlab("MCMC iterations")+ylab(TeX("$\\rho$"))

   grid.arrange(p1,p2,p3,p4,p5,p6,p7, ncol = 3)

}

#' Compute MCSE
#'
#' @param error_vec The vector composed of sum squared prediction error
#' @param nhours Number of simulation times
#' @param n.pred Number of unmonitored locations
#'
#' @return Monte Carlo Standard Error during simulation for nhours times

compute_mcse = function(error_vec, nhours = 500, n.pred = 100){
  mcse = sqrt(sum((na.omit(error_vec)[1:nhours]/n.pred-mean(na.omit(error_vec)[1:nhours]/n.pred))^2))/nhours
  return(mcse)
}


#' Combine results from 5 cases with different number of grid cells
#'
#'
#' @param resu_1 The first element from the output of the melding function when we consider 2 grid cells
#' @param resu_2 The first element from the output of the melding function when we consider 10 grid cells
#' @param resu_3 The first element from the output of the melding function when we consider 20 grid cells
#' @param resu_4 The first element from the output of the melding function when we consider 30 grid cells
#' @param resu_5 The first element from the output of the melding function when we consider 50 grid cells
#' @param nhours Number of simulation times
#' @param n.pred Number of unmonitored locations
#' @param nbeta Length of vector of regression parameter
#' @param ntheta Length of vector of covariance parameter
#' @param nab Length of vector of calibration parameter
#'
#' @return A list consisting of estimated parameter values and sum squared prediction error during simulation

#' @export
#'
#' @examples
compute_param_error = function(resu_1 ,resu_2,
                               resu_3 , resu_4 ,
                               resu_5,
                               nhours = 500, n.pred=100,nbeta=3,ntheta=2,nab=2){

  predictmat_2<-matrix(rep(NA,nhours*n.pred),ncol=nhours)
  betamat_2<-matrix(rep(NA,nbeta*nhours),ncol=nhours)
  betasdmat_2<-betamat_2
  thetamat_2<-matrix(rep(NA,ntheta*nhours),ncol=nhours)
  thetasdmat_2<-thetamat_2
  abmat_2<-matrix(rep(NA,nab*nhours),ncol=nhours)
  absdmat_2<-abmat_2

  predictmat_10<-matrix(rep(NA,nhours*n.pred),ncol=nhours)
  betamat_10<-matrix(rep(NA,nbeta*nhours),ncol=nhours)
  betasdmat_10<-betamat_2
  thetamat_10<-matrix(rep(NA,ntheta*nhours),ncol=nhours)
  thetasdmat_10<-thetamat_2
  abmat_10<-matrix(rep(NA,nab*nhours),ncol=nhours)
  absdmat_10<-abmat_2


  predictmat_20<-matrix(rep(NA,nhours*n.pred),ncol=nhours)
  betamat_20<-matrix(rep(NA,nbeta*nhours),ncol=nhours)
  betasdmat_20<-betamat_2
  thetamat_20<-matrix(rep(NA,ntheta*nhours),ncol=nhours)
  thetasdmat_20<-thetamat_2
  abmat_20<-matrix(rep(NA,nab*nhours),ncol=nhours)
  absdmat_20<-abmat_2

  predictmat_30<-matrix(rep(NA,nhours*n.pred),ncol=nhours)
  betamat_30<-matrix(rep(NA,nbeta*nhours),ncol=nhours)
  betasdmat_30<-betamat_2
  thetamat_30<-matrix(rep(NA,ntheta*nhours),ncol=nhours)
  thetasdmat_30<-thetamat_2
  abmat_30<-matrix(rep(NA,nab*nhours),ncol=nhours)
  absdmat_30<-abmat_2

  predictmat_50<-matrix(rep(NA,nhours*n.pred),ncol=nhours)
  betamat_50<-matrix(rep(NA,nbeta*nhours),ncol=nhours)
  betasdmat_50<-betamat_2
  thetamat_50<-matrix(rep(NA,ntheta*nhours),ncol=nhours)
  thetasdmat_50<-thetamat_2
  abmat_50<-matrix(rep(NA,nab*nhours),ncol=nhours)
  absdmat_50<-abmat_2

  error_2<-rep(NA,nhours)
  error_10<-rep(NA,nhours)
  error_20<-rep(NA,nhours)
  error_30<-rep(NA,nhours)
  error_50<-rep(NA,nhours)


  for(i in 1: nhours){
    predictmat_2[,i]<-resu_1[[1]][[i]]$prediction
    error_2[i]<-sum((resu_1[[2]]$true[,i]-predictmat_2[,i])^2)
    betamat_2[,i]<-resu_1[[1]][[i]]$beta.est
    betasdmat_2[,i]<-resu_1[[1]][[i]]$beta.est.sd
    thetamat_2[,i]<-resu_1[[1]][[i]]$theta.est
    thetasdmat_2[,i]<-resu_1[[1]][[i]]$theta.est.sd
    abmat_2[,i]<-resu_1[[1]][[i]]$ab.est
    absdmat_2[,i]<-resu_1[[1]][[i]]$ab.est.sd


    predictmat_10[,i]<-resu_2[[1]][[i]]$prediction
    error_10[i]<-sum((resu_2[[2]]$true[,i]-predictmat_10[,i])^2)
    betamat_10[,i]<-resu_2[[1]][[i]]$beta.est
    betasdmat_10[,i]<-resu_2[[1]][[i]]$beta.est.sd
    thetamat_10[,i]<-resu_2[[1]][[i]]$theta.est
    thetasdmat_10[,i]<-resu_2[[1]][[i]]$theta.est.sd
    abmat_10[,i]<-resu_2[[1]][[i]]$ab.est
    absdmat_10[,i]<-resu_2[[1]][[i]]$ab.est.sd


    predictmat_20[,i]<-resu_3[[1]][[i]]$prediction
    error_20[i]<-sum((resu_3[[2]]$true[,i]-predictmat_20[,i])^2)
    betamat_20[,i]<-resu_3[[1]][[i]]$beta.est
    betasdmat_20[,i]<-resu_3[[1]][[i]]$beta.est.sd
    thetamat_20[,i]<-resu_3[[1]][[i]]$theta.est
    thetasdmat_20[,i]<-resu_3[[1]][[i]]$theta.est.sd
    abmat_20[,i]<-resu_3[[1]][[i]]$ab.est
    absdmat_20[,i]<-resu_3[[1]][[i]]$ab.est.sd

    predictmat_30[,i]<-resu_4[[1]][[i]]$prediction
    error_30[i]<-sum((resu_4[[2]]$true[,i]-predictmat_30[,i])^2)
    betamat_30[,i]<-resu_4[[1]][[i]]$beta.est
    betasdmat_30[,i]<-resu_4[[1]][[i]]$beta.est.sd
    thetamat_30[,i]<-resu_4[[1]][[i]]$theta.est
    thetasdmat_30[,i]<-resu_4[[1]][[i]]$theta.est.sd
    abmat_30[,i]<-resu_4[[1]][[i]]$ab.est
    absdmat_30[,i]<-resu_4[[1]][[i]]$ab.est.sd

    predictmat_50[,i]<-resu_5[[1]][[i]]$prediction
    error_50[i]<-sum((resu_5[[2]]$true[,i]-predictmat_50[,i])^2)
    betamat_50[,i]<-resu_5[[1]][[i]]$beta.est
    betasdmat_50[,i]<-resu_5[[1]][[i]]$beta.est.sd
    thetamat_50[,i]<-resu_5[[1]][[i]]$theta.est
    thetasdmat_50[,i]<-resu_5[[1]][[i]]$theta.est.sd
    abmat_50[,i]<-resu_5[[1]][[i]]$ab.est
    absdmat_50[,i]<-resu_5[[1]][[i]]$ab.est.sd
  }

  return(list(error_2 = error_2,
              error_10 = error_10,
              error_20 = error_20,
              error_30 = error_30,
              error_50 = error_50,

              betamat_2 = betamat_2,
              betamat_10 = betamat_10,
              betamat_20 = betamat_20,
              betamat_30 = betamat_30,
              betamat_50 = betamat_50,

              betasdmat_2 = betasdmat_2,
              betasdmat_10 = betasdmat_10,
              betasdmat_20 = betasdmat_20,
              betasdmat_30 = betasdmat_30,
              betasdmat_50 = betasdmat_50,

              thetamat_2 = thetamat_2,
              thetamat_10 = thetamat_10,
              thetamat_20 = thetamat_20,
              thetamat_30 = thetamat_30,
              thetamat_50 = thetamat_50,

              thetasdmat_2 = thetasdmat_2,
              thetasdmat_10 = thetasdmat_10,
              thetasdmat_20 = thetasdmat_20,
              thetasdmat_30 = thetasdmat_30,
              thetasdmat_50 = thetasdmat_50,

              abmat_2 = abmat_2,
              abmat_10 = abmat_10,
              abmat_20 = abmat_20,
              abmat_30 = abmat_30,
              abmat_50 = abmat_50,

              absdmat_2 = absdmat_2,
              absdmat_10 = absdmat_10,
              absdmat_20 = absdmat_20,
              absdmat_30 = absdmat_30,
              absdmat_50 = absdmat_50))
}





#' Make Mean Squared Prediction Error table
#' @param nhours Number of simulation times
#' @param errorvec_1 The error vector from output of ``compute_param_error'' function with grid cell = 2
#' @param errorvec_2 The error vector from output of ``compute_param_error'' function with grid cell = 10
#' @param errorvec_3 The error vector from output of ``compute_param_error'' function with grid cell = 20
#' @param errorvec_4 The error vector from output of ``compute_param_error'' function with grid cell = 30
#' @param errorvec_5 The error vector from output of ``compute_param_error'' function with grid cell = 50
#' @param n.pred Number of unmonitored locations
#'
#' @return A dataframe summarizes the number of grid cells and average MSPE, MCSE
#' @export
#'
#' @examples mspe_tab(nhours = 500, errorvec_1 = compute_param_error()$error_2,errorvec_2 = compute_param_error()$error_10,errorvec_3 = compute_param_error()$error_20,errorvec_4 = compute_param_error()$error_30,errorvec_5 = compute_param_error()$error_50, n.pred = 100)
mspe_tab = function(nhours = 500, errorvec_1 = compute_param_error()$error_2,
                    errorvec_2 = compute_param_error()$error_10,
                    errorvec_3 = compute_param_error()$error_20,
                    errorvec_4 = compute_param_error()$error_30,
                    errorvec_5 = compute_param_error()$error_50, n.pred = 100){

  expo_summary = data.frame(Number.grid.cells = c(2, 10, 20, 30, 50),
                            MSPE = c(mean(errorvec_1/n.pred),
                                     mean(errorvec_2/n.pred),
                                     mean(errorvec_3/n.pred),
                                     mean(errorvec_4/n.pred),
                                     mean(errorvec_5/n.pred)),
                            MCSE = c(compute_mcse(errorvec_1, nhours, n.pred),
                                     compute_mcse(errorvec_2, nhours, n.pred),
                                     compute_mcse(errorvec_3, nhours, n.pred),
                                     compute_mcse(errorvec_4, nhours, n.pred),
                                     compute_mcse(errorvec_5, nhours, n.pred)))

 return(expo_summary)


}

#' Make Mean Squared Prediction Error Plot
#'
#' @param nhours Number of simulation times
#' @param errorvec_1 The error vector from output of ``compute_param_error'' function with grid cell = 2
#' @param errorvec_2 The error vector from output of ``compute_param_error'' function with grid cell = 10
#' @param errorvec_3 The error vector from output of ``compute_param_error'' function with grid cell = 20
#' @param errorvec_4 The error vector from output of ``compute_param_error'' function with grid cell = 30
#' @param errorvec_5 The error vector from output of ``compute_param_error'' function with grid cell = 50
#' @param n.pred Number of unmonitored locations
#'
#' @return The plot of the average of MSPE vs the number of grid cells
#' @export
#'
#' @examples mspe_tab(nhours = 500, errorvec_1 = compute_param_error()$error_2,errorvec_2 = compute_param_error()$error_10,errorvec_3 = compute_param_error()$error_20,errorvec_4 = compute_param_error()$error_30,errorvec_5 = compute_param_error()$error_50, n.pred = 100)
mspe_plt = function(nhours = 500,
                    errorvec_1 = compute_param_error()$error_2,
                    errorvec_2 = compute_param_error()$error_10,
                    errorvec_3 = compute_param_error()$error_20,
                    errorvec_4 = compute_param_error()$error_30,
                    errorvec_5 = compute_param_error()$error_50, n.pred = 100){

  expo_summary = data.frame(Number.grid.cells = c(2, 10, 20, 30, 50),
                            MSPE = c(mean(errorvec_1/n.pred),
                                     mean(errorvec_2/n.pred),
                                     mean(errorvec_3/n.pred),
                                     mean(errorvec_4/n.pred),
                                     mean(errorvec_5/n.pred)),
                            MCSE = c(compute_mcse(errorvec_1, nhours, n.pred),
                                     compute_mcse(errorvec_2, nhours, n.pred),
                                     compute_mcse(errorvec_3, nhours, n.pred),
                                     compute_mcse(errorvec_4, nhours, n.pred),
                                     compute_mcse(errorvec_5, nhours, n.pred)))
  ggplot(data = expo_summary) +
    geom_line(aes(x = Number.grid.cells,
                  y = MSPE))+theme_bw()+ xlab("Number of gird cells")
}

#' Generate plots to show the locations of monitored/unmonitored stations when grid cell = 2
#'
#' @param resu The output from the melding function when when grid cell = 2
#'
#' @return A ggplot with locations of monitored/unmonitored stations when grid cell = 2
#' @export
#'
#' @examples point_plt(compute_param_error())
point_plt = function(resu){
  ggplot()+
  geom_point(aes(x = resu[[2]]$sam.sloc[1:2,1], y = resu[[2]]$sam.sloc[1:2,2],
                 colour = "Sampled points from 2 grid cells",
                 size = "Sampled points from 2 grid cells"))+
  geom_point(aes(x = resu[[2]]$sam.sloc[3:22,1], y = resu[[2]]$sam.sloc[3:22,2],
                 colour = "Monitored locations",
                 size = "Monitored locations"))+
  geom_point(aes(x = resu[[2]]$sam.sloc1[,1], y = resu[[2]]$sam.sloc1[,2],
                 colour = "Unmonitored locations",
                 size = "Unmonitored locations"))+
  scale_colour_manual("",values = c("Sampled points from 2 grid cells" = "red",
                                 "Monitored locations" = "orange",
                                 "Unmonitored locations" = "black"))+
  scale_size_manual("",values = c("Sampled points from 2 grid cells" = 2,
                                 "Monitored locations" = 1.2,
                                 "Unmonitored locations" = 0.8))+
  theme(legend.position = "bottom")+
  xlim(-5,5)+ylim(-5,5)+theme_bw()+xlab("X-coordinate")+ylab("Y-coordinate")}



