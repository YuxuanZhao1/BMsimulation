% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{mspe_tab}
\alias{mspe_tab}
\title{Make Average Mean Squared Prediction Error table}
\usage{
mspe_tab(
  nhours = 500,
  errorvec_1 = compute_param_error()$error_2,
  errorvec_2 = compute_param_error()$error_10,
  errorvec_3 = compute_param_error()$error_20,
  errorvec_4 = compute_param_error()$error_30,
  errorvec_5 = compute_param_error()$error_50,
  n.pred = 100
)
}
\arguments{
\item{nhours}{Number of simulation times}

\item{errorvec_1}{The error vector from output of ``compute_param_error'' function with grid cell = 2}

\item{errorvec_2}{The error vector from output of ``compute_param_error'' function with grid cell = 10}

\item{errorvec_3}{The error vector from output of ``compute_param_error'' function with grid cell = 20}

\item{errorvec_4}{The error vector from output of ``compute_param_error'' function with grid cell = 30}

\item{errorvec_5}{The error vector from output of ``compute_param_error'' function with grid cell = 50}

\item{n.pred}{Number of unmonitored locations}
}
\value{
A dataframe summarizes the number of grid cells and average MSPE, MCSE
}
\description{
Make Average Mean Squared Prediction Error table
}
\examples{
mspe_tab(nhours = 500, errorvec_1 = compute_param_error()$error_2,errorvec_2 = compute_param_error()$error_10,errorvec_3 = compute_param_error()$error_20,errorvec_4 = compute_param_error()$error_30,errorvec_5 = compute_param_error()$error_50, n.pred = 100)
}
