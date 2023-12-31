% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simudata.R
\name{simulate_data}
\alias{simulate_data}
\title{Simulate realizations at locations of monitoring/unmonitored stations and within grid cells to construct 'nhours' independent data sets}
\usage{
simulate_data(m1, m2, n.pred, ran, nm, nhours, seed = 1234)
}
\arguments{
\item{m1}{Number of monitoring locations}

\item{m2}{Number of grid cells}

\item{n.pred}{Number of unmonitored locations for prediction}

\item{ran}{Range for the monitoring/unmonitored locations and sampling points within grid cells. Lower bound is -ran and upper bound is ran}

\item{nm}{Number of points with in each of grid cells}

\item{nhours}{Number of independent data sets during simulation}

\item{seed}{Random seed number}
}
\value{
A list including different relizations: Z(B), hat(Z)_u, Z_g, and locations
}
\description{
Simulate realizations at locations of monitoring/unmonitored stations and within grid cells to construct 'nhours' independent data sets
}
\examples{
simulate_data(20, 2, 100, 5, 1, 500, 1234)
}
