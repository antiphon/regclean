% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMC_UpdateZ.R
\name{MCMC_UpdateZ}
\alias{MCMC_UpdateZ}
\title{update Z}
\usage{
MCMC_UpdateZ(Z, theta, n, DistancesBetweenPoints, lambda_1, m)
}
\arguments{
\item{Z}{current vector of classification}

\item{theta}{current vector of parameters}

\item{n}{n total number of points in the restricted window}

\item{DistancesBetweenPoints}{matrix with distances between points}

\item{lambda_1}{estiamtion of intensity of the Strauss process}
}
\description{
update Z
}
