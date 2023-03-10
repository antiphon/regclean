% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VBclassify.R
\name{VBclassify}
\alias{VBclassify}
\title{Noise removal from regular point patterns}
\usage{
VBclassify(
  x,
  R,
  prior,
  eps = 0.01,
  rho,
  verb = FALSE,
  quads = 1,
  ninja = FALSE
)
}
\arguments{
\item{x}{Point pattern, as list with $x coordinate matrix $bbox bound box ranges (col-wise) (or a ppp)}

\item{R}{Range for Strauss model.}

\item{prior}{Priors, see Details.}

\item{eps}{Convergence criterion.}

\item{rho}{Intensity of dummies, default = 2*n/Volume. Two dummy sets will be generated.}

\item{verb}{Verbose output?}

\item{quads}{How many subwindows to use for dummy sampling stratification? Default is 1, no stratification.}

\item{ninja}{Instead of Strauss potential gamma^1(d < R), use a diminishing potential without extra parameter: (d/R)^1(d < R)}
}
\description{
Base process assumed to be the Strauss model, noise process assumed Poisson, processes independent.
}
\details{
Only rectangle windows supported at the moment.

The priors: It should be a list.
}
\references{
[1] A. Baddeley, J. Coeurjolly, E. Rubak, and R. Waagepetersen, "Logistic regression for spatial Gibbs point processes," Biometrika, 2013.
}
