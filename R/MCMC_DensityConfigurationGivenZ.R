#' computes p(y0,n0,y1,n1|theta)
#' @param Z vector of 0 and 1, labelled with 1 are Strauss points
#' @param ZrealE vector of 0 and 1, labelled with 1 are Strauss points in full window
#' @param theta current vector of parameters
#' @param volA volA volume of restricted window
#' @param DistancesBetweenPointE matrix with distances between points
#'
#' @export
MCMC_DensityConfigurationGivenZ<-function(ZE, theta, volAE,DistancesBetweenPointsE)
{
  
  ## we're evalueting p(y,Z| A,n,theta)=p(y0,n0,y1,n1| A,n,theta), for the indipendence of the two processes
  ## =p(y0,n0|A,n,theta)*p(y1,n1|A,n,theta)=density poisson*density strauss  
  
  n1E<-length(which(ZE == 1))#numeber of Strauss points in classification given by Z 
  n0E<-length(which(ZE == 0))#number of Poisson points in classification given by Z
  
  PoissonDensity <- exp(-(theta[1]-1)*volAE)*theta[1]^n0E  
  StraussDensity <- theta[2]^n1E * theta[3]^MCMC_RClose( ZE, theta[4],DistancesBetweenPointsE) #ignore normalizing constant
  if(is.nan(StraussDensity)) return(0)
  return(PoissonDensity * StraussDensity)
}
