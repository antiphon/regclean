#' Initialization of Z
#' @param nE total number of points
#' @param theta current vector of parameters
#' @param  lambda_1 current estiation of the intensity of the Strauss process
#' @param volAE volAE volume of full window
#' @param eps esp acceptance parameter for a configuration
#' @param DistancesBetweenPointsE  matrix with distances between points
#'
#' @export
MCMC_InitZ<-function(nE, theta, lambda_1, volAE, eps,DistancesBetweenPointsE)
  # simulate initial value of Z, we start with this initial value and then we update searching to make big posterior probabilities
{
  lambda_0 <- theta[1]
  
  #sample from the binomial distribution  
  n_0E<-rbinom(1,nE,lambda_0/(lambda_0+lambda_1))
  cont<-0
  # randomly pick n_0 noise points, value of the densities must be positive
  repeat
  {
    cont<-cont+1
    ZE<-rep(1,nE)
    ZE[sample(1:nE,n_0E,replace=F)] <- 0   
    
    #accept only configurations with density > 0      
    if(MCMC_DensityConfigurationGivenZ(ZE, theta, volAE, DistancesBetweenPointsE) > eps | cont > 100) {break} ##eps e grandissimo!!!
  }
  return(ZE)
} 
