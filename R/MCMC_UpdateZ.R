
#' update Z
#' @param Z current vector of classification
#' @param theta current vector of parameters
#' @param n n total number of points in the restricted window
#' @param DistancesBetweenPoints  matrix with distances between points
#' @param lambda_1 estiamtion of intensity of the Strauss process 
#'
#' @export

MCMC_UpdateZ<-function(Z, theta,n,DistancesBetweenPoints,lambda_1,m)
{
  n0<-length(which(Z==0));#number of Poisson points in actual classification
  n1<-length(which(Z==1));#number of Strauss points in actual classification
  
  # pick a random point and suggest to change its Z-value ( so if it was Strauss to Poisson and vice versa)
  Z_new <- Z
  Z_new[m] <- if (Z[m]) 0 else 1
  
  #apply MCMC algorithm
  if (Z[m] == 0 ) # if picked point was noise (Poisson)           
  {
    hastingsRatio <- (theta[2]*theta[3]^(MCMC_RCloseNeighbours( Z, theta[4], m,DistancesBetweenPoints))*lambda_1)/(theta[1]^2)
  }
  else # if picked point was true (Strauss)
  {
    hastingsRatio <- (theta[1]^2)/(lambda_1* theta[2]*theta[3]^(MCMC_RCloseNeighbours(Z_new, theta[4], m,DistancesBetweenPoints)))
  }
  
  if(runif(n=1) < hastingsRatio) # in this case we accept the proposal of new Z
  {    
    return(ZE=Z_new);
  }
  else # Z doesn't change
  {    
    return(ZE=Z);    
  }
}







MCMC_UpdateZ_Detailed<-function(Z, theta,n,DistancesBetweenPoints,lambda_1,m)
{
  n0<-length(which(Z==0));#number of Poisson points in actual classification
  n1<-length(which(Z==1));#number of Strauss points in actual classification
  
  # pick a random point and suggest to change its Z-value ( so if it was Strauss to Poisson and vice versa)
  Z_new <- Z
  Z_new[m] <- if (Z[m]) 0 else 1
  
  #apply MCMC algorithm
  if (Z[m] == 0 ) # if picked point was noise (Poisson)           
  {
    hastingsRatio <- (theta[2]*theta[3]^(MCMC_RCloseNeighbours( Z, theta[4], m,DistancesBetweenPoints))*lambda_1)/(theta[1]^2)    
    hastingsRatioZ<<-c(hastingsRatioZ,hastingsRatio)
  }
  else # if picked point was true (Strauss)
  {
    hastingsRatio <- (theta[1]^2)/(lambda_1* theta[2]*theta[3]^(MCMC_RCloseNeighbours(Z_new, theta[4], m,DistancesBetweenPoints)))
    hastingsRatioZ<<-c(hastingsRatioZ,hastingsRatio)
  }
  
  if(runif(n=1) < hastingsRatio) # in this case we accept the proposal of new Z
  {
    return(list(ZE=Z_new,accZ=1));
  }
  else # Z doesn't change
  {    
    return(list(ZE=Z,accZ=0));    
  }
}











