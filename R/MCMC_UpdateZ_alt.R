
#' update Z
#' @param Z current vector of classification
#' @param theta current vector of parameters
#' @param n n total number of points in the restricted window
#' @param DistancesBetweenPoints  matrix with distances between points
#' @param lambda_1 estiamtion of intensity of the Strauss process 
#'
#' @export

MCMC_UpdateZ_alt<-function(Z, theta,n,Phi,lambda_1)
{
  n0<-sum(1-Z)#number of Poisson points in actual classification
  n1<-sum(Z)#number of Strauss points in actual classification
  
  # pick a random point and suggest to change its Z-value ( so if it was Strauss to Poisson and vice versa)
  m <- sample(1:n, 1)
  Z_new <- Z
  Z_new[m] <- if (Z[m]) 0 else 1
  
  #apply MCMC algorithm
  S <- sum(Phi[m,]*Z)
  if (Z[m] == 0 ) # if picked point was noise (Poisson)           
  {
    hastingsRatio <- (theta[2]*theta[3]^(S)*lambda_1)/(theta[1]^2)
  }
  else # if picked point was true (Strauss)
  {
    hastingsRatio <- (theta[1]^2)/(lambda_1* theta[2]*theta[3]^(S))    
  }
  
  if(runif(n=1) < hastingsRatio) # in this case we accept the proposal of new Z
  {
    return(Z_new);
  }
  else # Z doesn't change
  {    
    return(Z);    
  }
}
