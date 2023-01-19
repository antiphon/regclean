
#' Number of r-close couples in the realization of the Strauss process
#'
#' @param Z vector of 0 and 1, labelled with 1 are Strauss points
#' @param r r interaction radius
#' @param DistancesBetweenPoints matrix with distances between points 
#
#' @export  
MCMC_RClose <- function(Z, r,DistancesBetweenPoints)
{
  s <- length(which(DistancesBetweenPoints[which(Z==1),which(Z==1)]< r & DistancesBetweenPoints [which(Z==1),which(Z==1)] > 0))/2 # the second condition is for not counting the point itself
  return(s)
}