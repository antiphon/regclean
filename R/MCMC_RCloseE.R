#' Number of r-close couples in the realization of the Strauss process where at least one point is in the restricted window
#'
#' @param Z vector of 0 and 1, labelled with 1 are Strauss points
#' @param binConsidered vector of 0 and 1, labelled with 1 are the points in the restricted window
#' @param r r interaction radius
#' @param DistancesBetweenPoints matrix with distances between all points 
#
#' @export  
MCMC_RCloseE  <-  function(Z,binConsidered, r,DistancesBetweenPoints)
{
  #number of pairs inside the considered window
  s <- length(which(DistancesBetweenPoints[which((Z==1)& (binConsidered==1)),which((Z==1)& (binConsidered==1))]< r & DistancesBetweenPoints [which((Z==1)& (binConsidered==1)),which((Z==1)& (binConsidered==1))] > 0))/2 
  #adding number of pairs where one point is inside considred window and one point outside 
  s=s+length(which(DistancesBetweenPoints[which((Z==1)& (binConsidered==0)),which((Z==1)& (binConsidered==1))] < r)) 
  return(s)
}