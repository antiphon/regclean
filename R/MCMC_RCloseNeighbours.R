
#' number of r-close neighbours of the m-th  point in the Strauss process
#' @param Z vector of 0 and 1, labelled with 1 are Strauss points
#' @param r interaction radius
#' @param m m identifies one particular point (the m-th point)
#' @param DistancesBetweenPoints  matrix with distances between points
#'
#' @export

MCMC_RCloseNeighbours <- function(Z, r, m,DistancesBetweenPoints)
{
  return ((length(which(DistancesBetweenPoints[m,which(Z==1)] <r & DistancesBetweenPoints[m,which(Z==1)] >0))))
}
