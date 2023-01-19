
#'compute indices of clusters in which to divide the posterior probabilities
#'
#' @import cluster
#' @param p posterior probabilities to be a true bubble computed by MCMC
#' @param nrcloseNeighbors vector whose components contains the number of r-close neighbors of the points
#
#' @export  

MCMC_computeClusters<-function(p, nrcloseNeighbors)
{  
  k<-min((max(nrcloseNeighbors)+1),3); # number of clusters: maximum 3, depends on the maximum number of r-close neighbors
  if(k==1) # in this case we have a realisation of an hardcore process
  {
    return(list(indUp=1:length(p)))
  }
  fa <- fanny(p, k, metric="SqEuclidean")
  fannyc<-fa$cluster
  
  ind<-list() 
  nClusters<-max(fannyc)
  for( h in 1:nClusters )
  {
    ind[[h]]<-which(fannyc==h)
  }
  if(k==3)
  {
  ordInd<-sort(c(p[ind[[1]]][1],p[ind[[2]]][1],p[ind[[3]]][1]),index.return=TRUE)$ix
  return(list(indUp=ind[[ordInd[3]]],indMiddle=ind[[ordInd[2]]],indLow=ind[[ordInd[1]]]))
  }
  if(k==2)
  {
    ordInd<-sort(c(p[ind[[1]]][1],p[ind[[2]]][1]),index.return=TRUE)$ix
    return(list(indUp=ind[[ordInd[2]]],indMiddle=ind[[ordInd[1]]]))
  }
  
  
}