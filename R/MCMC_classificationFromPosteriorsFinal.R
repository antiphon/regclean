#'classification from posteriors and some subproducts
#'
#' @param p posterior probability of points to be Strauss points
#' @param perc vector of probabilities to give to T-T T-F and F-F in couples
#' @param DistancesBetweenPointsE ditances between points of superposition pattern
#' @param nrcloseNeighbors, number of r-close neighbors for every point
#' @param priorExpectationPercentageNoise [TODO]
#
#' @import spatstat
#' @export      


MCMC_classificationFromPosteriors<-function(p,perc,DistancesBetweenPointsE,nrcloseNeighbors,priorExpectationPercentageNoise,r)
{
  lengthp<-length(p)# number of points we classify
  ZEstimated<-rep(0,lengthp); # vector which will contain our classification ( at starting all points are classified as noise)  
#   if(length(which(p==1.0))>=0.8*(length(p))) 
#   {
#     warning("need more iterations in MCMC algorithm since entries are more than 80 % equal to 1")   
#     return( list(ZEstimated=FALSE,percCouplesinMediumcluster=0,indclusterHigh=vector(),indclusterMedium=vector(),indclusterLow=vector(),indR=vector()))
#     
#   } 
  
  ##
  percCouplesinMediumcluster<-0
  indR<-NULL
  indclusters<-MCMC_computeClusters(p, nrcloseNeighbors)
  k<-length(indclusters)# number of clusters
  
  # up cluster we classify as true(note possibility to compute percentages also in this case)
  ZEstimated[indclusters$indUp]<-1  
if(k==1) {return(list(ZEstimated=ZEstimated,percCouplesinMediumcluster=NULL,indclusters=NULL,indR=NULL))}
  
  #points which corresponds to the medium cluster 
if(length(indclusters$indMiddle)>1)
{
  D<-DistancesBetweenPointsE[indclusters$indMiddle,indclusters$indMiddle]
  diag(D)<-Inf
  tot<-length(indclusters$indMiddle)

  #find the (real) couples inside the medium cluster and make for them classification using compute percentage
  indR<-vector() #indices points of couples
  indC<-vector() #respective indices ( so we know which point is with which)
  for(i in 1:tot)#loop on points of medium cluster 
  {
    
    v<-which(D[i,]<r)
    if((length(v)==1)&& length(which(D[v,]<r)==1))# &&  which(D[v,]<r)==i) )#I think this last part is already implied by first one
    {
      indR<-c(indR,i)
      indC<-c(indC,v)
      
    }  
    
  } 
  ## computing percentage of medium cluster that is covered by couples
  percCouplesinMediumcluster<-0
  if(length(indR>0))
  {   
    nCouples<-length(indR)/2;
    percCouplesinMediumcluster<-(nCouples*2)/tot;  
  }  
  
  if(percCouplesinMediumcluster>0)
  {       
    for(k in 1:(nCouples*2)) ## loop on the number of couples 
    {
      ra1<-runif(1)
      if(ra1<perc$percCouplesMixed) # one bubble at random classified as true other noise
      {
        ra2<-runif(1);
        if(ra2<0.5)
        {ZEstimated[indclusters$indMiddle[indR[k]]]<-1 ;}
        else {ZEstimated[indclusters$indMiddle[indC[k]]]<-1};
      }
      
      if((ra1>=perc$percCouplesMixed) & ra1<(perc$percCouplesMixed + perc$percCouplesTrue)) # both bubbles classified as true
      {
        ZEstimated[indclusters$indMiddle[indR[k]]]<-1
        ZEstimated[indclusters$indMiddle[indC[k]]]<-1
      }    
      
    }
  }
  ### not couples , most probably they belong to groups 
#  ZEstimated[indclusters$indMiddle[-indR]]<-rbinom(length(indclusters$indMiddle[-indR]),1,perc$percGroupsTrue)
 # classify teh bubbles in intermediate stripe that are not in couples
#  if(percCouplesinMediumcluster<1 )     
#  {
  ZEstimated[indclusters$indMiddle[-indR]]<-rbinom(length(indclusters$indMiddle[-indR]),1,(1-priorExpectationPercentageNoise)) # notice that if no couples are found in this way all the bubbles in the middle cluster are classified as true
#  }
}
  ## classify lower
  ZEstimated[indclusters$indLow]<-0#rbinom(length(indclusters$indLow),1,perc$percGroupsTrue)


 return( list(ZEstimated=ZEstimated,percCouplesinMediumcluster=percCouplesinMediumcluster,indclusters=indclusters,indR=indR))
  
}