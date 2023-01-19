#'percentages of TT TF FF in r-close couples of a Strauss process with parameters theta  
#'
#' @import rstrauss
#' @param theta vector of parameters of the process
#' @param win window of observation
#
#' @export 

MCMC_computePercentage<-function(theta,bbox)
{
  percCouplesNoise<-0;
  percCouplesTrue<-0;  
  percGroupsTrue<-0;
  
  percCouplesMixed<-0;
  n<-100# number of simulations on the base of which percMixed is estimated... before it was 200 I put 100 for 250 points
  for(j in 1:n)
  { 
    nCouplesNoise<-0;
    nCouplesTrue<-0;
    nCouplesMixed<-0;     
    
    
    V <- prod( apply(bbox, 2, diff)  )
    StG<- (1-theta[3])*pi*theta[4]^2; 
    lambda1<- lambert_W0(theta[2]*StG)/StG;
    nNoisePoints<-ceiling(theta[1]*V);
    #########################nTruePoints<-round(lambda1*V);
    #simulation Strauss=true points: I'm using statspat function because percentage seemes to be little better
#     truePoints<-rStrauss(theta[2],theta[3],theta[4])
#     nTruePoints<-truePoints$n
#     truePoints<- cbind(truePoints$x,truePoints$y)
    
    ## simulation using Tuomas function, in the case we have 250 points simulation with CFTP from spatstat doesn't converge
    nTruePoints<-round(lambda1*V);
    truePoints<-rstrauss(n=nTruePoints,gamma=theta[3],range=theta[4],toroidal=TRUE,iter=2*10^4,bbox=cbind(c(0,1),c(0,1)))$x   
    
    #################truePoints<-(rstrauss(n=nTruePoints, gamma = theta[3], range = theta[4], bbox = bbox))$x; 
    #plot(truePoints)
    #simulation Poisson=noise points
    noisePoints<-cbind(runif( nNoisePoints,bbox[1,1],bbox[2,1]),runif( nNoisePoints,bbox[1,2],bbox[2,2]))
#     points(noisePoints,pch=20)
    #superposition   
    datapp<-rbind(truePoints, noisePoints)  
    ## info from simulation
    tot<-nTruePoints+nNoisePoints; # total number of points in window of interest
    Zreal<-c(rep(1,nTruePoints), rep(0,nNoisePoints));
    DistanceBetweenPoints<-pairdist(datapp)
    diag(DistanceBetweenPoints)<-Inf # in the minimal distances we don't want to consider point itself   
    
#     
#     n_inGroups<-0 # number of points that have more then one r-close neighbor 
#     cont<-0
#     ## computing percentages in group
#     for(i in 1:tot)
#     {  
#       ind<-which(DistanceBetweenPoints[i,]<theta[4])      
#       if(length(ind)>1)
#       {        
#         n_inGroups<- n_inGroups+length(ind)
#         for(j in 1:length(ind))
#         {
#           cont<-cont+Zreal[ind[j]]          
#         }          
#             
#       }     
#     }
#     percGroupsTrue<-percGroupsTrue+cont/n_inGroups;
#     if(n_inGroups==0){ percGroupsTrue<-0}
#     
    
    
    
    
    
    
    
    
    indR<-vector()
    indC<-vector()
    
    for(i in 1:tot)#loop on points
    {
      
      v<-which(DistanceBetweenPoints[i,]<theta[4])
      if((length(v)==1)&& length(which(DistanceBetweenPoints[v,]<theta[4]))==1 &&  which(DistanceBetweenPoints[v,]<theta[4])==i)
      {
#         points(datapp[i,1],datapp[i,2],col="red",cex=1.5)
#         points(datapp[v,1],datapp[v,2],col="red",cex=1.5)
        indR<-c(indR,i)
        indC<-c(indC,v)
        
      }
      
      
    } 
    if(length(indR>0))
    {   
      nCouples<-length(indR)/2;
      for(k in 1:nCouples) ## loop on the number of couples 
      {
        if(Zreal[indR][k]==0 & Zreal[indC][k]==0){nCouplesNoise=nCouplesNoise+1};
        if(Zreal[indR][k]==1 & Zreal[indC][k]==1){nCouplesTrue=nCouplesTrue+1};
        if((Zreal[indR][k]==1 & Zreal[indC][k]==0)|(Zreal[indR][k]==0 & Zreal[indC][k]==1)){nCouplesMixed=nCouplesMixed+1};  
      }
      
      totCouples<-nCouplesNoise+nCouplesMixed+nCouplesTrue # same as nCouples  
      percCouplesNoise<-percCouplesNoise+(nCouplesNoise/totCouples)   
      percCouplesMixed<-percCouplesMixed+(nCouplesMixed/totCouples)
      percCouplesTrue<-percCouplesTrue+(nCouplesTrue/totCouples)
      
    }
  }
  # mean of percentages in the n simulations
  percCouplesNoise<-percCouplesNoise/n;
  percCouplesMixed<-percCouplesMixed/n; 
  percCouplesTrue<-percCouplesTrue/n;
  #percGroupsTrue<-percGroupsTrue/n;
  return(list(percCouplesMixed=percCouplesMixed,percCouplesTrue=percCouplesTrue,percCouplesNoise=percCouplesNoise))#,percGroupsTrue=percGroupsTrue));
}
