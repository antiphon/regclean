#'compute classification vector  MCMC with auxiliary method 1 means Strauss point 0 Poisson point
#'
#' @param r hardcore radius of the Strauss process
#' @param x list  x=pattern to classify bbox=observation window
#' @param prior information about parameters priors
#
#' @import coda spatstat
#' @export  
MCMCClassifyGaussianPriors<-function(x,r,prior,reps,repin,k)
{ 
  
  T0 <- Sys.time()
  # Check pattern  
  x <- VB_check_pattern(x) 
  # Check priors
  prior <- VB_check_prior(prior, x) 
  bboxE<-x$bbox; # big window as a matrix
  bbox<-cbind(c((bboxE[1,1]+r),(bboxE[2,1]-r)),c((bboxE[1,2]+r),(bboxE[2,2]-r))) 
  volAE<-bbox_volume(bboxE);# volume big window
  volA<-bbox_volume(bbox) #volume restricted window     
  totE<-nrow(x$x) # total number of points in the big window  
  nrcloseNeighbors<-rep(0,totE);# vector i-component will contain number of r-close neighbor of the i-th point
  
  
  priorExpectationPercentageNoise<-exp(prior$m[1])/(totE/volAE)
  
  
  
  
  DistancesBetweenPointsE<-pairdist(x$x)#matrix of the distances between points 
  for(i in 1:totE)
  {
    nrcloseNeighbors[i]<-length(which(DistancesBetweenPointsE[i,]>0 & DistancesBetweenPointsE[i,]<r) )
  }    
  
  
  indConsidered<-which(x$x[,1]>=bbox[1,1] & x$x[,1]<=bbox[2,1] & x$x[,2]>=bbox[1,2] &  x$x[,2]<=bbox[2,2]) #index of points in restricted window
  binConsidered<-rep(0,totE);  
  binConsidered[indConsidered]<-1  # ind ios 1 if point is in restricted window otherwise 0  
  
  simE<-x$x; # simulation in big window   
 
  tau<-c(0.1,0.1,0.1)#mcmcPar$tau;
  eps<-10^(-3)#mcmcPar$eps;   
  meanPriors<-prior$m;
  varPriors<-diag(prior$S);
  if(details==TRUE)
  {
    MCMC_computePosteriorsGaussianPriors<-MCMC_computePosteriorsGaussianPriors_Detailed
    
  } 
  list<-MCMC_computePosteriorsGaussianPriors(r,reps,repin,k,DistancesBetweenPointsE,indConsidered,binConsidered,totE,volA,volAE,bbox,meanPriors, varPriors,tau,eps)#,nrcloseNeighbors)
  theta<-list$theta;
  p<-list$p;  
  chain<-list$chain;
  postConv<-list$postConv;
  listDetails<-list$listDetails;
  
  ## classification using stripes
  perc<-MCMC_computePercentage(as.vector(theta),cbind(c(0,1),c(0,1)));
  ris<-MCMC_classificationFromPosteriors(p,perc,DistancesBetweenPointsE,nrcloseNeighbors,priorExpectationPercentageNoise,r);
  ZEstimated<-ris$ZEstimated;
  indclusters<-ris$indclusters;
  percCouplesinMediumcluster<-ris$percCouplesinMediumcluster;
  indR<-ris$indR;
  
  
  
  
  ## classification is done using directly the posterior probabilities
  #ZEstimated<-rbinom(length(p),1,p) 
  
  
  
  out<-list(ZEstimated=ZEstimated,            
            p=p,
            theta=theta,
            chain=chain,
            postConv=postConv,
            x=x,
            prior=prior,
            r=r,
            percCouplesinMediumcluster=percCouplesinMediumcluster,
            indclusters=indclusters,
            indR=indR,
            listDetails=listDetails,
            took=format(Sys.time()-T0))
  
  class(out) <- c("mcmcClassifyGaussianPriorsclassification", is(out))
  # Done!
  out
  
}


## theta fixed

MCMCClassifyGaussianPriorsFixedTheta<-function(x,theta,lambda_1,reps,repin,k)
{ 
  ## since it is detailed version
#   sumAccZ<<-0
#   listn0E<<-vector()
#   hastingsRatioZ<<-vector()
  r<-theta[4]
  eps<-10^(-3)
  T0 <- Sys.time()
  # Check pattern  
  x <- VB_check_pattern(x) 
  # Check priors  
  bboxE<-x$bbox; # big window as a matrix
  bbox<-cbind(c((bboxE[1,1]+r),(bboxE[2,1]-r)),c((bboxE[1,2]+r),(bboxE[2,2]-r))) 
  volAE<-bbox_volume(bboxE);# volume big window
  volA<-bbox_volume(bbox) #volume restricted window     
  totE<-nrow(x$x) # total number of points in the big window  
  priorExpectationPercentageNoise<-theta[1]/(theta[1]+lambda_1) ### based on the true value of theta
  nrcloseNeighbors<-rep(0,totE);# vector i-component will contain number of r-close neighbor of the i-th point
  
  DistancesBetweenPointsE<-pairdist(x$x)#matrix of the distances between points 
   for(i in 1:totE)
   {
    nrcloseNeighbors[i]<-length(which(DistancesBetweenPointsE[i,]>0 & DistancesBetweenPointsE[i,]<r) )
   }    
  simE<-x$x; # simulation in big window   
  
  
  
  
  list<-MCMC_computePosteriorsGaussianPriorsFixedTheta(reps,repin,k,DistancesBetweenPointsE,totE,volA,volAE,bbox,theta,lambda_1,eps)#,nrcloseNeighbors)
  p<-list$p;
  postConv<-list$postConv;


  ## classification using stripes
  perc<-MCMC_computePercentage(as.vector(theta),cbind(c(0,1),c(0,1)));  
  ris<-MCMC_classificationFromPosteriors(p,perc,DistancesBetweenPointsE,nrcloseNeighbors,priorExpectationPercentageNoise,r); 
  ZEstimated<-ris$ZEstimated;
  indclusters<-ris$indclusters;
  percCouplesinMediumcluster<-ris$percCouplesinMediumcluster;
  indR<-ris$indR;

  #classification using directly the postesteriors
#   ZEstimated<-rbinom(length(p),1,p) 
  
  out<-list(p=p,  
            postConv=postConv,
            ZEstimated=ZEstimated,
            x=x,           
            took=format(Sys.time()-T0))
  
  class(out) <- c("mcmcClassifyGaussianPriorsclassificationFIxedTheta", is(out))
  #' Done!
  out
  
}


## fixing Z and estimate only the paraemters theta


MCMCClassifyGaussianPriorsFixedZ<-function(x,n0E,r,prior,tau)
{
  ## since it is detalied version
  listProposedlambda0<<-vector()
  listProposedbeta<<-vector()
  listProposedgamma<<-vector()
  
  hastingsRatiolambda0<<-vector()
  hastingsRatiobeta<<-vector()
  hastingsRatiogamma<<-vector()
  
  acclambda0<<-0
  accbeta<<-0
  accgamma<<-0
  
  
  T0 <- Sys.time()
  # Check pattern  
  x <- VB_check_pattern(x) 
  # Check priors
  prior <- VB_check_prior(prior, x)
  
  bboxE<-x$bbox; # big window as a matrix
  bbox<-cbind(c((bboxE[1,1]+r),(bboxE[2,1]-r)),c((bboxE[1,2]+r),(bboxE[2,2]-r))) 
  volAE<-bbox_volume(bboxE);# volume big window
  volA<-bbox_volume(bbox) #volume restricted window     
  totE<-nrow(x$x) # total number of points in the big window  
  
  
  DistancesBetweenPointsE<-pairdist(x$x)#matrix of the distances between points
  indConsidered<-which(x$x[,1]>=bbox[1,1] & x$x[,1]<=bbox[2,1] & x$x[,2]>=bbox[1,2] &  x$x[,2]<=bbox[2,2]) #index of points in restricted window
  binConsidered<-rep(0,totE);  
  binConsidered[indConsidered]<-1  # ind ios 1 if point is in restricted window otherwise 0  
  
  simE<-x$x; # simulation in big window   
  
  #remember u used 10^6 total iterations in set1
  reps<-10^4#mcmcPar$reps; # total number of iterations MCMC algorithm
  repin<-10^3#10^4#mcmcPar$repin # burn-in iterations 
  k<-10#10^2#mcmcPar$k; # interval between saved iterations  
  
  tau<-c(0.1,0.1,0.1)
  
  #tau<-c(0.3,0.3,0.3)#mcmcPar$tau;
  eps<-10^(-3)#mcmcPar$eps;   
  meanPriors<-prior$m;
  varPriors<-diag(prior$S);
  list<-MCMC_computePosteriorsGaussianPriorsFixedZ(r,reps,repin,k,DistancesBetweenPointsE,indConsidered,binConsidered,totE,volA,volAE,bbox,meanPriors, varPriors,tau,eps,n0E)
  
  theta<-list$theta; 
  chain<-list$chain;  
  
  out<-list(
    theta=theta,
    chain=chain,
    x=x,
    prior=prior,
    r=r,
    took=format(Sys.time()-T0))
  
  class(out) <- c("mcmcClassifyGaussianPriorsclassificationFixedZ", is(out))
  #' Done!
  out
  
}


## old version (??)
MCMCClassifyGaussianPriors_Old<-function(x,r,prior)
{
  acceptClassification<-FALSE 
  diffPercNoise<-0; 
  T0 <- Sys.time()
  #' Check pattern  
  x <- VB_check_pattern(x) 
  #' Check priors
  prior <- VB_check_prior(prior, x)
  case<-1
  bboxE<-x$bbox; # big window as a matrix
  bbox<-cbind(c((bboxE[1,1]+r),(bboxE[2,1]-r)),c((bboxE[1,2]+r),(bboxE[2,2]-r))) 
  volAE<-bbox_volume(bboxE);# volume big window
  volA<-bbox_volume(bbox) #volume restricted window     
  totE<-nrow(x$x) # total number of points in the big window  
  
  
  DistancesBetweenPointsE<-pairdist(x$x)#matrix of the distances between points
  indConsidered<-which(x$x[,1]>=bbox[1,1] & x$x[,1]<=bbox[2,1] & x$x[,2]>=bbox[1,2] &  x$x[,2]<=bbox[2,2]) #index of points in restricted window
  binConsidered<-rep(0,totE);  
  binConsidered[indConsidered]<-1  # ind ios 1 if point is in restricted window otherwise 0  
  
  simE<-x$x; # simulation in big window   
  
  #remember u used 10^6 total iterations in set1
  reps<-10^5#mcmcPar$reps; # total number of iterations MCMC algorithm
  repin<-10^4#10^4#mcmcPar$repin # burn-in iterations 
  k<-10^2#10^2#mcmcPar$k; # interval between saved iterations  
  
  tau<-c(0.1,0.1,0.1)#mcmcPar$tau;
  eps<-10^(-3)#mcmcPar$eps;   
  meanPriors<-prior$m;
  varPriors<-diag(prior$S);
  list<-MCMC_computePosteriorsGaussianPriors(r,reps,repin,k,DistancesBetweenPointsE,indConsidered,binConsidered,totE,volA,volAE,bbox,meanPriors, varPriors,tau,eps)
  
  theta<-list$theta;
  p<-list$p;  
  chain<-list$chain;
  
  perc<-MCMC_computePercentage(theta,bboxE);   
  percRealNoise<-exp(meanPriors[1])/(totE/volAE)
  #   class<-MCMC_classificationFromPosteriorsFinal(p,perc,simE,r,percRealNoise) 
  #  
  #   ZEstimated<-class$ZEstimated;
  #   
  #   if(ZEstimated != FALSE)
  #   { 
  #     acceptClassification<-TRUE
  #     percEstimatedNoise<-(length(which(ZEstimated==0)))/totE  
  #     diffPercNoise<-abs(percEstimatedNoise-percRealNoise)
  #     if(diffPercNoise>0.1)
  #     {
  #       warning("too much difference between real and estimated number of noise points, probably there is no structure in the posteriors or the clustering was done in a uncorrect way")
  #       acceptClassification<-FALSE
  #     }
  #     
  #     if(acceptClassification==FALSE)
  #     { 
  #       case<-2
  #       acceptClassification<-TRUE
  #       class<-MCMC_classificationFromPosteriorsFinal2(p,perc,simE,r,percRealNoise)
  #       ZEstimated<-class$ZEstimated;
  #       percEstimatedNoise<-(length(which(ZEstimated==0)))/totE  
  #       diffPercNoise<-abs(percEstimatedNoise-percRealNoise)
  #       if(diffPercNoise>0.1)
  #       {
  #         warning("too much difference between real and estimated number of noise points,also with two groups division")
  #         acceptClassification<-FALSE
  #       }
  #     }
  #   } 
  #   percCouplesinMediumcluster<-class$percCouplesinMediumcluster;
  #   indclusterHigh<-class$indclusterHigh;
  #   indclusterMedium<-class$indclusterMedium;
  #   indclusterLow<-class$indclusterLow;
  #   indCouples<-indclusterMedium[class$indR]
  #   #'
  #   #' compile:
  #   out <- list(ZEstimated=ZEstimated, 
  #               acceptClassification=acceptClassification,
  #               p=p,
  #               perc=perc,
  #               theta=theta,              
  #               prior=prior,
  #               r=r, 
  #               case=case,
  #               percCouplesinMediumcluster=percCouplesinMediumcluster,
  #               diffPercNoise=diffPercNoise,
  #               indclusterHigh=indclusterHigh,
  #               indclusterMedium=indclusterMedium,
  #               indclusterLow=indclusterLow,
  #               indCouples=indCouples,
  #               data_facts = list(n=nrow(x$x), bbox=x$bbox),
  #               took=format(Sys.time()-T0))
  
  out<-list(p=p,
            theta=theta,
            chain=chain,
            x=x,
            prior=prior,
            r=r,
            took=format(Sys.time()-T0))
  
  class(out) <- c("mcmcClassifyGaussianPriorsclassification", is(out))
  # Done!
  out
  
}
