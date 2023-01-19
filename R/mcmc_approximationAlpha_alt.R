#' compute posteriors and estimation theta  , alt version
#' @import coda
#' @param r radius of hardcore process
#' @param reps total number of iterations of MCMC algorithm
#' @param repin burn-in iterations of MCMC algorithm
#' @param k gap between considered iterations in MCMC algorithm
#' @param DistancesBetweenPointsE matrix of the distances between points of the superposition
#' @param indConsidered subindices of points which belong to the restricted window
#' @param binConsidered vector of 0 and 1, 1 means that the point is inside the restricted window
#' @param volA volume of restricted window
#' @param volAE volume of big window
#' @param bbox restricted window in matrix form
#' @param supportPriorLeft vector of left limits prior support
#' @param supportPriorRight vector of right limits prior support
#
#' @export    
computePosteriorsMCMCApproximationAlpha_alt<-function(r, x, mcmcPar)
{ 
  bboxE <- x$bbox # original window
  bbox  <- bbox_dilate(bboxE, -r) # reduced window
  volAE <- bbox_volume(bboxE) # volume big window
  volA  <- bbox_volume(bbox) #volume restricted window     
  totE<-nrow(x$x) # total number of points in the big window  
  
  DistancesBetweenPointsE<-pairdist(x$x) #matrix of the distances between points
  bbdists <- bbox_distance(x$x, x$bbox)
  indConsidered <-  which(bbdists >= r) # within reduced window
  binConsidered<-rep(0, totE)
  binConsidered[indConsidered]<-1  # ind ios 1 if point is in restricted window otherwise 0  
  sim<-x$x[indConsidered,]# simulation in restricted window in matrix form  
  
  #'
  P <- mcmcPar
  #'
  #First guess of the parameters,  uniform number in the prior supports
  lambda_0  <-  runif(n=1, min=P$supportPriorLeft[1], max=P$supportPriorRight[1])
  beta  <- runif(n=1, min=P$supportPriorLeft[2], max=P$supportPriorRight[2])
  gamma <-  runif(n=1, min=P$supportPriorLeft[3], max=P$supportPriorRight[3])
  theta<- c(lambda_0, beta, gamma, r)
  #parameters for MCMC
  tau<-c(0.1,0.1,0.1)
  eps<-0.01  
  StG <- (1-theta[3])*pi*theta[4]^2 
  lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process   
  # first guess classification in big window and restricted window
  ZE <- Init( totE, theta, lambda_1, volAE, eps,DistancesBetweenPointsE)
  Z <- ZE[indConsidered]
  sumZ <- rep(0,length(indConsidered))
  
  #run MCMC method with approximation of the Strauss constant
  for(i in 1:repin)
  {
    StG <- (1-theta[3])*pi*theta[4]^2 
    lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process 
    ZE<-UpdateZ( ZE, theta,totE,DistancesBetweenPointsE,lambda_1);
    Z<-ZE[indConsidered];
    n1 <- length(which(Z==1)); # number of Strauss points
    n0<-length(Z)-n1; #number of Poisson points    
    theta <- UpdateThetaApproximationAlpha(ZE, theta, tau, supportPriorLeft, supportPriorRight,n0,n1,bbox,binConsidered,DistancesBetweenPointsE,StG,lambda_1,volA)  
  }
  
  for(i in (repin+1):reps)
  {
    StG <- (1-theta[3])*pi*theta[4]^2 
    lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process 
    ZE<-UpdateZ( ZE, theta,totE,DistancesBetweenPointsE,lambda_1);
    Z<-ZE[indConsidered]; # update parallely to Z also ZE 
    n1 <- length(which(Z==1)); # number of Strauss points
    n0<-length(Z)-n1;#number of Poisson points   
    theta <- UpdateThetaApproximationAlpha(ZE, theta, tau, supportPriorLeft, supportPriorRight,n0,n1,bbox,binConsidered,DistancesBetweenPointsE,StG,lambda_1,volA)  
    if(i%%k==0) # we consider values only every k iterations
    {
      sumZ<-sumZ+Z;   
      lambda_0<-append(lambda_0,theta[1]);
      beta<-append(beta,theta[2]);
      gamma<-append(gamma,theta[3]);
    }  
    
  }
  #posterior probabilities to be a Strauss point for every bubble in the restricted window
  p<-sumZ/((reps-repin)/k); 
  # estimation of theta from MCMC algorithms  
  chain<-mcmc(cbind(lambda_0,beta,gamma))
  summaryChain<-summary(chain)$statistics;
  theta<-c(summaryChain[1,1],summaryChain[2,1],summaryChain[3,1],r) 
  return(list(p=p,theta=theta,chain=chain));
} 


