
#' compute posteriors and estimation theta  
#'
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
MCMC_computePosteriorsApproximationAlpha_alt<-function(r, x, mcmcPar)
{ 
  
  bboxE<-x$bbox; # big window as a matrix
  bbox<-cbind(c((bboxE[1,1]+r),(bboxE[2,1]-r)),c((bboxE[1,2]+r),(bboxE[2,2]-r))) 
  volAE<-bbox_volume(bboxE);# volume big window
  volA<-bbox_volume(bbox) #volume restricted window     
  totE<-nrow(x$x) # total number of points in the big window  
  
  
  DistancesBetweenPointsE<-pairdist(x$x)#matrix of the distances between points
  indConsidered<-which(x$x[,1]>=bbox[1,1] & x$x[,1]<=bbox[2,1] & x$x[,2]>=bbox[1,2] &  x$x[,2]<=bbox[2,2]) #index of points in restricted window
  binConsidered<-rep(0,totE);
  binConsidered[indConsidered]<-1  # ind ios 1 if point is in restricted window otherwise 0  
  sim<-x$x[indConsidered,]# simulation in restricted window in matrix form  
  
  reps<-mcmcPar$reps; # total number of iterations MCMC algorithm
  repin<-mcmcPar$repin # burn-in iterations 
  k<-mcmcPar$k; # interval between saved iterations  
  supportPriorLeft<-mcmcPar$supportPriorLeft
  supportPriorRight<-mcmcPar$supportPriorRight
  
  tau<-mcmcPar$tau;
  eps<-mcmcPar$eps;  
  
  # Interaction matrix:
  Phi <- nl2spam(geomgraph(x$x, from=indConsidered, r=r))
  #First guess of the parameters,  uniform number in the prior supports
  lambda_0  <-  runif(n=1, min=supportPriorLeft[1], max=supportPriorRight[1]);
  beta  <- runif(n=1, min=supportPriorLeft[2], max=supportPriorRight[2]);
  gamma <-  runif(n=1, min=supportPriorLeft[3], max=supportPriorRight[3]);
  theta<- c(lambda_0, beta, gamma, r);  
  #parameters for MCMC  
  StG <- (1-theta[3])*pi*theta[4]^2 
  lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process   
  # first guess classification in big window and restricted window
  ZE<-MCMC_InitZ( totE, theta, lambda_1, volAE, eps,DistancesBetweenPointsE);  
  Z<-ZE[indConsidered]; 
  sumZ<-rep(0,length(indConsidered));
  
  #run MCMC method with approximation of the Strauss constant
  cat("Start:\n")
  for(i in 1:repin)
  {
    StG <- (1-theta[3])*pi*theta[4]^2 
    lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process 
    ZE<-MCMC_UpdateZ_alt( ZE, theta,totE, Phi, lambda_1)
    Z<-ZE[indConsidered];
    n1 <- length(which(Z==1)); # number of Strauss points
    n0<-length(Z)-n1; #number of Poisson points    
    theta <- MCMC_UpdateThetaApproximationAlpha_alt(ZE, theta, tau, supportPriorLeft, supportPriorRight,n0,n1,bbox,binConsidered,Phi,StG,lambda_1,volA)  
    #cat("  \r", i,"/", reps)
  }
  
  for(i in (repin+1):reps)
  {
    StG <- (1-theta[3])*pi*theta[4]^2 
    lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process 
    ZE<-MCMC_UpdateZ_alt( ZE, theta,totE,Phi,lambda_1);
    Z<-ZE[indConsidered]; # update parallely to Z also ZE 
    n1 <- length(which(Z==1)); # number of Strauss points
    n0<-length(Z)-n1;#number of Poisson points   
    theta <- MCMC_UpdateThetaApproximationAlpha_alt(ZE, theta, tau, supportPriorLeft, supportPriorRight,n0,n1,bbox,binConsidered,Phi,StG,lambda_1,volA)  
    if(i%%k==0) # we consider values only every k iterations
    {
      sumZ<-sumZ+Z;   
      lambda_0<-append(lambda_0,theta[1]);
      beta<-append(beta,theta[2]);
      gamma<-append(gamma,theta[3]);
    }  
    #cat("  \r", i,"/", reps)
  }
  #cat("\n")
  #posterior probabilities to be a Strauss point for every bubble in the restricted window
  p<-sumZ/((reps-repin)/k); 
  # estimation of theta from MCMC algorithms  
  #plot(p)  
  summaryChain<-summary(mcmc(cbind(lambda_0,beta,gamma)))$statistics;
  theta<-c(summaryChain[1,1],summaryChain[2,1],summaryChain[3,1],r) 
  cat("\n")
  return(list(p=p,theta=theta));
} 
