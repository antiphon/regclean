
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
#' @param meanPriors mean of the normal in the lognormal prior distributions 
#' @param varPriors variance of the normal in the lognormal prior distributions
#
#' @export    
MCMC_computePosteriorsGaussianPriors<-function(r,reps,repin,k,DistancesBetweenPointsE,indConsidered,binConsidered,totE,volA,volAE,bbox, meanPriors,varPriors,tau,eps)#,nrcloseNeighbors)
{ 
  ## If we want to make an update of Z which depends on the number of r-close neighbors
#   unique_nrcloseNeighbors<-unique(nrcloseNeighbors) # vector of possible number of r-close neighbors
#   indrcloseNeighbors<-vector(mode="list",length=length(unique_nrcloseNeighbors)) # list of vectors vector i contains indices of point which have i-1 r-close neighbor
#   for ( i in 1: length(unique_nrcloseNeighbors))
#   {
#     indrcloseNeighbors[[i]]<-which(nrcloseNeighbors==unique_nrcloseNeighbors[i])      
#   }
 
  #First guess of the parameters,  uniform number in the prior supports 
  lambda_0<-max(rlnorm(n=1,meanlog =meanPriors[1] ,sdlog =varPriors[1]),0.1)  
  beta  <- max(rlnorm(n=1,meanlog =meanPriors[2] ,sdlog =varPriors[2]),0.1)  
  gamma <-  min(rlnorm(n=1,meanlog =meanPriors[3] ,sdlog =varPriors[3]),0.99)    
  theta<- c(lambda_0, beta, gamma, r);  
 
  StG <- (1-theta[3])*pi*theta[4]^2 
  lambda_1 <- lambert_W0(theta[2]*StG)/StG; #estimation of the intensity of the Strauss process 
  
  # first guess classification in big window
  ZE<-MCMC_InitZ( totE, theta, lambda_1, volAE, eps,DistancesBetweenPointsE); 
  
  # to track the convergence of the posterior probablities
  postConv<-list() 
  cont<-0
 
  
  n0<-length(which(ZE==0))
  sumZE<-rep(0,totE);
  
  # MCMC ALGORITHM
  
  for(i in 1:repin)
  {   
    StG <- (1-theta[3])*pi*theta[4]^2 
    lambda_1 <- lambert_W0(theta[2]*StG)/StG; #estimation of the intensity of the Strauss process 
  
    ## update pf Z using the number of neighbors 
#     h<-sample(1:length(unique_nrcloseNeighbors),1) ## choose randomly how many r-close points the point that we update should have
#     ord<-sample(indrcloseNeighbors[[h]],l,replace=FALSE) # choose randomly one of teh point with k-1 r-close neighbors 

    #choose one point of Z to update
    ord<-sample(1:length(ZE),1,replace=FALSE)
    ZE<-MCMC_UpdateZ( ZE, theta,totE,DistancesBetweenPointsE,lambda_1,ord);
    Z<-ZE[indConsidered];
    n1 <- length(which(Z==1)); # number of Strauss points
    n0<-length(Z)-n1; #number of Poisson points
    theta <- MCMC_UpdateThetaGaussianPriors(ZE, theta, tau, meanPriors ,varPriors,n0,n1,bbox,binConsidered,DistancesBetweenPointsE,StG,lambda_1,volA)  
  }
  
  for(i in (repin+1):reps)
  {
    StG <- (1-theta[3])*pi*theta[4]^2 
    lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process 
    ## update pf Z using the number of neighbors 
#     h<-sample(1:length(unique_nrcloseNeighbors),1) ## choose randomly how many r-close points the point that we update should have
#     ord<-sample(indrcloseNeighbors[[h]],l,replace=FALSE) # choose randomly one of teh point with k-1 r-close neighbors 
   
    ord<-sample(1:length(ZE),1,replace=FALSE)
    ZE<-MCMC_UpdateZ( ZE, theta,totE,DistancesBetweenPointsE,lambda_1,ord);
    Z<-ZE[indConsidered]
    n1 <- length(which(Z==1)); # number of Strauss points
    n0<-length(Z)-n1;#number of Poisson points   
    theta <- MCMC_UpdateThetaGaussianPriors(ZE, theta, tau, meanPriors ,varPriors,n0,n1,bbox,binConsidered,DistancesBetweenPointsE,StG,lambda_1,volA)   
    if(i%%k==0) # we consider values only every k iterations
    {
      sumZE<-sumZE+ZE;   
      # add this line to check the convergence of the posterior probablities
      cont<-cont+1
      postConv[[cont]]<-sumZE/cont
      lambda_0<-append(lambda_0,theta[1]);
      beta<-append(beta,theta[2]);
      gamma<-append(gamma,theta[3]);
    }  
    
  }
  #posterior probabilities to be a Strauss point for every bubble in the restricted window
  p<-sumZE/((reps-repin)/k); 
  chain<-mcmc(cbind(lambda_0,beta,gamma))
  summaryChain<-summary(chain)$statistics;
  theta<-c(summaryChain[1,1],summaryChain[2,1],summaryChain[3,1],r) 
  return(list(p=p,theta=theta,chain=chain,postConv=postConv));
} 


## version with theta fixed 

MCMC_computePosteriorsGaussianPriorsFixedTheta<-function(reps,repin,k,DistancesBetweenPointsE,totE,volA,volAE,bbox,theta,lambda_1,eps)#,nrcloseNeighbors)
{ 
#   maxnrcloseNeighbors<-max(nrcloseNeighbors) # max number of r-close points in the superposition
#   indrcloseNeighbors<-vector(mode="list",length=(maxnrcloseNeighbors+1)) # list of vectors vector i contains indices of point which have i-1 r-close neighbor
#   for ( i in 1: (maxnrcloseNeighbors+1 ))
#   {
#     indrcloseNeighbors[[i]]<-which(nrcloseNeighbors==(i-1))      
#   }
  # first guess classification in big window and restricted window
  ZE<-MCMC_InitZ( totE, theta, lambda_1, volAE, eps,DistancesBetweenPointsE); 
 
  l<-1#ceiling((totE/100)*40)
  sumZE<-rep(0,totE); 
  # to track the convergence of the posterior probablities
  postConv<-list() 
  cont<-0
  #MCMC 
  for(i in 1:repin)
  {  
#     h<-sample(1:(maxnrcloseNeighbors+1),1) ## choose randomly how many r-close points the point that we update should have
#     ord<-sample(indrcloseNeighbors[[h]],l,replace=FALSE) # choose randomly one of teh point with k-1 r-close neighbors 
    ord<-sample(1:length(ZE),1,replace=FALSE)
    ZE<-MCMC_UpdateZ( ZE, theta,totE,DistancesBetweenPointsE,lambda_1,ord);
#     Z<-ZE[indConsidered];
#     n1 <- length(which(Z==1)); # number of Strauss points
#     n0<-length(Z)-n1;#number of Poisson points   
  }
  
  for(i in (repin+1):reps)
  {    
#     h<-sample(1:(maxnrcloseNeighbors+1),1) ## choose randomly how many r-close points the point that we update should have
#     ord<-sample(indrcloseNeighbors[[h]],l,replace=FALSE) # choose randomly one of teh point with k-1 r-close neighbors  
#     ord<-sample(indrcloseNeighbors[[h]],l,replace=FALSE) # choose randomly one of teh point with k-1 r-close neighbors 
    ord<-sample(1:length(ZE),1,replace=FALSE)
    ZE<-MCMC_UpdateZ( ZE, theta,totE,DistancesBetweenPointsE,lambda_1,ord);
#     Z<-ZE[indConsidered];
#     n1 <- length(which(Z==1)); # number of Strauss points
#     n0<-length(Z)-n1;#number of Poisson points 
    if(i%%k==0) # we consider values only every k iterations
    {
      sumZE<-sumZE+ZE;   
      cont<-cont+1
      postConv[[cont]]<-sumZE/cont
    }  
    
  }
  #posterior probabilities to be a Strauss point for every bubble in the restricted window
  p<-sumZE/((reps-repin)/k); 
  return(list(p=p,postConv=postConv));
} 




## considering Z known and estimating only the parameters

MCMC_computePosteriorsGaussianPriorsFixedZ<-function(r,reps,repin,k,DistancesBetweenPointsE,indConsidered,binConsidered,totE,volA,volAE,bbox,meanPriors, varPriors,tau,eps,n0E,theta)
{ 
  #First guess of the parameters,  uniform number in the prior supports
  lambda_0  <-  max(rlnorm(n=1,meanlog =meanPriors[1] ,sdlog =varPriors[1]),0.1)
  beta  <- max(rlnorm(n=1,meanlog =meanPriors[2] ,sdlog =varPriors[2]),0.1)
  gamma <-  min(rlnorm(n=1,meanlog =meanPriors[3] ,sdlog =varPriors[3]),0.99)
  theta<- c(lambda_0, beta, gamma, r);    
  ZE<-c(rep(1,totE-n0E),rep(0,n0E))
  Z<-ZE[indConsidered];
  n1 <- length(which(Z==1)); # number of Strauss points
  n0<-length(Z)-n1; #number of Poisson points
  
  #run MCMC method with approximation of the Strauss constant
  for(i in 1:repin)
  {    
    StG <- (1-theta[3])*pi*theta[4]^2 
    lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process     
    theta <- MCMC_UpdateThetaGaussianPriors_Detailed(ZE, theta, tau, meanPriors ,varPriors,n0,n1,bbox,binConsidered,DistancesBetweenPointsE,StG,lambda_1,volA)  
  }
  
  for(i in (repin+1):reps)
  {    
    StG <- (1-theta[3])*pi*theta[4]^2 
    lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process     
    theta <- MCMC_UpdateThetaGaussianPriors_Detailed(ZE, theta, tau, meanPriors ,varPriors,n0,n1,bbox,binConsidered,DistancesBetweenPointsE,StG,lambda_1,volA)   
    if(i%%k==0) # we consider values only every k iterations
    {
#       print(i)
#       cat(" ")      
      lambda_0<-append(lambda_0,theta[1]);
      beta<-append(beta,theta[2]);
      gamma<-append(gamma,theta[3]);
    }  
    
  }
  
  chain<-mcmc(cbind(lambda_0,beta,gamma))
  summaryChain<-summary(chain)$statistics;
  theta<-c(summaryChain[1,1],summaryChain[2,1],summaryChain[3,1],r) 
  return(list(theta=theta,chain=chain,postConv=postConv));
} 







## Detailed version


MCMC_computePosteriorsGaussianPriors_Detailed<-function(r,reps,repin,k,DistancesBetweenPointsE,indConsidered,binConsidered,totE,volA,volAE,bbox, meanPriors,varPriors,tau,eps)#,nrcloseNeighbors)
{ 
  ## details we're going to compute
  listProposedlambda0<<-vector()
  listProposedbeta<<-vector()
  listProposedgamma<<-vector()
  hastingsRatiolambda0<<-vector()
  hastingsRatiobeta<<-vector()
  hastingsRatiogamma<<-vector()
  acclambda0<<-0
  accbeta<<-0
  accgamma<<-0
  sumAccZ<<-0
  hastingsRatioZ<<-vector()
  
  
  
#   maxnrcloseNeighbors<-max(nrcloseNeighbors) # max number of r-close points in the superposition
#   indrcloseNeighbors<-vector(mode="list",length=(maxnrcloseNeighbors+1)) # list of vectors vector i contains indices of point which have i-1 r-close neighbor
#   for ( i in 1: (maxnrcloseNeighbors+1 ))
#   {
#     indrcloseNeighbors[[i]]<-which(nrcloseNeighbors==(i-1))   
#     
#   }
  #First guess of the parameters,  uniform number in the prior supports
  lambda_0  <-  max(rlnorm(n=1,meanlog =meanPriors[1] ,sdlog =varPriors[1]),0.1)
  beta  <- max(rlnorm(n=1,meanlog =meanPriors[2] ,sdlog =varPriors[2]),0.1)
  gamma <-  min(rlnorm(n=1,meanlog =meanPriors[3] ,sdlog =0.000000001),0.99)
  theta<- c(lambda_0, beta, gamma, r);  
  #parameters for MCMC  
  StG <- (1-theta[3])*pi*theta[4]^2 
  lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process   
  # first guess classification in big window and restricted window
  ZE<-MCMC_InitZ( totE, theta, lambda_1, volAE, eps,DistancesBetweenPointsE);   
  n0<-length(which(ZE==0))
  sumZE<-rep(0,totE);
  #l<-1
  # MCMC
  for(i in 1:repin)
  {   
    StG <- (1-theta[3])*pi*theta[4]^2 
    lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process 
    ## update pf Z using the number of neighbors 
#     h<-sample(1:(maxnrcloseNeighbors+1),1) ## choose randomly how many r-close points the point that we update should have
#     ord<-sample(indrcloseNeighbors[[h]],l,replace=FALSE) # choose randomly one of teh point with k-1 r-close neighbors 

    ord<-sample(1:length(ZE),1,replace=FALSE)
    risE<-MCMC_UpdateZ_Detailed( ZE, theta,totE,DistancesBetweenPointsE,lambda_1,ord);
    ZE<-risE$ZE;    
    sumAccZ<-sumAccZ+risE$accZ;    
    Z<-ZE[indConsidered];
    n1 <- length(which(Z==1)); # number of Strauss points
    n0<-length(Z)-n1; #number of Poisson points    
    theta <- MCMC_UpdateThetaGaussianPriors_Detailed(ZE, theta, tau, meanPriors ,varPriors,n0,n1,bbox,binConsidered,DistancesBetweenPointsE,StG,lambda_1,volA)  
  }
  
  for(i in (repin+1):reps)
  {
    StG <- (1-theta[3])*pi*theta[4]^2 
    lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process 
#     h<-sample(1:(maxnrcloseNeighbors+1),1) ## choose randomly how many r-close points the point that we update should have
#     ord<-sample(indrcloseNeighbors[[h]],l,replace=FALSE) # choose randomly one of teh point with k-1 r-close neighbors 

    ord<-sample(1:length(ZE),1,replace=FALSE)
    risE<-MCMC_UpdateZ_Detailed( ZE, theta,totE,DistancesBetweenPointsE,lambda_1,ord);
    ZE<-risE$ZE;
    sumAccZ<-sumAccZ+risE$accZ;
    Z<-ZE[indConsidered]
    n1 <- length(which(Z==1)); # number of Strauss points
    n0<-length(Z)-n1;#number of Poisson points   
    theta <- MCMC_UpdateThetaGaussianPriors_Detailed(ZE, theta, tau, meanPriors ,varPriors,n0,n1,bbox,binConsidered,DistancesBetweenPointsE,StG,lambda_1,volA)   
    if(i%%k==0) # we consider values only every k iterations
    {
      sumZE<-sumZE+ZE;   
      lambda_0<-append(lambda_0,theta[1]);
      beta<-append(beta,theta[2]);
      gamma<-append(gamma,theta[3]);
    }  
    
  }
  #posterior probabilities to be a Strauss point for every bubble in the restricted window
  p<-sumZE/((reps-repin)/k); 
  # estimation of theta from MCMC algorithms  
  chain<-mcmc(cbind(lambda_0,beta,gamma))
  summaryChain<-summary(chain)$statistics;
  theta<-c(summaryChain[1,1],summaryChain[2,1],summaryChain[3,1],r) 
  
  ## in this case we also return list with all details
  listDetails<-list(listProposedlambda0=listProposedlambda0,
                    listProposedbeta=listProposedbeta,
                    listProposedgamma=listProposedgamma,
                    hastingsRatiolambda0=hastingsRatiolambda0,
                    hastingsRatiobeta=hastingsRatiobeta,
                    hastingsRatiogamma=hastingsRatiogamma,
                    acclambda0=acclambda0,
                    accbeta=accbeta,
                    accgamma=accgamma,
                    sumAccZ=sumAccZ,
                    hastingsRatioZ=hastingsRatioZ)
  
  return(list(p=p,theta=theta,chain=chain,listDetails=listDetails));
} 

## old version (update of Z simply takes one random point)



MCMC_computePosteriorsGaussianPriors_Old<-function(r,reps,repin,k,DistancesBetweenPointsE,indConsidered,binConsidered,totE,volA,volAE,bbox, meanPriors,varPriors,tau,eps)
{ 
  #First guess of the parameters,  uniform number in the prior supports
  lambda_0  <-  max(rlnorm(n=1,meanlog =meanPriors[1] ,sdlog =varPriors[1]),0.1)
  beta  <- max(rlnorm(n=1,meanlog =meanPriors[2] ,sdlog =varPriors[2]),0.1)
  gamma <-  min(rlnorm(n=1,meanlog =meanPriors[3] ,sdlog =varPriors[3]),0.99)
  theta<- c(lambda_0, beta, gamma, r);  
  #parameters for MCMC  
  StG <- (1-theta[3])*pi*theta[4]^2 
  lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process   
  # first guess classification in big window and restricted window
  ZE<-MCMC_InitZ( totE, theta, lambda_1, volAE, eps,DistancesBetweenPointsE); 
  l<-1#ceiling((totE/100)*40)
  sumZE<-rep(0,totE);
  
  #run MCMC method with approximation of the Strauss constant
  for(i in 1:repin)
  {
    StG <- (1-theta[3])*pi*theta[4]^2 
    lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process 
    ord<-sample(1:totE,l,replace=FALSE)
    ZE<-MCMC_UpdateZ( ZE, theta,totE,DistancesBetweenPointsE,lambda_1,ord);
    Z<-ZE[indConsidered];
    n1 <- length(which(Z==1)); # number of Strauss points
    n0<-length(Z)-n1; #number of Poisson points 
    theta <- MCMC_UpdateThetaGaussianPriors(ZE, theta, tau, meanPriors ,varPriors,n0,n1,bbox,binConsidered,DistancesBetweenPointsE,StG,lambda_1,volA)  
  }
  
  for(i in (repin+1):reps)
  {   
    StG <- (1-theta[3])*pi*theta[4]^2 
    lambda_1 <- lambert_W0(theta[2]*StG)/StG;#estimation of the intensity of the Strauss process 
    ord<-sample(1:totE,l,replace=FALSE)
    ZE<-MCMC_UpdateZ( ZE, theta,totE,DistancesBetweenPointsE,lambda_1,ord);
    Z<-ZE[indConsidered]
    n1 <- length(which(Z==1)); # number of Strauss points
    n0<-length(Z)-n1;#number of Poisson points   
    listn0<<-c(listn0,n0)
    theta <- MCMC_UpdateThetaGaussianPriors(ZE, theta, tau, meanPriors ,varPriors,n0,n1,bbox,binConsidered,DistancesBetweenPointsE,StG,lambda_1,volA)   
    if(i%%k==0) # we consider values only every k iterations
    {
#       print(i)
#       cat(" ")
      sumZE<-sumZE+ZE;   
      lambda_0<-append(lambda_0,theta[1]);
      beta<-append(beta,theta[2]);
      gamma<-append(gamma,theta[3]);
    }  
    
  }
  #posterior probabilities to be a Strauss point for every bubble in the restricted window
  p<-sumZE/((reps-repin)/k); 
  # estimation of theta from MCMC algorithms  
  #plot(p)  
  chain<-mcmc(cbind(lambda_0,beta,gamma))
  summaryChain<-summary(chain)$statistics;
  theta<-c(summaryChain[1,1],summaryChain[2,1],summaryChain[3,1],r) 

  return(list(p=p,theta=theta,chain=chain));
} 





