
#' Update of theta in MCMC algorithm
#' @import rstrauss
#' @param ZE vector of classification in the big window
#' @param theta vector of actual parameters
#' @param tau variances for determine new proposals of parameters
#' @param meanPriors mean of the prior lognormal distribution of the parameters
#' @param varPriors variance of the prior lognormal distribution of the parameters
#' @param n0 total number of Poisson points
#' @param n1 total number of Strauss points
#' @param bbox retricted window in matrix form
#' @param binConsidered vector of 0 and 1, 1 means that the point is inside the restricted window
#' @param DistancesBetweenPointsE  matrix with distances between points
#' @param StG parameter for estimation of lambda_1
#' @param lambda_1 current estiation of the intensity of the Strauss process
#' @param volA volume of restricted window
#'
#' @export
MCMC_UpdateThetaGaussianPriors<-function(ZE, theta, tau,meanPriors , varPriors,n0,n1,bbox,binConsidered,DistancesBetweenPointsE,StG,lambda_1,volA) 
{
  #choose randomly one of the parameters to update
  M<-sample(1:3,1, replace=F)  
  if (M == 1)
  {      
    lambda<- theta[1]
    t<-rnorm(n=1, mean=0, sd=tau[1])
    e<-exp(t);
    lambda_new<-lambda*e    
    ratio_likelihood<-e^(n0)*exp(-theta[1]*(e-1)*volA) 
    ratio_proposal<-e
    ratio_Z<-e^n0*((lambda+lambda_1)/(lambda_new+lambda_1))^(n0+n1); 
    ratio_priors<-(1/e)*exp(((-t)*(2*log(lambda)+t-2*meanPriors[1]))/(2*varPriors[1]))
    hastingsRatio<-ratio_likelihood*ratio_proposal*ratio_Z*ratio_priors;
    if((runif(n=1)) < hastingsRatio)
    {     
     theta[1] <- lambda_new
    }    
  }
  
  if (M == 2)
  {     
    beta<- theta[2]
    t<-rnorm(n=1, mean=0, sd=tau[2])
    e<-exp(t); 
    beta_new<-beta*e;   
    lambda_1_new <- lambert_W0(beta_new*StG)/StG;
    ratio_likelihood<-e^n1*(exp(approximate_strauss_constant(beta, theta[3], theta[4], bbox)-approximate_strauss_constant(beta_new, theta[3], theta[4], bbox)));
    ratio_proposal<-e
    ratio_Z<-(lambda_1_new/lambda_1)^n1*((theta[1]+lambda_1)/(theta[1]+lambda_1_new))^(n0+n1); 
    ratio_priors<-(1/e)*exp(((-t)*(2*log(beta)+t-2*meanPriors[2]))/(2*varPriors[2]))
    hastingsRatio<-ratio_likelihood*ratio_proposal*ratio_Z*ratio_priors;
    if((runif(n=1)) < hastingsRatio)
    {    
      theta[2] <- beta_new
    }
     
  }
  
  if (M == 3) ## proposal for gamma has different form with respect to the proposals for beta and lambda0
  {
    gamma<- theta[3]
    gamma_new<-min(rlnorm(n=1,meanlog =meanPriors[3] ,sdlog =varPriors[3]),0.99);
    e<-gamma_new/gamma;  
    t<-log(e);
    StG_new <- (1-gamma_new)*pi*theta[4]^2; 
    lambda_1_new <- lambert_W0(theta[2]*StG_new)/StG_new;
    ratio_likelihood<-e^( MCMC_RCloseE(ZE,binConsidered,theta[4],DistancesBetweenPointsE))*(exp(approximate_strauss_constant(theta[2], gamma, theta[4], bbox)-approximate_strauss_constant(theta[2], gamma_new, theta[4], bbox)));
    ratio_Z<-(lambda_1_new/lambda_1)^n1*((theta[1]+lambda_1)/(theta[1]+lambda_1_new))^(n0+n1);   
    hastingsRatio<-ratio_likelihood*ratio_Z; 
    
  if((runif(n=1)) < hastingsRatio)
    {  
      theta[3] <-  gamma_new
    }
    
  }

### update of gamma with the same form as for beta and lambda_0  
  
#   if (M == 3)
#   {    
#     gamma<- theta[3]
#     t<-rnorm(n=1, mean=0, sd=tau[3])
#     e<-exp(t);
#     gamma_new<-gamma*e;   
#     if(gamma_new>=1) gamma_new=0.99;
#   
#     StG_new <- (1-gamma_new)*pi*theta[4]^2; 
#     lambda_1_new <- lambert_W0(theta[2]*StG_new)/StG_new;
#     
#     ratio_likelihood<-e^( MCMC_RCloseE(ZE,binConsidered,theta[4],DistancesBetweenPointsE))*(exp(approximate_strauss_constant(theta[2], gamma, theta[4], bbox)-approximate_strauss_constant(theta[2], gamma_new, theta[4], bbox)));
#     ratio_proposal<-e    
#     ratio_Z<-(lambda_1_new/lambda_1)^n1*((theta[1]+lambda_1)/(theta[1]+lambda_1_new))^(n0+n1);   
#     ratio_priors<-(1/e)*exp(((-t)*(2*log(gamma)+t-2*meanPriors[3]))/(2*varPriors[3]))    
#     hastingsRatio<-ratio_likelihood*ratio_proposal*ratio_Z*ratio_priors;  
#     if((runif(n=1)) < hastingsRatio)
#     {      
#       theta[3] <-  gamma_new
#     }
#     
#   }
  return(theta);
}




## detailed version

MCMC_UpdateThetaGaussianPriors_Detailed<-function(ZE, theta, tau,meanPriors , varPriors,n0,n1,bbox,binConsidered,DistancesBetweenPointsE,StG,lambda_1,volA) 
{   
  M<-sample(1:3,1, replace=F)  
  if (M == 1)
  {      
    lambda<- theta[1]
    t<-rnorm(n=1, mean=0, sd=tau[1])
    e<-exp(t);
    lambda_new<-lambda*e  
    listProposedlambda0<<-c(listProposedlambda0,lambda_new)
    ratio_likelihood<-e^(n0)*exp(-theta[1]*(e-1)*volA) 
    ratio_proposal<-e^2 ## still to check wether is correct
    ratio_Z<-e^n0*((lambda+lambda_1)/(lambda_new+lambda_1))^(n0+n1); 
    ratio_priors<-(1/e)*exp(((-t)*(2*log(lambda)+t-2*meanPriors[1]))/(2*varPriors[1]))
    hastingsRatio<-ratio_likelihood*ratio_proposal*ratio_Z*ratio_priors;
    hastingsRatiolambda0<<-c(hastingsRatiolambda0,hastingsRatio)
    if((runif(n=1)) < hastingsRatio)
    {  
      acclambda0<<-acclambda0+1
      theta[1] <- lambda_new
    }    
  }
  
  if (M == 2)
  {     
    beta<- theta[2]
    t<-rnorm(n=1, mean=0, sd=tau[2])
    e<-exp(t); 
    beta_new<-beta*e; 
    listProposedbeta<<-c(listProposedbeta,beta_new)
    lambda_1_new <- lambert_W0(beta_new*StG)/StG;
    ratio_likelihood<-e^n1*(exp(approximate_strauss_constant(beta, theta[3], theta[4], bbox)-approximate_strauss_constant(beta_new, theta[3], theta[4], bbox)));
    ratio_proposal<-e^2 ## still to check wether is correct
    ratio_Z<-(lambda_1_new/lambda_1)^n1*((theta[1]+lambda_1)/(theta[1]+lambda_1_new))^(n0+n1); 
    ratio_priors<-(1/e)*exp(((-t)*(2*log(beta)+t-2*meanPriors[2]))/(2*varPriors[2]))
    hastingsRatio<-ratio_likelihood*ratio_proposal*ratio_Z*ratio_priors;
    hastingsRatiobeta<<-c(hastingsRatiobeta,hastingsRatio)
    if((runif(n=1)) < hastingsRatio)
    {   
      accbeta<<-accbeta+1
      theta[2] <- beta_new
    }
    
  }
  
  if (M == 3) ## proposal for gamma has different form with respect to the proposals for beta and lambda0
  {
    gamma<- theta[3]
    gamma_new<-min(rlnorm(n=1,meanlog =meanPriors[3] ,sdlog =varPriors[3]),0.99);
    listProposedgamma<<-c(listProposedgamma,gamma_new)
    e<-gamma_new/gamma;  
    t<-log(e);
    StG_new <- (1-gamma_new)*pi*theta[4]^2; 
    lambda_1_new <- lambert_W0(theta[2]*StG_new)/StG_new;
    ratio_likelihood<-e^( MCMC_RCloseE(ZE,binConsidered,theta[4],DistancesBetweenPointsE))*(exp(approximate_strauss_constant(theta[2], gamma, theta[4], bbox)-approximate_strauss_constant(theta[2], gamma_new, theta[4], bbox)));
    ratio_Z<-(lambda_1_new/lambda_1)^n1*((theta[1]+lambda_1)/(theta[1]+lambda_1_new))^(n0+n1);   
    hastingsRatio<-ratio_likelihood*ratio_Z;   
    hastingsRatiogamma<<-c(hastingsRatiogamma,hastingsRatio)
    if((runif(n=1)) < hastingsRatio)
    {  
      accgamma<<-accgamma+1
      theta[3] <-  gamma_new
    }
    
  }
  
  ### update of gamma with the same form as for beta and lambda_0  
  
#     if (M == 3)
#     {      
#       gamma<- theta[3]
#       t<-rnorm(n=1, mean=0, sd=tau[3])
#       e<-exp(t);
#       gamma_new<-gamma*e;
#       listProposedgamma<<-c(listProposedgamma,gamma_new)
#       if(gamma_new>=1) gamma_new=0.99;
#     
#       StG_new <- (1-gamma_new)*pi*theta[4]^2; 
#       lambda_1_new <- lambert_W0(theta[2]*StG_new)/StG_new;
#       
#       ratio_likelihood<-e^( MCMC_RCloseE(ZE,binConsidered,theta[4],DistancesBetweenPointsE))*(exp(approximate_strauss_constant(theta[2], gamma, theta[4], bbox)-approximate_strauss_constant(theta[2], gamma_new, theta[4], bbox)));
#       ratio_proposal<-e      
#       ratio_Z<-(lambda_1_new/lambda_1)^n1*((theta[1]+lambda_1)/(theta[1]+lambda_1_new))^(n0+n1);      
#       ratio_priors<-(1/e)*exp(((-t)*(2*log(gamma)+t-2*meanPriors[3]))/(2*varPriors[3]))     
#       hastingsRatio<-ratio_likelihood*ratio_proposal*ratio_Z*ratio_priors;
#       hastingsRatiogamma<<-c(hastingsRatiogamma,hastingsRatio)
#       if((runif(n=1)) < hastingsRatio)
#       {         
#         theta[3] <-  gamma_new
#       }
#       
#     }
  return(theta);
}





