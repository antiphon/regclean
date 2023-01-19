#' Update of theta in MCMC algorithm
#' 
#' @param ZE vector of classification in the big window
#' @param theta vector of actual parameters
#' @param tau variances for determine new proposals of parameters
#' @param supportPriorLeft left limits of prior support
#' @param supportPriorRight right limits of prior support
#' @param n0 total number of Poisson points
#' @param n1 total number of Strauss points
#' @param bbox retricted window in matrix form
#' @param binConsidered vector of 0 and 1, 1 means that the point is inside the restricted window
#' @param DistancesBetweenPointsE  matrix with distances between points
#' @param StG parameter for estimation of lambda_1
#' @param lambda_1 current estiation of the intensity of the Strauss process
#' @param volA volume of restricted window
#'
#' @import rstrauss
#' @export

MCMC_UpdateThetaApproximationAlpha_alt<-function(ZE, theta, tau, supportPriorLeft, supportPriorRight,n0,n1,bbox,binConsidered,Phi,StG,lambda_1,volA) 
{   
  M<-sample(1:3,1, replace=F)
  if (M == 1)
  {   
    lambda<- theta[1]
    t<-rnorm(n=1, mean=0, sd=tau[1])
    e<-exp(t);
    lambda_new<-lambda*e    
    if (lambda_new < supportPriorRight[1] && lambda_new > supportPriorLeft[1])
    {    
      new_term<-e^n0*((lambda+lambda_1)/(lambda_new+lambda_1))^(n0+n1);
      hastingsRatio<-e^(n0+1)*exp(-theta[1]*(e-1)*volA)*new_term;       
      if((runif(n=1)) < hastingsRatio)
      {      
        theta[1] <- lambda_new
      }
    }
  }
  
  if (M == 2)
  {  
    beta<- theta[2]
    t<-rnorm(n=1, mean=0, sd=tau[2])
    e<-exp(t); 
    beta_new<-beta*e;    
    lambda_1_new <- lambert_W0(beta_new*StG)/StG;
    if (beta_new < supportPriorRight[2] && beta_new > supportPriorLeft[2])
    { 
      new_term<-(lambda_1_new/lambda_1)^n1*((theta[1]+lambda_1)/(theta[1]+lambda_1_new))^(n0+n1);  
      hastingsRatio<-new_term*e^(n1+1)*(exp(approximate_strauss_constant(beta, theta[3], theta[4], bbox)-approximate_strauss_constant(beta_new, theta[3], theta[4], bbox)));
      if((runif(n=1)) < hastingsRatio)
      {       
        theta[2] <- beta_new
      }
    } 
  }
  
  if (M == 3)
  {
    gamma<- theta[3]
    t<-rnorm(n=1, mean=0, sd=tau[3])
    e<-exp(t);
    gamma_new<-gamma*e;
    StG_new <- (1-gamma_new)*pi*theta[4]^2; 
    lambda_1_new <- lambert_W0(theta[2]*StG_new)/StG_new;
    
    if (gamma_new < supportPriorRight[3] && gamma_new > supportPriorLeft[3]) {
      new_term<-(lambda_1_new/lambda_1)^n1*((theta[1]+lambda_1)/(theta[1]+lambda_1_new))^(n0+n1);   
      S <- sum(ZE*Phi%*%ZE)/2
      #S<-MCMC_RCloseE(ZE,binConsidered,theta[4],DistancesBetweenPointsE)
      hastingsRatio<-new_term*e^( S +1)*exp( approximate_strauss_constant(theta[2], gamma, theta[4], bbox)-
                                             approximate_strauss_constant(theta[2], gamma_new, theta[4], bbox) ) 
      
      if((runif(n=1)) < hastingsRatio)
      {  
        theta[3] <-  gamma_new
      }
    }
  }
  return(theta);
}

