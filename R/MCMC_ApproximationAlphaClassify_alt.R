#'compute classification vector  MCMC with auxiliary method 1 means Strauss point 0 Poisson point
#'
#' @param r hardcore radius of the Strauss process
#' @param datappE simulated superposition, NEED SIMULATION IN THE WINDOW [-r,1+r]^2
#
#' @import coda spatstat
#' @export  
MCMC_ApproximationAlphaClassify_alt<-function(x,r,mcmcPar)
{
  cat("alt version\n")
  T0 <- Sys.time()
  # Check pattern
  x <- VB_check_pattern(x)
  # Check mcmcPar
  mcmcPar <- MCMC_check_mcmcPar(mcmcPar, x,r)
  #
  # Run estimation:
  rlist <- MCMC_computePosteriorsApproximationAlpha_alt(r, x, mcmcPar)
  #
  theta<-rlist$theta
  p<-rlist$p
  
  perc<-MCMC_computePercentage(theta,bboxE)
  
  ZEstimated<-MCMC_classificationFromPosteriors(p,perc,sim,r)
  #
  # compile:
  out <- list(ZEstimated=ZEstimated,
              p=p,
              theta=theta,              
              mcmcPar=mcmcPar,
              r=r,              
              data_facts = list(n=nrow(x$x), bbox=x$bbox),
              took=format(Sys.time()-T0))
  class(out) <- c("mcmcApproximationAlphaclassification", is(out))
  # Done!
  out
  
}
