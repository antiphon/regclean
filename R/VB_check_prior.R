#' Check priors for VBclassify
#' 
#' @param prior list of priors
#' @param x point pattern, list with $x and $bbox
#'
#' @details 
#' 
#' Priors:
#' p0 ~ scalar[0,1], prior for P(z_i=1) for all i. Default: 0.5
#' m  ~ numeric(3),  mean for theta=log(lambda,beta,gamma) 
#' S  ~ Covariance matrix fof theta
#'
#' Default prior variance is 10000, flat priors.

#' @import spam
#' @export

VB_check_prior <- function(prior, x, ninja = FALSE){
  #'
  n <- nrow(x$x)
  V <- bbox_volume(x$bbox)
  if(missing(prior)) prior <- NULL
  if(is.null(prior)){
    p0 <- 0.5
    beta <- lambda <- n/V * p0
    S <- spam::diag.spam( x = c(1, 1, 1)*1e4)
    m <- log(c(lambda, beta, 0.2))
  }
  else if(!is.list(prior)) stop("'prior' should be a list.")
  else{
    p0 <- if(is.null(prior$p0)) 0.5 else prior$p0
    m <- if(is.null(prior$m)) { # mean
            beta <- lambda <- n/V * p0
            log(c(lambda, beta, 0.2))
          }
          else prior$m
    S <- if(is.null(prior$S)){
            spam::diag.spam( x = c(1, 1, 1)*1e4 )
          } else prior$S
  }
  if(ninja){
    S[3,3] <- 1e-9
    m[3] <- 1
  }
  list(m=m, S=S, p0=p0)
}
