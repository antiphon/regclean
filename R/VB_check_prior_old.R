#' Check priors for VBclassify
#' 
#' @param prior list of priors
#' @param x point pattern, list with $x and $bbox
#'
# #' @import Matrix
# #' @export
# VB_check_prior <- function(prior, x){
#   #'
#   n <- nrow(x$x)
#   V <- bbox_volume(x$bbox)
#   if(missing(prior)) prior <- NULL
#   if(is.null(prior)){
#     p0 <- 0.5
#     beta <- lambda <- n/V * p0
#     S <- Diagonal( x = c(1, 1, 1))
#     m <- log(c(lambda, beta, 0.2))
#     
#   }
#   else if(!is.list(prior)) stop("'prior' should be a list.")
#   else{
#     p0 <- if(is.null(prior$p0)) 0.5 else prior$p0
#     m <- if(is.null(prior$m)) { # mean
#             beta <- lambda <- n/V * p0
#             log(c(lambda, beta, 0.2))
#           }
#           else prior$m
#     S <- if(is.null(prior$S)){
#             Diagonal( x = c(1, 1, 1))
#           } else prior$S
#   }
#   list(m=m, S=S, p0=p0)
# }
