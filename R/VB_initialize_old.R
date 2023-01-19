#' Create the initial state object for VB strauss
# #' @import cutgeom
# 
# VB_initial_state <- function(x, R, prior, rho, eps){
#   V <- bbox_volume(x$bbox)
#   n <- nrow(x$x)
#   #' Create dummy configuration
#   rho <- if(missing(rho)) 2*n/V else rho
#   nd <- 2 * round( rho * V )
#   dummies <- apply(x$bbox, 2, function(ab) runif(nd, ab[1], ab[2]) )
#   #'
#   #' Interaction:
#   loc <- rbind(as.matrix(x$x), dummies)
#   Phi <- nl2Matrix( geom(loc, to = 1:n, r = R) )
#   #' auxiliary observations:
#   y <- rep(c(1,0), c(n, nd))
#   #'
#   #' offset
#   offset <- rep(-log(rho), n+nd)
#   #' Initial z:
#   z <- rep(c(prior$p0, 0, 1), c(n, nd/2, nd/2) )
#   #'
#   #' Prior precalculated
#   p00 <- log(prior$p0 * V/(1-prior$p0))
#   #' compile
#   
#   list(Phi=Phi, 
#        y=y, 
#        z=z, 
#        z_hist=z[1:n], 
#        m=prior$m, 
#        m_hist=prior$m,
#        S=prior$S, 
#        m0 = prior$m, 
#        S0i=solve(prior$S), 
#        Sm0 = (prior$m)%*%solve(prior$S),
#        offset=offset,
#        xi = rep(2, n+nd),
#        p00 = p00,
#        eps=eps,
#        verb=FALSE,
#        L = NULL)
# }