#' Create the initial state object for VB strauss
#' 
#' @param quads number of quads, shortest edge. Longer edge count according to aspect ratio.
#' 
#' spam version
#' 
VB_initial_state <- function(x, R, prior, rho, eps, quads, keep.Q = FALSE, ninja = FALSE){
  V <- bbox_volume(x$bbox)
  n <- nrow(x$x)
  d <- ncol(x$x)
  #
  # Create dummy configuration
  rho <- if(missing(rho)) rep(n/V, 2) else rho
  if(length(rho)<2) rho <- c(rho,rho)[1:2]
  
  # sample: 
  el <- apply(x$bbox, 2, diff)
  n_w <- quads * el/min(el)  # grid
  # each window will hold
  wn <- round(rho * V / prod(n_w))
  if(any(wn<1)) stop("Reduce 'quads' or increase 'rho', automatically getting 0 dummy points per quandrant.")
  # 
  step <- el/n_w
  lefts <- lapply(1:d, function(i) seq(x$bbox[1,i], x$bbox[2,i]-step[i], length=n_w[i]))
  grid <- expand.grid(lefts)
  
  dummies0 <- do.call(rbind, lapply(split(grid, 1:nrow(grid)), 
                                    function(ab) apply(rbind(ab, ab+step), 2, 
                                                       function(g) runif(wn[1], g[1], g[2]) ) ))
  dummies1 <- do.call(rbind, lapply(split(grid, 1:nrow(grid)), 
                                    function(ab) apply(rbind(ab, ab+step), 2, 
                                                       function(g) runif(wn[2], g[1], g[2]) ) ))
  nd <- c(nrow(dummies0), nrow(dummies1))
  dummies <- rbind(dummies0, dummies1)
  
  #dummies0 <- apply(grid, 1, function(ab) apply(rbind(ab, ab+step), 2, function(g) runif(wn[1], g[1], g[2]) ) )
  #
  # Interaction:
  loc <- rbind(as.matrix(x$x), dummies)
  N <- nrow(loc)
  # no interaction from Poisson dummies to anywhere
  from <- setdiff(1:N, if(nd[1]>0)(n+1):(n+1 + nd[1]) else NULL ) 
  #
  # @TODO: Generalise this bit.
  #
  # Strauss:
  Phi <- nl2spam( geomgraph(loc, from=from, to = 1:n, r = R) )
  #  model in development
  if(ninja){
    D <- as.matrix(dist(loc))
    E <- (-(1-D/R)/(D/R))
    E[is.na(E)] <- 0
    diag(E) <- 0
    Phi <- Phi * as.spam(E)
  }
  #
  # auxiliary observations:
  y <- rep(c(1,0), c(n, sum(nd) ))
  #
  # Initial z:
  z <- rep(c(prior$p0, 0, 1), c(n, nd[1], nd[1]*0+nd[2]) )
  #
  # offset, will be recalculated later? (TODO)
  offset <- rep(-log(sum(rho)), length(z))    #-log(rho[1] * (1-z) + rho[2] * z)

  # Prior precalculated
  #p00 <- log(prior$p0 * V/(1-prior$p0))
  p00 <- log(prior$p0/(1-prior$p0))
  # starting values for xi:
  xi <- abs( offset )
  # compile
  list(Phi=Phi, 
       y=y, 
       z=z, 
       z_hist=z[1:n], 
       m=prior$m, 
       m_hist=prior$m,
       S=prior$S, 
       m0 = prior$m, 
       S0i=solve(prior$S), 
       Sm0 = (prior$m)%*%solve(prior$S),
       offset=offset,
       xi = xi,
       p00 = p00,
       eps=eps,
       rho=rho,
       verb=FALSE,
       L = NULL,
       Q = if(keep.Q) dummies else NULL)
}