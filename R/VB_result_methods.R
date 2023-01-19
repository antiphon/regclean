#' Print method for VB classification
#'
#' @param x object from VBclassify
#'
#' @export
print.VBclassification <- function(x, ...) {
  cat("Result of noise classification using VB\nCall: ")
  print(x$call)
}

#' Plot method for VB classification
#' 
#' Plots the z-trajectories.
#' @export
plot.VBclassification <- function(x, ...) {
  it <- nrow(x$z_hist)
  plot(NA, xlab="iteration", ylab="P(z=1)", ylim=0:1, xlim=c(0,it), ...)
  for(i in 1:ncol(x$z_hist))
    lines(x$z_hist[,i], col=i)
}





#' Summary method for VB classification
#' 
#' @export
summary.VBclassification <- function(x){
  
  post <- x$posterior
  pri <- x$prior
  s <- sqrt(diag(post$S))
  sp <- sqrt(diag(pri$S))
  m <- post$m[,1]
  mp <- pri$m
  nam <- names(x$coef)
  #est <- data.frame(mean=m, sd=s, ci5=m-1.95*s,ci95=m+1.96*s, m.prior=mp, sd.prior=sp)
  #rownames(est) <- nam
  est <- m
  names(est) <- nam
  n <- x$data_facts$n
  o<-list(ndata = n,
          W = x$data_facts$bbox,
          nest_true = sum(x$pred),
          noise_prop = round(100*sum(1-x$pred)/n, 3),
          estimates = est,
          coef = x$coef,
          pred=x$pred,
          prob=x$prob,
          x=x)
  class(o) <- "VBclassification_summary"
  o
}

#' Print the summary of VB classification result
#' 
#' @param x VB classification result
#' 
#' @export
print.VBclassification_summary <- function(x, ...){
  print(x$x)
  bb <- x$W
  w <- paste(apply(bb, 2, function(v)paste0("[",v[1],",",v[2],"]")), collapse=" x ")
  cat("Pattern:",x$ndata, "points in",w," \n")
  cat("Estimated true points:", x$nest_true,"\n")
  cat("Noise proportion:", x$noise_prop, "%\n")
  cat("R given:", x$x$R, "\n")
  cat("\nMean estimates in typical parametrization:\n")
  print(x$coef)
  
  # estimates
  cat("\nEstimates in exponential-form parametrization:\n")
  print(x$estimates)
  #
  # Diagnostics
  x<-x$x
  it <- nrow(x$z_hist)-1
  cat("\nConvergence eps:", x$eps)
  cat("\nIterations:", it)
  cat("\nTime used:", x$took)
  cat("\n")
}






