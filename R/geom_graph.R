#' Geometric Graph
#' 
#' @param x matrix
#' @param from indices
#' @param to indices
#' @param r radius
#'
#' @useDynLib regclean
geomgraph <- function(x, from, to, r){
  x <- as.matrix(x)
  if(missing(from)) from <- 1:nrow(x)
  if(missing(to)) to <- 1:nrow(x)
  if(missing(r)) stop("r needed")
  c_geom(x, from, to, r)
}

#' Transform neighbourhood list to a Matrix
#' 
#' @param x Neighbourhood list as given by \code{\link{geom}}
#' 
#' @import spam
nl2spam <- function(x){
  n <- length(x)
  ij<-do.call(rbind, sapply(1:n, function(i) if(length(x[[i]])) cbind(i, x[[i]])   ))
  if(is.null(ij)) spam(0, ncol=n, nrow=n)
  else spam(list(indices=ij, rep(1, nrow(ij))), ncol=n, nrow=n)
}