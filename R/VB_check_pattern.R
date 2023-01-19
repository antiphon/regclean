#' Check input pattern
#' 
#' @param x Point pattern candidate.
#' 
#' @details 
#' 
#' We use the simple format of a list with
#' $x = Matrix of size n x d, where d is dimension.
#' $bbox = 2 x d matrix, with ranges of bounding box in columns.
#' 
#' @export

VB_check_pattern <- function(x){
  # if ppp:
  if("ppp"%in% is(x)){
    x <- list(x=cbind(x$x, x$y), bbox=cbind(x$window$xrange, x$window$yrange))
  }
  # if list or matrix:
  x <- if(is.list(x)) x else if(!is.null(ncol(x))) list(x=as.matrix(x)) else stop("Can't handle the format of x.")
  x$bbox <- if(!is.null(x$bbox)) x$bbox else apply(x$x, 2, range)
  #'
  x
}