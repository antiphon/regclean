#' Bounding box volume
#' @export
bbox_volume <- function(bb){
  prod(apply(bb, 2, diff))
}

#' Expand box with a constant
#' @export
bbox_expand <- function(bb, r){
  apply(bb, 2, function(v)v+c(-1,1)*r)  
}
