#' @title Calculate the netted matrix
#'
#' @description Subtract the transpose and keep positive elements of a weighted
#' adjacency matrix.
#'
#' @param x A square numeric matrix or an object of class unsumnet
#' @return A numeric matrix with the positive elements of x-x' (or x$AW - x$AW' fir unsumnet object)
#' @export
#'
nettedMatrix <- function(x) UseMethod("nettedMatrix")

#' @rdname nettedMatrix
#' @S3method nettedMatrix default
#' @export
nettedMatrix.default <- function(x) {
  
  if(!is.matrix(x) | !is.numeric(x)) stop("Must be a numeric matrix")
  if(nrow(x) != ncol(x)) stop("Must be square matrix")
    
  wOut <- x - t(x)
  
  return(wOut * (wOut>0))
}

#' @rdname nettedMatrix
#' @S3method nettedMatrix unsumnet
#' @export
nettedMatrix.unsumnet <- function(x) {
  nettedMatrix.default(x$AW)
}