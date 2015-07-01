#' @title Calculate the netted matrix
#'
#' @description Subtract the transpose and set negative elements to zero for a weighted
#' adjacency matrix.
#'
#' @param x A square numeric matrix or an object of class \code{unsumnet}
#' @return A numeric matrix with the positive elements of x-x' (or x$AW - x$AW'
#' for \code{unsumnet} object) and zero for negative elements.
#' @export
#'
nettedMatrix <- function(x) UseMethod("nettedMatrix")

#' @rdname nettedMatrix
#' @export
#' @examples
#' nettedMatrix(matrix(rnorm(25), nrow=5))
nettedMatrix.default <- function(x) {
  
  if(!is.matrix(x) | !is.numeric(x)) stop("Must be a numeric matrix")
  if(nrow(x) != ncol(x)) stop("Must be square matrix")
    
  wOut <- x - t(x)
  
  # set negative elements to zero and return
  return(wOut * (wOut>0))
}

#' @rdname nettedMatrix
#' @export
nettedMatrix.unsumnet <- function(x) {
  nettedMatrix.default(x$AW)
}