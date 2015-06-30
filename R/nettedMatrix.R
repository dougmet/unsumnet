#' @title Calculate the netted matrix
#'
#' @description Subtract the transpose and keep positive elements of a weighted
#' adjacency matrix.
#'
#' @param w A square numeric matrix
#'
#' @return A numeric matrix with the positive elements of w-w'
#' @export
#'
nettedMatrix <- function(w) UseMethod("nettedMatrix")

#' @title Calculate the netted matrix (default)
#'
#' @description Subtract the transpose and keep positive elements of a weighted
#' adjacency matrix.
#'
#' @param w A square numeric matrix
#'
#' @return A numeric matrix with the positive elements of w-w'
#' @export
#'
nettedMatrix.default <- function(w) {
  
  if(!is.matrix(w) | !is.numeric(w)) stop("Must be a numeric matrix")
  if(nrow(w) != ncol(w)) stop("Must be square matrix")
    
  wOut <- w - t(w)
  
  return(wOut * (wOut>0))
}

#' @title Calculate the netted matrix (default)
#'
#' @description Subtract the transpose and keep positive elements of a weighted
#' adjacency matrix.
#'
#' @param x An \code{unsumnet} object usually from \code{\link{unsum}}.
#'
#' @return A numeric matrix with the positive elements of aw-aw'
#' @export
#'
nettedMatrix.unsumnet <- function(x) {
  nettedMatrix.default(x$AW)
}