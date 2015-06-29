#' @title Calculate the netted matrix
#'
#' @description Subtract the transpose and keep positive elements of a weighted
#' adjacency matrix.
#'
#' @param w A square numeric matrix
#'
#' @return
#' @export
#'
nettedMatrix <- function(w) {
  
  if(!is.matrix(w) | !is.numeric(w)) stop("Must be a numeric matrix")
  if(nrow(w) != ncol(w)) stop("Must be square matrix")
    
  wOut <- w - t(w)
  
  return(wOut * (wOut>0))
}