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
  
  if(is.matrix(w)) {
    
    if(nrow(w) != ncol(w)) stop("Must be square matrix")
    
    wOut <- w - t(w)
  } else {
    stop("Must be numeric matrix")
  }
  
  return(wOut * (wOut>0))
}