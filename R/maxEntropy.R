
#' Find the so-called maximum entropy solution 
#'
#' @param rs NumericVector the row sums of the matrix
#' @param cs NumericVector the column sums of the matrix
#'
#' @return
#' @export
#'
maxEntropy <- function(rs, cs) {
  
  if (any(is.na(rs))) stop("maxEntropy can't yet cope with missing values in rs")
  if (any(is.na(cs))) stop("maxEntropy can't yet cope with missing values in cs")
  
  if(length(rs)!=length(cs)) stop("rs and cs must be same length")
  
  if(sum(rs) != sum(cs)) stop("Sum of row and column sums should be the same")
  
  aw <- matrix(1, nrow=length(rs), ncol=length(cs)) - diag(1, nrow=length(rs))
  
  #err <- sum((csaw - cs)^2 + (rsaw - rs)^2)
  
  
}