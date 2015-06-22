
#' Find the so-called maximum entropy solution 
#'
#' @description The maximum entropy solution is obtained by spreading the edge
#' weights as evenly as possible throught the weighted adjacency matrix. This
#' version forces the diagonal elements to be zero along with any rows (columns)
#' that have a zero row (column) sum in the input. Missing values are no
#'
#' @param rs NumericVector the row sums of the matrix
#' @param cs NumericVector the column sums of the matrix
#'
#' @return A matrix that satisfies the row and column sum constraints or FALSE
#'  if the algorithm failed to converge
#' @export
#'
maxEntropy <- function(rs, cs, minError=1e-18) {
  
  if (any(is.na(rs))) stop("maxEntropy can't yet cope with missing values in rs")
  if (any(is.na(cs))) stop("maxEntropy can't yet cope with missing values in cs")
  if (any(rs<0)) stop("maxEntropy can't yet cope with negative values in rs")
  if (any(cs<0)) stop("maxEntropy can't yet cope with negative values in cs")
  
  if(length(rs)!=length(cs)) stop("rs and cs must be same length")
  
  if(sum(rs) != sum(cs)) {
    stop(paste("Sum of row and column sums should be the same. sum(rs)-sum(cs) =",
               sum(rs)-sum(cs)))
  }
  
  # Create matrix of ones and zero the diagonal
  aw <- matrix(1, nrow=length(rs), ncol=length(cs)) - diag(1, nrow=length(rs))
  
  # Call the C++ function to do the iteration
  aw <- maxEntropyCpp(aw, rs, cs, minError);
  
  if(any(is.nan(aw))) {
    warning("maxEntropy created NaNs")
    return(FALSE)
  }
    
  
  # Calculate the final row/col sums
  csaw <- colSums(aw)
  rsaw <- rowSums(aw)
  
  # Compare with input to check the error value
  err <- sum((csaw - cs)^2 + (rsaw - rs)^2)
  
  # If the error is too big then it must have failed.
  if (err>minError)
    return(FALSE)
  
  
  return(aw)
}