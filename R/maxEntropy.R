
#' Find the so-called maximum entropy solution 
#'
#' @description The maximum entropy solution is obtained by spreading the edge
#' weights as evenly as possible throught the weighted adjacency matrix. This
#' version forces the diagonal elements to be zero along with any rows (columns)
#' that have a zero row (column) sum in the input.
#' 
#' Missing values are allowed but will result in a warning. The usefulness of 
#' such results is up to the user.
#' 
#' @param data a data.frame object containing row sums (first) and column sums
#' (second) and optionally a vector of node names in any column. The row and 
#' column sums are extracted and passed to maxEntropy.numeric
#' @param rs NumericVector the row sums of the matrix.
#' @param cs NumericVector the column sums of the matrix.
#'
#' @return A matrix that satisfies the row and column sum constraints or FALSE
#'  if the algorithm failed to converge. Dimension names will be pulled through
#'  if available from the \code{data} or from the names of \code{rs}.
#' @export
#'
maxEntropy <- function(data, ...) UseMethod("maxEntropy")

#' @rdname maxEntropy
#' @export
#' @examples
#' maxEntropy(neast)
maxEntropy.data.frame <- function(data, ...) {
  # Clean the input
  constraints <- processInput(data)
  
  aw <- maxEntropy(constraints[,1], constraints[,2], ...)
  
  if(!is.null(dimnames(constraints))) {
    dimnames(aw) <- list(dimnames(constraints)[[1]], dimnames(constraints)[[1]])
  }
  
  return(aw)
}

#' @rdname maxEntropy
#' @export
#' @examples
#' maxEntropy(neast$outSums, neast$inSums)
maxEntropy.numeric <- function(rs, cs, minError=1e-18) {

  if(length(rs)!=length(cs)) stop("rs and cs must be same length")
  
  if (any(is.na(rs)) | any(rs<0, na.rm=TRUE) | 
      any(is.na(cs)) | any(cs<0, na.rm=TRUE))
    warning("maxEntropy is not well defined with missing values")
  else {
    if(sum(rs) != sum(cs)) {
      stop(paste("Sum of row and column sums should be the same when there are no missings. sum(rs)-sum(cs) =",
                 sum(rs)-sum(cs)))
    }
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
  err <- sum((csaw - cs)^2 + (rsaw - rs)^2, na.rm=TRUE)
  
  # If the error is too big then it must have failed.
  if (err>minError) {
    return(FALSE)
  }
  
  # Try to pull through dimnames
  if(!is.null(names(rs))) {
    dimnames(aw) <- list(names(rs), names(rs))
  }
  
  return(aw)
}