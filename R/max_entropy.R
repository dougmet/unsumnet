
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
#' @param rs Either a data.frame object containing row sums (first) and column sums
#' (second) and optionally a vector of node names in any column. Or a vector of row
#' sums only if the column sums, \code{cs}, are also supplied. For a data.frame 
#' the row and column sums are extracted and passed to max_entropy.numeric
#' @param cs NumericVector the column sums of the matrix.
#' @param minError Numeric. The algorithm will keep iterating until the mean squared 
#' error against the constraints drops below this value.
#' @param ... extra arguments passed to \code{\link{max_entropy.numeric}}
#'
#' @return A matrix that satisfies the row and column sum constraints or FALSE
#'  if the algorithm failed to converge. Dimension names will be pulled through
#'  if available from the \code{data} or from the names of \code{rs}.
#' @export
#'
max_entropy <- function(rs, ...) UseMethod("max_entropy")

#' @rdname max_entropy
#' @export
#' @examples
#' max_entropy(neast)
max_entropy.data.frame <- function(rs, ...) {
  # Clean the input
  constraints <- process_input(rs)
  
  aw <- max_entropy(constraints[,1], constraints[,2], ...)
  
  if(!is.null(dimnames(constraints))) {
    dimnames(aw) <- list(dimnames(constraints)[[1]], dimnames(constraints)[[1]])
  }
  
  return(aw)
}

#' @rdname max_entropy
#' @export
#' @examples
#' max_entropy(neast$outSum, neast$inSum)
max_entropy.numeric <- function(rs, cs, minError=1e-18, ...) {

  if(length(rs)!=length(cs)) stop("rs and cs must be same length")
  
  if (any(is.na(rs)) | any(rs<0, na.rm=TRUE) | 
      any(is.na(cs)) | any(cs<0, na.rm=TRUE))
    warning("max_entropy is not well defined with missing values")
  else {
    if(sum(rs) != sum(cs)) {
      stop(paste("Sum of row and column sums should be the same when there are no missings. sum(rs)-sum(cs) =",
                 sum(rs)-sum(cs)))
    }
  }
  
  # Create matrix of ones and zero the diagonal
  aw <- matrix(1, nrow=length(rs), ncol=length(cs)) - diag(1, nrow=length(rs))
  
  # Call the C++ function to do the iteration
  aw <- max_entropyCpp(aw, rs, cs, minError);
  
  if(any(is.nan(aw))) {
    warning("max_entropy created NaNs")
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