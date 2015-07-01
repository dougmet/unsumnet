#' Summarising Unsummed Networks
#'
#' @description \code{\link{summary}} method for class "unsumnet"
#'
#' @param x an object of class "unsumnet", usually from a call to
#' \code{\link{unsum}}.
#'
#' @export
#'
summary.unsumnet <- function(x) {
  
  nNodes <- nrow(x$A)
  cat("Reconstructed matrix for", nNodes, "nodes.\n\n")
  cat("Target number of edges: ", x$targetEdges, "\n")
  cat("Actual number of edges: ", x$nEdges, "\n\n")
  
  edgeWeights <- x$AW[x$A>0]
  cat("Edge weight summary\n")
  summary(edgeWeights)
  
  if (ncol(x$AW)<=10) {
    print(x$AW)
  } else {
    cat("\nTop 10 rows/cols")
    print(x$AW[1:10,1:10])
  }
  
}

#' Print method for \code{unsumnet}
#'
#' @param x and object, usually generated from \code{\link{unsum}}.
#'
#' @export
print.unsumnet <- function(x) {
  nNodes <- nrow(x$A)
  cat("Reconstructed matrix for", nNodes, "nodes.\n\n")
  cat("Number of edges: ", x$targetEdges, "\n")

  if (ncol(x$AW)<=10) {
    print(x$AW)
  } else {
    cat("\nTop 10 rows/cols")
    print(x$AW[1:10,1:10])
  }
}