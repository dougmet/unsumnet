#' Summarising Unsummed Networks
#'
#' @description \code{\link{summary}} method for class "unsumnet"
#'
#' @param object object of class "unsumnet", usually from a call to
#' \code{\link{unsum}}.
#' @param ... passed on to summary for numeric vectors
#'
#' @export
#'
summary.unsumnet <- function(object, ...) {
  
  nNodes <- nrow(object$A)
  cat("Reconstructed matrix for", nNodes, "nodes.\n\n")
  cat("Target number of edges: ", object$targetEdges, "\n")
  cat("Actual number of edges: ", object$nEdges, "\n\n")
  
  edgeWeights <- object$AW[object$A>0]
  cat("Edge weight summary\n")
  summary(edgeWeights, ...)
  
  od <- options(digits=3)
  on.exit(options(od))
  
  if (ncol(object$AW)<=20) {
    print(object$AW)
  } else {
    cat("\nFirst 20 rows/cols")
    print(object$AW[1:20,1:20])
  }
  
}

#' Print method for \code{unsumnet}
#'
#' @param x an object, usually generated from \code{\link{unsum}}.
#' @param ... passed on to other print methods.
#'
#' @export
print.unsumnet <- function(x, ...) {
  nNodes <- nrow(x$A)
  cat("Reconstructed matrix for", nNodes, "nodes.\n\n")
  cat("Number of edges: ", x$targetEdges, "\n")

  
  od <- options(digits=3)
  on.exit(options(od))
  
  if (ncol(x$AW)<=20) {
    print(x$AW, ...)
  } else {
    cat("\nFirst 20 rows/cols")
    print(x$AW[1:20,1:20], ...)
  }
  
}