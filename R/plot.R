
#' Plot method for an \code{unsumnet} network
#'
#' @description Passes the weighted adjacency matrix to \code{\link{plotUnsum}}
#' which is a wrapper around the igraph plot. You will need to install
#' \code{igraph} to use this plot method.
#' @param x an \code{unsumnet} object
#' @param ... extra options passed through to \code{\link{plot.igraph}}
#'
#' @export
#'
plot.unsumnet <- function(x, ...) {
  # Plot the weighted adjacency matrix
  plotUnsum(x$AW, ...)
}

#' Plot function for a weighted adjacency matrix
#'
#' @description This will create a directed \code{igraph} network and call its
#' plotting methods
#' @param aw a matrix
#' @param ... arguments passed through to \code{\link{plot.igraph}}
#'
#' @export
#'
#' @examples
#' plotUnsum(rbind(c(0,1,2),c(0,0,1), c(1,0,0)))
plotUnsum <- function(aw, ...) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  op<-par(mar=c(0,0,1.5,0)) # Change margin settings on R side
  g <- igraph::graph.adjacency(aw, weighted = TRUE)
  igraph::V(g)$size <- 50
  igraph::V(g)$color <- 'gold'
  if(!is.null(dimnames(aw)))
    igraph::V(g)$label <- dimnames(aw)[[1]]
  lo <- igraph::layout.circle(g)
  sw <- mean(aw)/2
  igraph::plot.igraph(g, edge.width=igraph::E(g)$weight/sw,
                      edge.curved=igraph::is.mutual(g),
                      layout=lo, 
                      vertex.frame.color='slateblue',
                      edge.arrow.size=0.7,
                      edge.color=rgb(0,0,0,0.5),
                      ...)
  par(op) # put the settings back
}