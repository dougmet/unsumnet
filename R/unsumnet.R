#' unsumnet: Unsumming summed networks
#'
#' The unsumnet package uses simulated annealing to reconstruct networks
#' from aggregated data, specifically the row and column sums of the
#' weighted adjacency matrix.
#' 
#' 
#' @section Primary functionality:
#' \code{\link{unsum}} sample from the set of networks satisfying the row and 
#' column sum constraints. Extra control over the sparsity (number of edges) is
#' provided.
#' 
#' @section Extra functions:
#' \code{\link{maxEntropy}} calculates the commonly used "maximum entropy" 
#' (outer product) solution by iteratively forcing row and column constraints.
#' \code{\link{nettedMatrix}} provides the netted exposures from the gross matrix.
#'
#' @section Data sets:
#' \code{\link{neast}} a fictional banking system to demonstrate funcionality.
#'
#' @docType package
#' @name unsumnet
#' @importFrom Rcpp evalCpp
#' @useDynLib unsumnet
NULL
