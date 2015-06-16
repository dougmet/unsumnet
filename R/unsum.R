#' Generate unsummed matrix
#' 
#' @param constraints A matrix or data frame containing the row and column sums
#'
#' @description Given the row and column sums for a square matrix as constraints,
#' this function will use simulated generate candidate matrices that are sampled
#' uniformly from the space of possible matrices. This is done using simulated
#' annealing.
#' 
#' @return a matrix who's row and column sums match constraints
#' @author Douglas Ashton

unsum <- function(constraints) {
  
}