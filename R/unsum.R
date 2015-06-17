#' Generate unsummed adjacency matrix
#' 
#' @param constraints A matrix or data frame containing the row and column sums.
#' NA or negative values are interpretted as unconstrained.
#' @description Given the row and column sums for a positive weighted
#' adjacency matrix, this function will generate candidate matrices that are
#' sampled uniformly from the space of possible matrices. This is done using
#' simulated annealing.
#' 
#' @return a matrix who's row and column sums match constraints
#' @author Douglas Ashton
#' @export
unsum <- function(constraints,
                  nedges=NULL,
                  maxEdges=FALSE,
                  noReturn=FALSE,
                  mct_schedule=100,
                  hot_time=2000,
                  beta0=1e-2,
                  betamax=1e5,
                  mu0=0.1,
                  cooling_rate=1.02,
                  max_time=1e6,
                  cgmax=1e-3) {
  
  # Clean the input
  constraints <- processInput(constraints)
  
  # Quick length/type check of other inputs
  stopifnot(is.logical(maxEdges), length(maxEdges)==1)
  stopifnot(is.logical(noReturn), length(noReturn)==1)
  stopifnot(is.numeric(mct_schedule), length(mct_schedule)==1)
  stopifnot(is.numeric(hot_time), length(hot_time)==1)
  stopifnot(is.numeric(beta0), length(beta0)==1)
  stopifnot(is.numeric(betamax), length(betamax)==1)
  stopifnot(is.numeric(mu0), length(mu0)==1)
  stopifnot(is.numeric(cooling_rate), length(cooling_rate)==1)
  stopifnot(is.numeric(max_time), length(max_time)==1)
  stopifnot(is.numeric(cgmax), length(cgmax)==1)
  
  # I think something bad happens. Stop it.
  if(maxEdges & noReturn)
    stop("Can't have maxEdges=TRUE and noReturn=TRUE together")
  
  if(!is.null(nedges)) {
    stopifnot(is.numeric(nedges), length(nedges)==1)
    if(nedges > nrow(constraints)*(nrow(constraints)-1))
      stop("Too many edges specified for matrix size")
    
    if(nedges > calcMaxEdges(constraints[,1], constraints[,2])) {
      warning("More edges specified than possible for these row/col sums")
    }
    
    if(maxEdges)
      stop("Can't specify maxEdges and nedges together")
  }
  
  # Specify the maximum number of edges if asked for
  if(maxEdges) nedges <- calcMaxEdges(constraints[,1], constraints[,2])
  
  unsumcpp(constraints,
           nedges,
           maxEdges,
           noReturn,
           mct_schedule,
           hot_time,
           beta0,
           betamax,
           mu0,
           cooling_rate,
           max_time,
           cgmax)
  
}

#' Calculate maximum possible edges
#' 
#' @param rs numeric vector, row sums
#' @param cs numeric vector, col sums
#' @description For a set of row and column sums compute maximum number of
#' edges assuming no self-loops
#' 
#' @return maximum possible number of edges
#' @author Douglas Ashton
calcMaxEdges <- function(rs, cs) {
  n1 <- sum(rs>0)
  n2 <- sum(cs>0)
  nb <- sum(rs>0 & cs>0)
  n1*n2 - nb
}