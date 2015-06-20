#' Generate unsummed adjacency matrix
#' 
#' @param constraints A matrix or data frame containing the row and column sums.
#' NA or negative values are interpretted as unconstrained.
#' @param nEdges Integer The desired number of edges (non-zero elements) the 
#' adjacency matrix should have.
#' @param maxEdges Logical if set to TRUE then all possible edges will be used.
#' The number of possible edges is given by \code{\link{calcMaxEdges}}. Default
#' is FALSE.
#' @param noReturn If set to TRUE then only one directed edge can go between two
#' nodes. This would be appropriate if only modelling netted positions for example.
#' Default is FALSE.
#' @param mctSchedule Sets the number of MC loops between measurements. Also effects
#' the quench speed as beta is only changed every mctSchedule loops.
#' @param hotTime How many MC loops to run at the high temperature before starting
#' the quench. The bigger the number the more time the network has to truly
#' randomise before quenching.
#' @param beta0 Starting inverse temperature. The smaller the number the hotter,
#' or more random, the starting configuration is.
#' @param betaMax The maximum beta value.
#' @param mu0 The field that couples to the number of edges. The bigger this number
#' the more strict the edge restriction is during the quench.
#' @param coolingRate A number greater than, but close to, 1 that multiplies by \code{beta}
#' every \code{mctSchedule} loops.
#' @param maxTime The maximum number of loops before giving up on that quench and
#' starting again.
#' @param minError When the error drops below this threshold it is decided that it
#' is close enough to a solution. The remaining error is resolved using the row/col
#' iterator. A bigger number will give better performance but risks skewing the
#' distirbution of networks.
#'
#' @description Given the row and column sums for a positive weighted
#' adjacency matrix, this function will generate candidate matrices that are
#' sampled uniformly from the space of possible matrices. This is done using
#' simulated annealing.
#' 
#' @return a matrix who's row and column sums match constraints
#' @author Douglas Ashton
#' @export
unsum <- function(constraints,
                  nEdges=NULL,
                  maxEdges=FALSE,
                  noReturn=FALSE,
                  mctSchedule=100,
                  hotTime=2000,
                  beta0=1e-2,
                  betaMax=1e5,
                  mu0=0.1,
                  coolingRate=1.02,
                  maxTime=1e6,
                  minError=1e-3) {
  
  # Clean the input
  constraints <- processInput(constraints)
  
  # Quick length/type check of other inputs
  stopifnot(is.logical(maxEdges), length(maxEdges)==1)
  stopifnot(is.logical(noReturn), length(noReturn)==1)
  stopifnot(is.numeric(mctSchedule), length(mctSchedule)==1)
  stopifnot(is.numeric(hotTime), length(hotTime)==1)
  stopifnot(is.numeric(beta0), length(beta0)==1)
  stopifnot(is.numeric(betaMax), length(betaMax)==1)
  stopifnot(is.numeric(mu0), length(mu0)==1)
  stopifnot(is.numeric(coolingRate), length(coolingRate)==1)
  stopifnot(is.numeric(maxTime), length(maxTime)==1)
  stopifnot(is.numeric(minError), length(minError)==1)
  
  # I think something bad happens. Stop it.
  if(maxEdges & noReturn)
    stop("Can't have maxEdges=TRUE and noReturn=TRUE together")
  
  if(!maxEdges & is.null(nEdges))
    stop("Must supply nEdges or set maxEdges=TRUE")
  
  if(!is.null(nEdges)) {
    stopifnot(is.numeric(nEdges), length(nEdges)==1)
    if(nEdges > nrow(constraints)*(nrow(constraints)-1))
      stop("Too many edges specified for matrix size")
    
    if(nEdges > calcMaxEdges(constraints[,1], constraints[,2])) {
      warning("More edges specified than possible for these row/col sums")
    }
    
    if(maxEdges)
      stop("Can't specify maxEdges and nEdges together")
  }
  
  # Specify the maximum number of edges if asked for
  if(maxEdges) nEdges <- calcMaxEdges(constraints[,1], constraints[,2])
  
  if(coolingRate<=1) stop("coolingRate must be greater than 1")
  
  unsumcpp(constraints,
           nEdges,
           maxEdges,
           noReturn,
           mctSchedule,
           hotTime,
           beta0,
           betaMax,
           mu0,
           coolingRate,
           maxTime,
           minError)
  
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
#' @export
calcMaxEdges <- function(rs, cs) {
  n1 <- sum(rs>0)
  n2 <- sum(cs>0)
  nb <- sum(rs>0 & cs>0)
  n1*n2 - nb
}