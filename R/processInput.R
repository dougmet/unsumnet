#' Process the input constraints
#' 
#' @param constraints A matrix or data frame containing the row and column sums
#' @param useRowNames Logical if TRUE it the vertex names will be taken from 
#' the data frame's row names. This is not used for matrices.
#' @description This function checks and prepares the network constraints before
#' calling the C++ code.
#' 
#' @return an object ready for sending to the C++ functions as input
#' @author Douglas Ashton
#' 
processInput <- function(constraints, useRowNames=FALSE) {
  
  # Lots and lots of checking
  
  # So you only need to change it in one place
  defaultNames <- function() paste0("V",1:nrow(constraints))
  
  if (is.data.frame(constraints)) {
    # For data frames the idea is to get the names out and ensure it's numeric
    
    # Try to get the names, first row names
    charCols <- which(sapply(constraints, class) %in% c("character","factor"))
    
    if (useRowNames) {
      vertexNames <- row.names(constraints)
    } else if(length(charCols) > 0) {
      # look for character or factor columns
      vertexNames <- as.character(constraints[,charCols[1]])
    } else {
      vertexNames <- defaultNames()
    }
    
    # Warning or stop? Can't have duplicates.
    if (any(duplicated(vertexNames))) {
      warning("Duplicates in vertex names: ",
              vertexNames[duplicated(vertexNames)], ". Using defaults.")
      vertexNames <- defaultNames()
    }
    
    # Minimum row and column constraints
    numCols <- which(sapply(constraints, is.numeric))
    if (length(numCols)<2) {
      stop("Need at least two numeric columns in constraints data frame")
    }
    
    # Now what we really want, a matrix
    constraints <- as.matrix(constraints[,numCols])
    dimnames(constraints)[[1]] <- vertexNames
    
  } else if (is.matrix(constraints)) {
    
    if(!is.numeric(constraints)) {
      stop("Input matrix must be numeric")
    }
    
    # Minimum row and column constraints
    if (ncol(constraints<2)) {
      stop("Need at least two columns in constraints matrix")
    }
    
    # Vertex names from row names
    if(is.null(dimnames(constraints)[[1]])) {
      vertexNames <- defaultNames()
    } else {
      vertexNames <- dimnames(constraints)[[1]]
    }
    
  } else {
    stop("constraints must be a matrix or a data frame")
  }
  
  # On the C++ side negative values signify unconstrained. Replace the NAs.
  constraints <- replace(constraints, is.na(constraints), -1)
  
  return(constraints)
  
}
