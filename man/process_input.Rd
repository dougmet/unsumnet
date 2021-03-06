% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_input.R
\name{process_input}
\alias{process_input}
\title{Process the input constraints}
\usage{
process_input(constraints, useRowNames = FALSE)
}
\arguments{
\item{constraints}{A matrix or data frame containing the row and column sums}

\item{useRowNames}{Logical if TRUE it the vertex names will be taken from 
the data frame's row names. This is not used for matrices.}
}
\value{
an object ready for sending to the C++ functions as input
}
\description{
This function checks and prepares the network constraints before
calling the C++ code.
}
\author{
Douglas Ashton
}

