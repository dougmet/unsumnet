% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/max_entropy.R
\name{max_entropy}
\alias{max_entropy}
\alias{max_entropy.data.frame}
\alias{max_entropy.numeric}
\title{Find the so-called maximum entropy solution}
\usage{
max_entropy(rs, ...)

\method{max_entropy}{data.frame}(rs, ...)

\method{max_entropy}{numeric}(rs, cs, minError = 1e-18, ...)
}
\arguments{
\item{rs}{Either a data.frame object containing row sums (first) and column sums
(second) and optionally a vector of node names in any column. Or a vector of row
sums only if the column sums, \code{cs}, are also supplied. For a data.frame 
the row and column sums are extracted and passed to max_entropy.numeric}

\item{...}{extra arguments passed to \code{\link{max_entropy.numeric}}}

\item{cs}{NumericVector the column sums of the matrix.}

\item{minError}{Numeric. The algorithm will keep iterating until the mean squared 
error against the constraints drops below this value.}
}
\value{
A matrix that satisfies the row and column sum constraints or FALSE
 if the algorithm failed to converge. Dimension names will be pulled through
 if available from the \code{data} or from the names of \code{rs}.
}
\description{
The maximum entropy solution is obtained by spreading the edge
weights as evenly as possible throught the weighted adjacency matrix. This
version forces the diagonal elements to be zero along with any rows (columns)
that have a zero row (column) sum in the input.

Missing values are allowed but will result in a warning. The usefulness of 
such results is up to the user.
}
\examples{
max_entropy(neast)
max_entropy(neast$outSum, neast$inSum)
}

