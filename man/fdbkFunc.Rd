\name{fdbkFunc}
\alias{fdbkFunc}
\title{Feedback Function}
\description{
Helper function for \code{\link{ebdbn}} when feedback is included in a network.
}
\usage{
fdbkFunc(y)
}
\arguments{
  \item{y}{A list of R (PxT) matrices of observed time course profiles}
}
\details{
This function puts the observed variables (e.g., genes) in the correct format to be
incorporated as feedback in \code{\link{ebdbn}}.  Loosely speaking, if the matrices in the list of observed data
are of the form \eqn{y_1, ..., y_T}{y(1), ..., y(T)}, the function returns a list
of matrices equal of the form \eqn{0, y_1, ..., y_{T-1}}{0, y(1), ..., y(T-1)} which
are then used as inputs in the state space model.
}
\value{
A list of R (PxT) matrices of feedback profiles for the observed data.
}
\author{Andrea Rau}
\seealso{\code{\link{ebdbn}}}
\examples{
library(ebdbNet)
tmp <- runif(1) ## Initialize random number generator
set.seed(125214) ## Save seed

## Simulate autoregressive data
## 5 replicates, 10 time points, 10 genes
y <- simFunc(R = 5, T = 10, P = 10, v = rep(10,10), perc = 0.15)$y

## Create feedback profiles
u <- fdbkFunc(y)
}
\keyword{methods}
