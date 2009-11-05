\name{sumfunc}
\alias{sumfunc}
\title{Sum Function}
\description{
Helper function for \code{\link{hankel}} to determine what percentage of variance is explained
by the singular values of the block-Hankel matrix.
}
\usage{
sumfunc(x, cutoff)
}
\arguments{
  \item{x}{Vector of singular values from singular value decomposition of block-Hankel matrix}
  \item{cutoff}{Value to determine cutoff to be considered for singular values (e.g., 0.90)}
}
\details{
This function is used to help determine the number of large singular values, which in turn
determines the hidden state dimension to be estimated.
}
\value{
The number of large singular values, as defined by the user-supplied cutoff.
}
\author{Andrea Rau}
\seealso{\code{\link{hankel}}}
\examples{

library(ebdbNet)
tmp <- runif(1)		## Initialize random number generator
set.seed(1241901)		## Set seed

x <- matrix(rnorm(100), nrow = 10, ncol = 10)
svd <- svd(x)$d
dim <- sumfunc(svd, 0.90) ## 7 singular values contribute 90% of total variance
}
\keyword{methods}
