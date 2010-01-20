\name{zCutoff}
\alias{zCutoff}
\title{Calculate Z-Scores Based on Posterior Means and Variances}
\description{
Function to generate z-scores based on posterior means and variances of network parameters.
}
\usage{
zCutoff(Post, varPost)
}
\arguments{
  \item{Post}{Posterior mean of interaction matrix}
  \item{varPost}{Posterior variance of interaction matrix}
}
\details{
Calculate the z-scores of a matrix based on a posterior mean and variance from \code{\link{ebdbn}}. Note that
the posterior distributions from the linear feedback state space model are all Gaussian.  
}
\value{
  \item{z }{Matrix of z-scores}
  \item{z95 }{Matrix indicating edges significant at a 95\% level of significance}
  \item{z99 }{Matrix indicating edges significant at a 99\% level of significance}
  \item{z99.9 }{Matrix indicating edges significant at a 99.9\% level of significance}
}
\author{Andrea Rau}
\seealso{\code{\link{ebdbn}}}
\examples{
library(ebdbNet)
tmp <- runif(1) ## Initialize random number generator
set.seed(125214) ## Save seed

## 10 observed variables
P <- 10 

## Create artificial posterior mean and covariance matrix
DPost <- matrix(rnorm(P*P), nrow = P, ncol = P)
DvarPost <- vector("list", P)
for(i in 1:P) {
	DvarPost[[i]] <- diag(0.5, P)
}

# Use zscore significance level of 95%
z <- zCutoff(DPost, DvarPost)$z95 ## 17 edges with z-scores significant at 95%
}
\keyword{methods}
