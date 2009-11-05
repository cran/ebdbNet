\name{zcutoff}
\alias{zcutoff}
\title{Calculate Z-Scores Based on Posterior Means and Variances}
\description{
Function to generate z-scores based on posterior means and variances of network parameters.
}
\usage{
zcutoff(Post, varPost, R, cutoff = 3.89, varlist = TRUE)
}
\arguments{
  \item{Post}{Posterior mean of interaction matrix}
  \item{varPost}{Posterior variance of interaction matrix}
  \item{R}{Number of replicates}
  \item{cutoff}{Optional additional cutoff value for z-score significance}
  \item{varlist}{Indicator for format of varPost}
}
\details{
Calculate the z-scores of a matrix based on a posterior mean and variance from \code{\link{ebdbn}}.  Note that 
the posterior distributions from the linear feedback state space model are all Gaussian.  If the user would
like a different level of significance than 95\%, 99\%, or 99.9\%, this can be specified in the argument
\code{cutoff}.  If the specified posterior variance is in list format (i.e., one covariance matrix for each
row of the posterior mean matrix), the argument \code{varlist} should be equal to TRUE.  Otherwise, if the
posterior variance is in matrix format (of the same dimensions as the posterior mean), the argument
\code{varlist} should be equal to FALSE.
}
\value{
  \item{z }{Matrix of z-scores}
  \item{z95 }{Matrix indicating edges significant at a 95\% level of significance}
  \item{z99 }{Matrix indicating edges significant at a 99\% level of significance}
  \item{z99.9 }{Matrix indicating edges significant at a 99.9\% level of significance}
  \item{zchoice }{Matrix indicating edges significant at a user supplied level of significance}
}
\author{Andrea Rau}
\seealso{\code{\link{ebdbn}}}
\examples{

library(ebdbNet)
tmp <- runif(1) 	## Initialize random number generator
set.seed(125214)	## Save seed

R <- 5 		## 5 replicates
T <- 10 		## 10 time points
P <- 10		## 10 observed variables
v <- rep(10, P)	## Observed variable precisions

## Create artificial posterior mean matrix
DPost <- matrix(rnorm(P*P), nrow = P, ncol = P)

## Create artificial posterior variance matrix
DvarPost <- vector("list", P)
for(i in 1:P) {
	DvarPost[[i]] <- diag(0.5, P)
}

# Use zscore significance level of 95%
z <- zcutoff(DPost, DvarPost, R)$z95 ## 17 edges with z-scores significant at 95%
}
\keyword{methods}

