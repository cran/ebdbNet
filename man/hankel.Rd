\name{hankel}
\alias{hankel}
\title{Perform Singular Value Decomposition of Block-Hankel Matrix}
\description{
This function constructs a block-Hankel matrix based on time-course data, performs the
subsequent singular value decomposition (SVD) on this matrix, and returns the number
of large singular values as defined by a user-supplied cutoff criterion.
}
\usage{
hankel(y, R, T, P, lag, cutoff, type)
}
\arguments{
  \item{y}{A list of R (PxT) matrices of observed time course profiles}
  \item{R}{Number of replicates}
  \item{T}{Number of time points}
  \item{P}{Number of observed variables (e.g., genes)}
  \item{lag}{Maximum relevant time lag to be used in constructing the block-Hankel matrix}
  \item{cutoff}{Cutoff to be used, determined by desired percent of total variance explained}
  \item{type}{Method to combine results across replicates ("median" or "mean")}
}
\details{
Constructs the block-Hankel matrix \eqn{H} of autocovariances of time series observations is constructed (see
references for additional information), where the maximum relevant time lag must be specified as \code{lag}.
In the context of gene networks, this corresponds to the maximum relevant biological time lag between a gene
and its regulators.  This quantity is experiment-specific, but will generally be small for gene expression
studies (on the order of 1, 2, or 3).

The singular value decomposition of \eqn{H} is performed, and the singular values are ordered by size and scaled
by the largest singular value.  Note that if there are T time points in the data, only the first (T - 1)
singular values will be non-zero.

To choose the number of large singular values, we wish to find the point at which the inclusion of an additional
singular value does not increase the amount of explained variation enough to justify its inclusion (similar to
choosing the number of components in a Principal Components Analysis). The user-supplied value of \code{cutoff}
gives the desired percent of variance explained by the first set of K principal components.  The algorithm
returns the value of K, which may subsequently be used as the dimension of the hidden state in \code{\link{ebdbn}}.

The argument \code{type} takes the value of "median" or "mean", and is used to determine how results from
replicated experiments are combined (i.e., median or mean of the
per-replicate final hidden state dimension).
}
\value{
  \item{svs }{Vector of singular values of the block-Hankel matrix \eqn{H}}
  \item{dim }{Number of large singular values, as determined by the user-supplied cutoff}
}
\references{
Masanao Aoki and Arthur Havenner (1991). State space modeling of multiple time series. 
\emph{Econometric Reviews} 10(1), 1-59.

Martina Bremer (2006). \emph{Identifying regulated genes through the correlation structure of time dependent 
microarray data}. Ph. D. thesis, Purdue University.

Andrea Rau, Florence Jaffrezic, Jean-Louis Foulley, and R. W. Doerge (2009). An Empirical Bayesian
Method for Estimating Biological Networks from Temporal Microarray Data. 
\emph{Statistical Applications in Genetics and Molecular Biology} (submitted).
}
\author{Andrea Rau}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

library(ebdbNet)
tmp <- runif(1) 	## Initialize random number generator
set.seed(125214)	## Save seed

R <- 5 		## 5 replicates
T <- 10 		## 10 time points
P <- 10		## 10 observed variables
v <- rep(10, P)	## Observed variable precisions

## Create D matrix with 10% of possible edges non-null
## Non-null edges are either U(0.2, 1) or U(-1, -0.2) with equal probability
D <- matrix(0, nrow = P, ncol = P)        
index <- expand.grid(seq(1:P),seq(1:P))
selected.index <- sample(seq(1:(P*P)), ceiling(0.10 * P * P))
selected.edges <- index[selected.index,]   
for(edge in 1:ceiling(0.10 * P * P)) {
	tmp <- runif(1)
	if(tmp > 0.5) {
    		D[selected.edges[edge,1], selected.edges[edge,2]] <- 
			runif(1, min = 0.2, max = 1)
	}
	else {
    		D[selected.edges[edge,1], selected.edges[edge,2]] <- 
			runif(1, min = -1, max = -0.2)
	}	
}
Dtrue <- abs(sign(D))
     
## Simulate data
y <- vector("list", R)
for(r in 1:R) {
	y[[r]] <- matrix(NA, nrow = P, ncol = T)
      y[[r]][,1] <- rnorm(P, mean = 0, sd = sqrt(v^(-1)))
      for(t in 2:T) {
       	y[[r]][,t] <- D \%*\% y[[r]][, t - 1] + 
			as.matrix(rnorm(P,mean = 0, sd = sqrt(v^(-1))))
      }
}

## Determine the number of hidden states to be estimated (with lag = 1)
lag1 <- hankel(y, R, T, P, lag = 1, cutoff = 0.90, type = "median")$dim
K <- lag1  ## K = 5

}
\keyword{methods}
