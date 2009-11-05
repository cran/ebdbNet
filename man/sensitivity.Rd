\name{sensitivity}
\alias{sensitivity}
\title{Calculate Sensitivity and Specificity of a Network}
\description{
Function to calculate the true positives (TP), true negatives (TN), false positives (FP), and 
false negatives (FN) of an estimated network, given the structure of the true network.
}
\usage{
sensitivity(trueMatrix, estMatrix, P)
}
\arguments{
  \item{trueMatrix}{Structure of the true network}
  \item{estMatrix}{Structure of the estimated network}
  \item{P}{Number of observed variables}
}
\details{
The matrices \code{trueMatrix} and \code{estMatrix} must both be of dimension (PxP).
Both matrices should have 0's for entries where no gene-gene interaction is present,
and non-zero entries where edges are non-null.
}
\value{
  \item{TP }{Number of true positives}
  \item{FP }{Number of false positives}
  \item{FN }{Number of false negatives}
  \item{TN }{Number of true negatives}
}
\author{Andrea M. Rau}
\seealso{\code{\link{calcAUC}}}
\examples{

library(ebdbNet)
tmp <- runif(1)	## Initialize random number generator
set.seed(16933)	## Set seed
P <- 10		## 10 genes

## Create artificial true D matrix
Dtrue <- matrix(0, nrow = P, ncol = P)   
index <- expand.grid(seq(1:P),seq(1:P))
selected.index <- sample(seq(1:(P*P)), ceiling(0.25 * P * P))
selected.edges <- index[selected.index,]   
for(edge in 1:ceiling(0.25 * P * P)) {
	tmp <- runif(1)
	if(tmp > 0.5) {
    		Dtrue[selected.edges[edge,1], selected.edges[edge,2]] <- 
			runif(1, min = 0.2, max = 1)
	}
	else {
    		Dtrue[selected.edges[edge,1], selected.edges[edge,2]] <- 
			runif(1, min = -1, max = -0.2)
	}	
}

## Create artificial estimated D matrix
Dest <- matrix(0, nrow = P, ncol = P)   
index <- expand.grid(seq(1:P),seq(1:P))
selected.index <- sample(seq(1:(P*P)), ceiling(0.25 * P * P))
selected.edges <- index[selected.index,]   
for(edge in 1:ceiling(0.25 * P * P)) {
	tmp <- runif(1)
	if(tmp > 0.5) {
    		Dest[selected.edges[edge,1], selected.edges[edge,2]] <- 
			runif(1, min = 0.2, max = 1)
	}
	else {
    		Dest[selected.edges[edge,1], selected.edges[edge,2]] <- 
			runif(1, min = -1, max = -0.2)
	}	
}
     
check <- sensitivity(Dtrue, Dest, P)
check$TP	## 5 True Positives
check$FP	## 20 False Positives
check$TN	## 55 True Negatives
check$FN	## 20 False Negatives

}
\keyword{methods}

