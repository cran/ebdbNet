\name{ebdbNet-package}
\alias{ebdbNet-package}
\alias{ebdbNet}
\docType{package}
\title{
Empirical Bayes Dynamic Bayesian Network (EBDBN) Inference
}
\description{
This package is used to infer the adjacency matrix of a network from time course data using an
empirical Bayes estimation procedure based on Dynamic Bayesian Networks.
}
\details{
\tabular{ll}{
Package: \tab ebdbNet\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-11-05\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
Posterior distributions (mean and variance) of network parameters are estimated using time-course data based on
a linear feedback state space model that allows for a set of hidden states to be incorporated.
The algorithm is composed of three principal parts: choice of hidden state dimension (see \code{\link{hankel}}),
estimation of hidden states via the Kalman filter and smoother, and calculation of posterior distributions
based on the empirical Bayes estimation of hyperparameters in a hierarchical Bayesian framework
(see \code{\link{ebdbn}}).
}
\author{
Andrea Rau

Maintainer: Andrea Rau <arau@stat.purdue.edu>
}
\references{
Andrea Rau, Florence Jaffrezic, Jean-Louis Foulley, and R. W. Doerge (2009). An Empirical Bayesian
Method for Estimating Biological Networks from Temporal Microarray Data. 
\emph{Statistical Applications in Genetics and Molecular Biology} (submitted).
}
\keyword{methods}
\examples{


library(ebdbNet)
library(GeneNet) 	## Load GeneNet package to use T-cell activation data

data(tcell)		## Load T-cell activation data
## Use full dataset (all 44 replicates)
tc44 <- combine.longitudinal(tcell.10, tcell.34)
## Descriptions of gene names and functions	
desc <- tcell.gene.description				

tmp <- runif(1) 	## Initialize random number generator
set.seed(4568818)	## Set seed
R <- 44		## 44 replicates
T <- 10 		## 10 time points		
P <- 58		## 58 observed variables

## Put data into correct format for algorithm 
## (List, with one matrix per replicate (P rows and T columns)
tcell.dat <- vector("list", R)
for(r in 1:R) {
	tcell.dat[[r]] <- t(tc44[seq(from = r, to = R*T, by = R),])
}

## Use only subset of T-cell data for faster example
R <- 5 		## 5 replicates
T <- 10		## 10 time points
P <- 15		## 15 genes
tcell.sub.dat <- vector("list", R)
for(r in 1:R) {
	tcell.sub.dat[[r]] <- tcell.dat[[r]][1:P, 1:T]
}

####################################################
# Run EBDBN 
####################################################

K <- 0

## Choose alternative value of K using hankel if hidden states to be estimated
# lag1 <- henkel(tcell.sub.dat, R, T, P, lag = 1, cutoff = 0.50, 
#	type = "median")$svs
# K <- 3

## Choose appropriate initial values for delta and v
## Pick any initial values for x, hyperparameters not involved in simplified model
x.0 <- alpha.0 <- beta.0 <- gamma.0 <- mu.0 <- sigma.0 <- 0
delta.0 <- runif(P)
v.0 <- runif(P)

## If K > 0, set appropriate initial values for x and hyperparameters
# x.0 <- vector("list", R)
# for(i in 1:R) {x.0[[i]] <- matrix(rnorm(K*T, mean = 0, sd = 1),nrow = K, ncol = T)}
# alpha.0 <- runif(K)
# beta.0 <- runif(P)
# gamma.0 <- runif(K)
# delta.0 <- runif(P)
# v.0 <- runif(P)
# mu.0 <- rep(0, K)
# sigma.0 <- rep(1, K)

## Run algorithm
net <- ebdbn(tcell.sub.dat, P, K, R, T, x.0, alpha.0, beta.0, gamma.0,
	delta.0, v.0, mu.0, sigma.0, conv.1 = 0.10,
    	conv.2 = 0.01, conv.3 = 0.01, maxiter = 1)

## Use z-score significance level of 99%
## (24 edges identified at significance level of 99%)
z <- zcutoff(net$DPost, net$DvarPost, R)$z99 

}
