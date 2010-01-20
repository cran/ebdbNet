\name{visualize}
\alias{visualize}
\title{Visualize Network Structure}
\description{
This function allows the user to visualize the network structure determined
by the EBDBN, either by appropriately formatting results to be imported in
Cytoscape, or invoking the Rgraphviz package to visualize network directly
in R.
}
\usage{
visualize(zscores, sig, format = "Cytoscape", type)
}
\arguments{
  \item{zscores}{The matrix of z scores obtained from \code{\link{ebdbn}}.}
  \item{sig}{The 0-1 matrix indicating absence or presence of an edge obtained from \code{\link{zCutoff}}.}
  \item{format}{Format of visualization}
  \item{type}{Type of network}
}
\details{
Allows visualization of the network structure, as determined by the z-scores
from the posterior distribution of the network (calculated in \code{\link{ebdbn}}),
and the 0-1 matrix identifying the absence or presence of an edge 
(calculated in \code{\link{zCutoff}}).

The argument 'format' takes the value of "Rgraphviz" or 
"Cytoscape", and is used to determine how results should be formatted.
For the former, the package Rgraphviz is invoked to draw a directed
graph. For the latter, a table is created where the first column contains the
identification number of "from" nodes, the second column contains the
identification number of "to" nodes, and the third column contains
the interaction type (1 = activation, -1 = inhibition). This table can
be input into Cytoscape to visualize the network. 

The argument 'type' takes the value of feedback" or "input"
depending on the type of network inferred. 
}
\value{
  \item{edges }{If 'type' = "Cytoscape", a matrix of appropriately formatted values}
}
\author{Andrea Rau}
\examples{
library(ebdbNet)
tmp <- runif(1) ## Initialize random number generator
set.seed(125214) ## Save seed

## Simulate data
simData <- simFunc(R = 5, T = 10, P = 10, v = rep(10, 10), perc = 0.10)
Dtrue <- simData$Dtrue
y <- simData$y
 
####################################################
## Run EB-DBN without hidden states
####################################################
net <- ebdbn(input = "feedback", y, K = 0, 
	conv.1 = 0.15, conv.2 = 0.10, conv.3 = 0.10)

## Calculate sensitivities, specificities, and precisions of D matrix
z <- zCutoff(net$DPost, net$DvarPost)

## Create results matrix for Cytoscape, based on 95 perc significance
cytoscape.results <- visualize(z$z, z$z95, format = "Cytoscape", type = "feedback")

## Visualize using Rgraphviz, if appropriately installed
## library(Rgraphviz)
## visualize(z$z, z$z95, format = "Rgraphviz", type = "feedback")
}
\keyword{methods}
