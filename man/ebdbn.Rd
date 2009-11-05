\name{ebdbn}
\alias{ebdbn}
\title{Empirical Bayes Dynamic Bayesian Network (EBDBN) Estimation}
\description{
A function to infer the posterior mean and variance of network parameters
using an empirical Bayes estimation procedure for a Dynamic Bayesian Network (DBN).
}
\usage{
ebdbn(y, P, K, R, T, x.0, alpha.0, beta.0, gamma.0, delta.0, v.0, 
mu.0, sigma.0, conv.1 = 0.15, conv.2 = 0.05, conv.3 = 0.01, maxiter = 100)
}
\arguments{
  \item{y}{A list of R (PxT) matrices of observed time course profiles}
  \item{P}{Number of observed variables (e.g., genes)}
  \item{K}{Number of hidden states}
  \item{R}{Number of replicates}
  \item{T}{Number of time points}
  \item{x.0}{A list of R (KxT) matrices of initial values for hidden state profiles}
  \item{alpha.0}{(Kx1) vector of initial values for \eqn{\alpha}{alpha}}
  \item{beta.0}{(Px1) vector of initial values for \eqn{\beta}{beta}}
  \item{gamma.0}{(Kx1) vector of initial values for \eqn{\gamma}{gamma}}
  \item{delta.0}{(Px1) vector of initial values for \eqn{\delta}{delta}}
  \item{v.0}{(Px1) vector of initial values for \eqn{V} (e.g., gene precisions)}
  \item{mu.0}{(Kx1) vector of initial values for \eqn{\mu}{mu} (hidden state mean)}
  \item{sigma.0}{(Kx1) vector of initial values for the diagonal of the covariance matrix \eqn{\Sigma}{Sigma} (hidden state covariance)}
  \item{conv.1}{Value of convergence criterion 1}
  \item{conv.2}{Value of convergence criterion 2}
  \item{conv.3}{Value of convergence criterion 3}
  \item{maxiter}{Value of the maximum number of overall iterations to be allowed before the algorithm terminates}
}
\details{
This function infers the parameters of a network, based on the state space model
\deqn{x_t = Ax_{t-1} + By_{t-1} + w_t}{x(t) = Ax(t-1) + By(t-1) + w(t)}
\deqn{y_t = Cx_t + Dy_{t-1} + z_t}{y(t) = Cx(t) + Dy(t-1) + z(t)}
where \eqn{x_t}{x(t)} represents the expression of K hidden states at time \eqn{t}, \eqn{y_t}{y(t)} represents the 
expression of P observed states (e.g., genes) at time \eqn{t}, \eqn{w_t \sim MVN(0,I)}{w(t) ~ MVN(0,I)}, and 
\eqn{z_t \sim MVN(0, V^{-1})}{z(t) ~ MVN(0, V^(-1))}, with \eqn{V = diag(v_1, ..., v_P)}{V = diag(v_1, ..., v_P)}.  
Note that the dimensions of the matrices \eqn{A}, \eqn{B}, \eqn{C}, and \eqn{D} are
(KxK), (KxP), (PxK), and (PxP), 
respectively.  The parameters of greatest interest are contained in the matrix
\eqn{D}, which encodes the direct interactions among observed variables from one time to the next.

The value of K is chosen prior to running the algorithm by using \code{\link{hankel}}.  The hidden states are
estimated using the classic Kalman filter.  Posterior distributions of \eqn{A}, \eqn{B}, \eqn{C}, and \eqn{D} are estimated using an
empirical Bayes procedure based on a hierarchical Bayesian structure defined over the parameter set.  Namely, if
\eqn{a_{(j)}}{a(j)}, \eqn{b_{(j)}}{b(j)}, \eqn{c_{(j)}}{c(j)}, \eqn{d_{(j)}}{d(j)}, denote vectors made up
of the rows of matrices \eqn{A}, \eqn{B}, \eqn{C}, and \eqn{D} respectively, then
\deqn{a_{(j)}\vert \alpha \sim N(0, diag(\alpha)^{-1})}{a(j)|\alpha ~ N(0, diag(alpha)^(-1))} 
\deqn{b_{(j)}\vert \beta \sim N(0, diag(\beta)^{-1})}{b(j)|\beta ~ N(0, diag(beta)^(-1))} 
\deqn{c_{(j)}\vert \gamma \sim N(0, diag(\gamma)^{-1})}{c(j)|\gamma ~ N(0, diag(gamma)^(-1))} 
\deqn{d_{(j)}\vert \delta \sim N(0, diag(\delta)^{-1})}{d(j)|\delta ~ N(0, diag(delta)^(-1))} 
where \eqn{\alpha = (\alpha_1, ..., \alpha_K)}{alpha = (alpha_1, ..., alpha_K)},
\eqn{\beta = (\beta_1, ..., \beta_P)}{beta = (beta_1, ..., beta_P)},
\eqn{\gamma = (\gamma_1, ..., \gamma_K)}{gamma = (gamma_1, ..., gamma_K)}, and
\eqn{\delta = (\delta_1, ..., \delta_K)}{delta = (delta_1, ..., delta_K)}.

An EM-like algorithm is used to estimate the hyperparameters in an iterative procedure conditioned on current
estimates of the hidden states.

\code{conv.1}, \code{conv.2}, and \code{conv.3} correspond to convergence criteria \eqn{\Delta_1}, 
\eqn{\Delta_2}, and \eqn{\Delta_3} in the reference below, respectively.  After terminating the algorithm, the z-scores of the
\eqn{C} and \eqn{D} matrices can be calculated by inputting \code{CPost} and \code{CvarPost} or \code{DPost} and
\code{DvarPost}, respectively, into \code{\link{zcutoff}}.  This in turn determines the presence or absence
of edges in the network.

See the reference below for additional details about the implementation of the algorithm.
}
\value{
  \item{APost }{Posterior mean of matrix \eqn{A}}
  \item{BPost }{Posterior mean of matrix \eqn{B}}
  \item{CPost }{Posterior mean of matrix \eqn{C}}
  \item{DPost }{Posterior mean of matrix \eqn{D}}
  \item{CvarPost }{Posterior variance of matrix \eqn{C}}
  \item{DvarPost }{Posterior variance of matrix \eqn{D}}
  \item{xPost }{Posterior mean of hidden states \eqn{x}}
  \item{alphaEst }{Estimated value of \eqn{\alpha}{alpha}}
  \item{betaEst }{Estimated value of \eqn{\beta}{beta}}
  \item{gammaEst }{Estimated value of \eqn{\gamma}{gamma}}
  \item{deltaEst }{Estimated value of \eqn{\delta}{delta}}
  \item{vEst }{Estimated value of precisions \eqn{V}}
  \item{muEst }{Estimated value of \eqn{\mu}{mu} (hidden state means)}
  \item{sigmaEst }{Estimated value of \eqn{\Sigma}{Sigma} (hidden state covariance matrix)}
  \item{alliterations }{Total number of iterations run}
}
\references{Andrea Rau, Florence Jaffrezic, Jean-Louis Foulley, and R. W. Doerge (2009). An Empirical Bayesian
Method for Estimating Biological Networks from Temporal Microarray Data. 
\emph{Statistical Applications in Genetics and Molecular Biology} (submitted).}
\author{Andrea Rau}
\seealso{\code{\link{hankel}}, \code{\link{zcutoff}}}
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


####################################################
## Run EB-DBN
####################################################

K <- 0

## Choose alternative value of K using hankel if hidden states are to be estimated
# lag1 <- hankel(y, R, T, P, lag = 1, cutoff = 0.90, type = "median")$dim
# K <- lag1

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

## Run algorithm, if K > 0 set maxiter > 1 (say 100 or 200)
net <- ebdbn(y, P, K, R, T, x.0, alpha.0, beta.0, gamma.0,
        delta.0, v.0, mu.0, sigma.0, conv.1 = 0.15,
        conv.2 = 0.10, conv.3 = 0.10, maxiter = 1)

## Calculate sensitivities, specificities, and precisions of D matrix
## Use z-score significance level of 95%

z <- zcutoff(net$DPost, net$DvarPost, R)
sens.95 <- sensitivity(Dtrue, z$z95, P)  ## 7 TP, 0 FP, 3 FN, 90 TN
}

\keyword{methods}
\keyword{models}
\keyword{ts}
