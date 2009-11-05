`ebdbn` <-
function(y, P, K, R, T, x.0, alpha.0, beta.0, gamma.0, delta.0, v.0,
mu.0, sigma.0, conv.1 = .15, conv.2 = .05, conv.3 = .01, maxiter = 100) 
{
if(is.list(y) != TRUE) {
stop("Error: ", paste(sQuote("y"), sep = ""), " should be a list.")
}
if(K > 0) {
if(length(alpha.0) != K | length(beta.0) != P | length(gamma.0) != K | length(delta.0) != P |
length(v.0) != P) {
stop("Error: Initial values of hyperparameters are not the correct length.")
}
}
if(K == 0) {
if(length(delta.0) != P | length(v.0) != P) {
stop("Error: Initial values of hyperparameters are not the correct length.")
}
}
APost <- rep(0, K*K)
BPost <- rep(0, P*K)
CPost <- rep(0, P*K)
CvarPost <- rep(0, P*K*K)
DPost <- rep(0, P*P)
DvarPost <- rep(0, P*P*P)
if(K>0) {x0 <- as.vector(t(do.call(rbind, x.0)))}
if(K==0) {x0 <- 0}
yorig <- as.vector(t(do.call(rbind, y)))

## Run Full Algorithm and Find Posterior Mean and Variance

## Load "RunWrap" C code
## Use RunWrap.dll for Windows, RunWrap.so for Unix
## dyn.load("RunWrap.dll")
test = .C("RunWrap", R = as.integer(R), P = as.integer(P), 
T = as.integer(T), K = as.integer(K), xx = as.double(x0), 
yy = as.double(yorig), alpha = as.double(alpha.0), 
beta = as.double(beta.0), gamma = as.double(gamma.0), 
delta = as.double(delta.0), v = as.double(v.0),
mu = as.double(mu.0), sigma = as.double(sigma.0),
conv1 = as.double(conv.1), conv2 = as.double(conv.2),
conv3 = as.double(conv.3), APost = as.double(APost),
BPost = as.double(BPost), CPost = as.double(CPost),
DPost = as.double(DPost), CvarPost = as.double(CvarPost),
DvarPost = as.double(DvarPost),
alliterations = as.integer(0), 
maxiterations = as.integer(maxiter),
subiterations = as.integer(200), PACKAGE = "ebdbNet")

## Estimated Posterior Mean and Variance

APost <- matrix(test$APost, nrow = K, ncol = K, byrow = TRUE)
BPost <- matrix(test$BPost, nrow = K, ncol = P, byrow = TRUE)
CPost <- matrix(test$CPost, nrow = P, ncol = K, byrow = TRUE)
DPost <- matrix(test$DPost, nrow = P, ncol = P, byrow = TRUE)
CvarPost <- vector("list", P)
DvarPost <- vector("list", P)
for(i in 1:P) {
DvarPost[[i]] <- matrix(test$DvarPost[(P*P*(i-1)+1):(P*P*i)], nrow = P, ncol = P, byrow = TRUE)
CvarPost[[i]] <- matrix(test$CvarPost[(K*K*(i-1)+1):(K*K*i)], nrow = K, ncol = K, byrow = TRUE)
}

## Estimated hidden states
xPost <- vector("list", R)
for(r in 1:R) {
xPost[[r]] <- matrix(test$xx[(K*T*(r-1)+1):(K*T*r)], nrow = K, ncol = T, byrow = TRUE)
}

## Estimated hyperparameters
alphaEst <- test$alpha
betaEst <- test$beta
gammaEst <- test$gamma
deltaEst <- test$delta
vEst <- test$v
muEst <- test$mu
sigmaEst <- test$sigma

## Total number of iterations needed
alliterations <- test$alliterations

return(list(APost = APost, BPost = BPost, CPost = CPost, DPost = DPost, CvarPost = CvarPost, 
DvarPost = DvarPost, xPost = xPost, alphaEst = alphaEst, betaEst = betaEst, 
gammaEst = gammaEst, deltaEst = deltaEst, vEst = vEst, muEst = muEst, sigmaEst = sigmaEst,
alliterations = alliterations))
}

