`zCutoff` <-
function(Post, varPost) 
{
    	if (is.matrix(Post) != TRUE) {
        	stop("Error: ", paste(sQuote("Post"), sep = ""), " must be a matrix.")
    	}
	dim1 <- dim(Post)[1]
	dim2 <- dim(Post)[2]
	P <- length(varPost)

	var <- matrix(0, nrow = dim1, ncol = dim2)
	for(i in 1:P) {
		var[i,] <- diag(varPost[[i]])
	}
	sd <- sqrt(var)
	z <- Post / sd
	z95 <- z99 <- z99.9 <- matrix(0, nrow = dim1, ncol = dim2)
	z95[which(abs(z) > 1.96, arr.ind = TRUE)] = 1
	z99[which(abs(z) > 2.58, arr.ind = TRUE)] = 1
	z99.9[which(abs(z) > 3.30, arr.ind = TRUE)] = 1

	return(list(z = z, z95 = z95, z99 = z99, z99.9 = z99.9))
}

