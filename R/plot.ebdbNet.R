`plot.ebdbNet` <-
function (x, sig.level, interactive = "FALSE", ...) 
{
	
	ebdbn <- x
    	if (class(ebdbn) != "ebdbNet") {
        	stop("Error: ", paste(sQuote("ebdbn"), sep = ""), " must be of class ", 
            paste(dQuote("ebdbNet"), sep = ""), sep = "")
    	}

	## Choose edges based on user-defined significance level
	ebdbn.z <- ebdbn$z
	ebdbn.net <- matrix(0, nrow = nrow(ebdbn.z), ncol = ncol(ebdbn.z))
	cutoff <- qnorm((1+sig.level)/2)
	ebdbn.net[which(abs(ebdbn.z) > cutoff, arr.ind = TRUE)] <- 1 

	## Feedback graphs
	if(dim(ebdbn.net)[1] == dim(ebdbn.net)[2]) {
		if(is.null(rownames(ebdbn.net)) == TRUE) rownames(ebdbn.net) <- 1:dim(ebdbn.net)[1];
		ebdbn.igraph <- graph.adjacency(ebdbn.net, mode=c("directed"), add.rownames = TRUE)
		if(interactive == FALSE) {
			plot(ebdbn.igraph, layout = layout.circle, vertex.label = c(rownames(ebdbn.net)), ...)
		}
		if(interactive == TRUE) {
			tkplot(ebdbn.igraph, vertex.label = c(rownames(ebdbn.net)), ...)
		}
	}


	## Input graphs
	if(dim(ebdbn.net)[1] != dim(ebdbn.net)[2]) {
		if(is.null(rownames(ebdbn.net)) == TRUE) rownames(ebdbn.net) <- 1:dim(ebdbn.net)[1];
		if(is.null(colnames(ebdbn.net)) == TRUE) colnames(ebdbn.net) <- paste("u", 1:dim(ebdbn.net)[2], sep = "");

		tmp <- which(ebdbn.net != 0, arr.ind = TRUE)[,2:1]
		tmp[,1] <- tmp[,1] + dim(ebdbn.net)[1]
		ebdbn.igraph <- graph.bipartite(types = c(rep(0, nrow(ebdbn.net)), rep(1, ncol(ebdbn.net))),
			edges = (matrix(t(tmp), nrow = 1)-1), directed = TRUE)


		if(interactive == FALSE) {
			plot(ebdbn.igraph, layout = layout.circle, vertex.color = c(rep("SkyBlue2", nrow(ebdbn.net)),
				rep("green", ncol(ebdbn.net))), vertex.label = c(rownames(ebdbn.net), 
				colnames(ebdbn.net)), ...)
		}
		if(interactive == TRUE) {
			tkplot(ebdbn.igraph, vertex.label = c(rownames(ebdbn.net), colnames(ebdbn.net)), ...)
		}

	}

}
