`visualize` <-
function (zscores, sig, format = "Cytoscape", type = "feedback") 
{
    if (format != "Cytoscape" & format != "Rgraphviz") {
        stop("Error: ", paste(sQuote("format"), sep = ""), " must be equal to ", 
            paste(dQuote("Cytoscape"), sep = ""), " or ", paste(dQuote("Rgraphviz"), 
                sep = ""))
    }

    if (type != "feedback" & type != "input") {
	stop("Error: ", paste(sQuote("type"), sep = ""), " must be equal to ",
	    paste(dQuote("feedback"), sep = ""), " or ", paste(dQuote("input"),
		sep = ""))
    }	


    if (format == "Cytoscape") {
        row <- dim(zscores)[1]
        col <- dim(zscores)[2]
        network <- matrix(0, nrow = row, ncol = col)
        network <- sign(zscores) * sig
        pos <- matrix(which(network == 1, arr.ind = TRUE)[, 2:1], 
            ncol = 2)
        neg <- matrix(which(network == -1, arr.ind = TRUE)[, 
            2:1], ncol = 2)
        pos.edges <- cbind(pos, rep(1, length(which(network == 
            1))))
        neg.edges <- cbind(neg, rep(-1, length(which(network == 
            -1))))
        edges <- rbind(pos.edges, neg.edges)
	if(type == "input") {
		edges <- cbind(paste("TF", edges[,1], sep = ""),
			paste("G", edges[,2], sep = ""),
			edges[,3])
	}
        return(edges)
    }

  if (format == "Rgraphviz") {
	print("Rgraphviz package has been removed from CRAN repository.")
#        require(Rgraphviz)
#        graph.par(list(nodes = list(fill = "lightgray", lwd = 2, 
#            fontsize = 10, textCol = "black")))
#        edges.full <- which(sig != 0, arr.ind = TRUE)
# 	if(type == "input") {
# 		edges <- cbind(paste("TF", edges.full[,2], sep = ""),
# 			paste("G", edges.full[,1], sep = ""))
# 	}
# 
#        V <- as.character(unique(as.vector(edges)))
#        gR <- new("graphNEL", nodes = V, edgemode = "directed")
#        gX <- addEdge(edges[,1], edges[,2], gR, abs(zscores[edges.full]))
#        nd <- nodes(gX)[degree(gX)[[1]] > 0 | degree(gX)[[2]] > 0]
#        gX <- subGraph(nd, gX)
#        edgeRenderInfo(gX) <- list(arrowhead = "normal", arrowtail = "none", 
#             lwd = 2)
#        gX2 <- layoutGraph(gX)
#        renderGraph(gX2)
     	}
}
