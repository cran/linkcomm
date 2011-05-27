# Function(s) to extract link communitites in directed, undirected, weighted, or unweighted networks so that we can cluster nodes, allowing them to belong to multiple different communities.
#
# Author: Alex T. Kalinka (alex.t.kalinka@gmail.com)
#
# See: Ahn et al. (2010). Link communities reveal multiscale complexity in networks. Nature 466:761-765.


.onLoad <- function(lib, pkg) 
	{
	require(utils)
	packageStartupMessage("\nWelcome to linkcomm version ",packageDescription(pkg)$Version,"\n\nFor a step-by-step guide to using linkcomm functions:\n   > vignette(topic = \"linkcomm\", package = \"linkcomm\")\nTo run an interactive demo:\n   > demo(topic = \"linkcomm\", package = \"linkcomm\")\nTo cite, see:\n   > citation(\"linkcomm\")\nNOTE: To use linkcomm, you require read and write permissions in the current directory (see: help(\"getwd\"), help(\"setwd\"))\n")
	}


.getLoopInds <- function(x) # Returns 1 if an edge is a loop, 2 otherwise. 
	{
	return(length(unique(x)))
	}


.removeEdgeDuplicates <- function(x, weights = NA) # x is an edge list.
	{ 
	x <- cbind(as.character(x[,1]),as.character(x[,2]))
	inds <- NULL
	if(length(which(duplicated(x)==TRUE))>0){
		inds <- append(inds,which(duplicated(x)==TRUE))
		}
	if(length(inds)>0){
		x <- x[-inds,]
		if(length(weights)>1){
			weights <- weights[-inds]
			}
		inds <- NULL
		}
	# Also remove duplicates that are in the opposite orientation.
	xx <- rbind(x,cbind(x[,2],x[,1]))
	if(length(which(duplicated(xx)==TRUE))>0){
		dd <- which(duplicated(xx)==TRUE) - nrow(x)
		inds <- append(inds,dd)
		}
	if(length(inds)>0){
		x <- x[-inds,]
		if(length(weights)>1){
			weights <- weights[-inds]
			}
		}
	if(length(weights)>1){
		ret <- list()
		ret$edges <- x
		ret$weights <- weights
		return(ret)
	}else{
		return(x)
		}
	}


getLinkCommunities <- function(x, hcmethod = "average", edglim = 10^4, directed = FALSE, dirweight = 0.5, plot = TRUE, removetrivial = TRUE, verbose = TRUE) 
	# x is an edge list. Nodes can be ASCII names or integers, but are always treated as character names in R.
	# If plot is true (default), a dendrogram and partition density score as a function of dendrogram height are plotted side-by-side.
	# When there are more than "edglim" edges, hierarchical clustering is carried out via temporary files written to disk using compiled C++ code.
	{		
	if(ncol(x)==3){
		xx <- cbind(as.character(x[,1]),as.character(x[,2]))
		wt <- as.numeric(as.character(x[,3]))
	}else if(ncol(x)==2){
		xx <- cbind(as.character(x[,1]),as.character(x[,2]))
		wt <- NA
	}else{
		stop("\nInput data must be an edge list.\n")
		}
	if(verbose){
		cat("   Removing loops...\n")
		}
	loops <- apply(xx,1,.getLoopInds)
	if(length(which(loops==1))>0){
		xx <- xx[-which(loops==1),]
		if(length(wt)>1){
			wt <- wt[-which(loops==1)]
			}
		}
	if(verbose){
		cat("   Removing edge duplicates...\n")
		}
	if(length(wt)>1){
		yy <- .removeEdgeDuplicates(xx, weights = wt)
		xx <- yy$edges
		wt <- yy$weights
	}else{
		xx <- .removeEdgeDuplicates(xx)
		}

	el <- xx # Modified edge list returned to user.
	len <- nrow(xx) # Number of edges.
	nnodes <- length(unique(c(as.character(xx[,1]),as.character(xx[,2])))) # Number of nodes.
	myedges <- xx # To be used later.

	xx <- graph.edgelist(xx, directed = directed) # Creates "igraph" object.
	edges <- cbind(xx[[3]],xx[[4]]) # Edges with numerical node IDs

	# Can we use carriage returns in our progress indicators?
	if(.Platform$OS.type == "unix"){
		carriageret <- TRUE
	}else{
		carriageret <- FALSE
		}

	# Switch depending on size of network.
	if(len <= edglim){
		disk <- FALSE
		emptyvec <- rep(1,(len*(len-1))/2)
		if(length(wt)>1){ weighted <- TRUE}else{ wt <- 0; weighted <- FALSE}
		dissvec <- .C("getEdgeSimilarities",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(len),rowlen=integer(1),weights=as.double(wt),as.logical(directed),as.double(dirweight),as.logical(weighted),as.logical(disk), dissvec = as.double(emptyvec), as.logical(carriageret), as.logical(verbose))$dissvec
		distmatrix <- matrix(1,len,len)
		distmatrix[lower.tri(distmatrix)] <- dissvec
		colnames(distmatrix) <- 1:len
		rownames(distmatrix) <- 1:len
		distobj <- as.dist(distmatrix) # Convert into 'dist' object for hclust.
		rm(distmatrix)
		if(verbose){
			cat("\n   Hierarchical clustering of edges...")
			}
		hcedges <- hclust(distobj, method = hcmethod)
		hcedges$order <- rev(hcedges$order)
		if(verbose){cat("\n")}
		#return(hcedges)
	}else{
		disk <- TRUE
		if(length(wt)>1){ weighted <- TRUE}else{ wt <- 0; weighted <- FALSE}
		rowlen <- .C("getEdgeSimilarities",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(len),rowlen=integer(len-1),weights=as.double(wt),as.logical(directed),as.double(dirweight),as.logical(weighted),as.logical(disk), dissvec = double(1), as.logical(carriageret), as.logical(verbose))$rowlen
		#return(rowlen)
		if(verbose){cat("\n")}
		hcobj <- .C("hclustLinkComm",as.integer(len),as.integer(rowlen),heights = single(len-1),hca = integer(len-1),hcb = integer(len-1), as.logical(carriageret), as.logical(verbose))
		if(verbose){cat("\n")}
		hcedges<-list()
		hcedges$merge <- cbind(hcobj$hca, hcobj$hcb)
		hcedges$height <- hcobj$heights
		#return(hcedges)
		hcedges$order <- .C("hclustPlotOrder",as.integer(len),as.integer(hcobj$hca),as.integer(hcobj$hcb),order=integer(len))$order
		hcedges$order <- rev(hcedges$order)
		hcedges$method <- "single"
		class(hcedges) <- "hclust"
		#return(hcedges)
		}

	hh <- unique(round(hcedges$height, digits = 5)) # Round to 5 digits to prevent numerical instability affecting community formation.
	countClusters <- function(x,ht){return(length(which(ht==x)))}
	clusnums <- sapply(hh, countClusters, ht = round(hcedges$height, digits = 5)) # Number of clusters at each height.
	ldlist <- .C("getLinkDensities",as.integer(hcedges$merge[,1]), as.integer(hcedges$merge[,2]), as.integer(edges[,1]), as.integer(edges[,2]), as.integer(len), as.integer(clusnums), pdens = double(length(hh)), heights = as.double(hh), pdmax = double(1), csize = integer(1), as.logical(removetrivial), as.logical(carriageret), as.logical(verbose))
	pdens <- c(0,ldlist$pdens)
	heights <- c(0,hh)
	pdmax <- ldlist$pdmax
	csize <- ldlist$csize

	if(csize == 0){
		stop("\nThere are no clusters appearing in this network; maybe try a larger network.\n")
		}

	if(verbose){
		cat("\n   Partition density maximum height = ",pdmax,"\n")
		}

	# Read in optimal clusters from a file.
	clus <- list()
	for(i in 1:csize){
		if(verbose){
			mes<-paste(c("   Finishing up...1/4... ",floor((i/csize)*100),"%"),collapse="")
			cat(mes,"\r")
			flush.console()
			}
		clus[[i]] <- scan(file = "linkcomm_clusters.txt", nlines = 1, skip = i-1, quiet = TRUE)
		}

	file.remove("linkcomm_clusters.txt")

	# Extract nodes for each edge cluster.
	ecn <- data.frame()
	ee <- data.frame()
	for(i in 1:(length(clus))){
		if(verbose){
			mes<-paste(c("   Finishing up...2/4... ",floor((i/length(clus))*100),"%"),collapse="")
			cat(mes,"\r")
			flush.console()
			}
		ee <- rbind(ee,cbind(myedges[clus[[i]],],i))
		nodes <- V(xx)$name[(unique(c(edges[clus[[i]],]))+1)]
		both <- cbind(nodes,rep(i,length(nodes)))
		ecn <- rbind(ecn,both)
		}
	colnames(ecn) <- c("node","cluster")

	# Extract the node-size of each edge cluster and order largest to smallest.
	ss <- NULL
	unn <- unique(ecn[,2])
	for(i in 1:length(unn)){
		if(verbose){
			mes<-paste(c("   Finishing up...3/4... ",floor((i/length(unn))*100),"%"),collapse="")
			cat(mes,"\r")
			flush.console()
			}
		ss[i] <- length(which(ecn[,2]==unn[i]))
		}
	names(ss) <- unn
	ss <- sort(ss,decreasing=T)

	# Extract the number of edge clusters that each node belongs to.
	oo<-NULL
	unn <- unique(ecn[,1])
	for(i in 1:length(unn)){
		if(verbose){
			mes<-paste(c("   Finishing up...4/4... ",floor((i/length(unique(ecn[,1])))*100),"%"),collapse="")
			cat(mes,"\r")
			flush.console()
			}
		oo[i]<-length(which(ecn[,1]==unn[i]))
		}
	if(verbose){cat("\n")}
	names(oo)<-unn
	
	pdplot <- cbind(heights,pdens)

	# Add nodeclusters of size 0.
	missnames <- setdiff(V(xx)$name,names(oo))
	m <- rep(0,length(missnames))
	names(m) <- missnames
	oo <- append(oo,m)

	all <- list()

	all$numbers <- c(len,nnodes,length(clus)) # Number of edges, nodes, and clusters.
	all$hclust <- hcedges # Return the 'hclust' object. To plot the dendrogram: 'plot(lcobj$hclust,hang=-1)'
	all$pdmax <- pdmax # Partition density maximum dendrogram height.
	all$pdens <- pdplot # Add data for plotting Partition Density as a function of dendrogram height.
	all$nodeclusters <- ecn # n*2 character matrix of node names and the cluster ID they belong to.
	all$clusters <- clus # Clusters of edge IDs arranged as a list of lists.
	all$edges <- ee # Edges and the clusters they belong to, arranged so we can easily put them into an edge attribute file for Cytoscape.
	all$numclusters <- sort(oo,decreasing=TRUE) # The number of clusters that each node belongs to (named vector where the names are node names).
	all$clustsizes <- ss # Cluster sizes sorted largest to smallest (named vector where names are cluster IDs).
	all$igraph <- xx # igraph graph.
	all$edgelist <- el # Edge list.
	all$directed <- directed # Logical indicating if graph is directed or not.

	class(all) <- "linkcomm"

	if(plot){
		if(verbose){
			cat("   Plotting...\n")
			}
		if(len < 1500){ # Will be slow to plot dendrograms for large networks.
			if(len < 500){
				all <- plot(all, type="summary", verbose = verbose)
			}else{
				all <- plot(all, type="summary", right = FALSE, verbose = verbose) # Slow to reverse order of large dendrograms.
				}
		}else if(len <= edglim){
			par(mfrow=c(1,2), mar=c(5.1,4.1,4.1,2.1))
			plot(hcedges,hang=-1,labels=FALSE)
			abline(pdmax,0,col='red',lty=2)
			plot(pdens,heights,type='n',xlab='Partition Density',ylab='Height')
			lines(pdens,heights,col='blue',lwd=2)
			abline(pdmax,0,col='red',lty=2)
		}else{
			plot(heights,pdens,type='n',xlab='Height',ylab='Partition Density')
			lines(heights,pdens,col='blue',lwd=2)
			abline(v = pdmax,col='red',lwd=2)
			}
		}

	return(all)
	
	}



