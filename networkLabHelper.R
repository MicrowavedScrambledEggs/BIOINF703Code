library(network)
library(sna)

localclustering <- function(y,low.deg.exclude = TRUE,low.deg.value = 0) {
	## calculates the clustering coefficent for a network y
	N = network.size(y)
	nhoodsize = rep(0,N)
	nedgeinhood = rep(0,N)
	el = as.matrix(y,matrix.type = "edgelist")
	for (i in 1:N) {
		n = get.neighborhood(y,i,"combined")
		nhoodsize[i] = length(n)
		nedgeinhood[i] = sum(el[,1] %in% n &  el[,2] %in% n) #+ nhoodsize[i]
	}
	##print(nhoodsize)
	## determine treatment of degree 0 and 1 nodes
	if (low.deg.exclude){
		## get rid of low degree nodes from calculation
		nedgeinhood = nedgeinhood[nhoodsize > 1]
		nhoodsize = nhoodsize[nhoodsize > 1]
	} else {
		## make size 0/1 n'hoods have local clustering given by low.deg.value
		i = nhoodsize < 2
		nhoodsize[i] = 2
		nedgeinhood[i] = low.deg.value ## change this to 1 to get local clustering 1 for these neighbourhoods
	}
	return(mean(2*nedgeinhood/(nhoodsize *(nhoodsize - 1))))
}


readGraphFromFile <- function(file,returnNames = FALSE){
	
	## read in table
 a = read.table(file)
 
 ## get number of edges
 m = nrow(a)
 
 ## flatten edges
 b = unlist(a)
 
 ## get number of vertices
 n = length(levels(b))
 
 ## make a edge lists matrix
 el <- matrix(c(b),nrow = m,ncol = 2)
 
 emat = matrix(0,nrow = n,ncol = n)
 for (i in 1:m){
 	emat[el[i,1],el[i,2]] = emat[el[i,2],el[i,1]] = 1
  }
 diag(emat) <- 0
 
 g<- network(emat, directed = FALSE)
 
 
 if (returnNames == FALSE){
  return(g)
 } else {
 	return(list(graph = g, vnames = levels(b)))
 }
}

 sampleER <- function(n,p){
   mat = rgraph(n, tprob = p,mode = "graph")
   return(network(mat,directed = FALSE))
 }
 
 rgraphFromDegreeDist <- function(g){
	# constructs a random undirected graph with n = order(g) nodes by resampling from the degree dist of g
	N =  get.network.attribute(g,"n")
	dd = 1
	while (sum(dd) %% 2 == 1) dd = sample(degree(g,gmode = "graph"),N,replace = T)
	nedge = sum(dd)/2
	## build an edgelist
	el = matrix(sample(rep(1:N,times = dd),2*nedge, replace = F),ncol = 2)
	## remove loops
	el = el[el[,1]!=el[,2],]
	## remove repeats - first need common order for all edges
	outorder = el[,1] > el[,2]
	el[outorder,] = el[outorder,2:1]
	return(network(el,directed = FALSE))
}