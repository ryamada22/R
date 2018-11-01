library(igraph)
g <- make_tree(10, 3, mode = "undirected")
n.edge <- length(E(g))
e.len <- runif(n.edge)

g$weight <- e.len

nv <- length(V(g))
paths <- list()
for(i in 1:nv){
	paths[[i]] <- all_shortest_paths(g,from = i)
}

dep.array <- array(0,rep(nv,3))
paths.mat <- matrix(0,0,nv)

for(i in 1:(nv-1)){
	for (j in (i+1):nv){
		tmp <- all_shortest_paths(g,from=i,to=j)[[1]]
		dep.array[i,j,unlist(tmp)] <- 1
		if(length(unlist(tmp)>2)){
			tmp2 <- rep(0,nv)
			tmp2[unlist(tmp)] <- 1
			paths.mat <- rbind(paths.mat,tmp2)
		}
	}
}

incl <- matrix(0,length(paths.mat[,1]),length(paths.mat[,1]))

for(i in 1:length(incl[,1])){
	for(j in 1:length(incl[,1])){
		tmp <- paths.mat[i,] - paths.mat[j,]
		if(prod(tmp <= 0)){
			incl[i,j] <- 1
		}
	}
}
diag(incl) <- 0
independent.paths <- paths.mat[which(apply(incl,1,sum) == 0),]
independent.paths

# ‚±‚ê‚Ì•¡ŽG‚³‚Æ‚©‚ð
