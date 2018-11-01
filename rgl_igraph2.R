library(rgl) # package for 3d object handling
# reading the bunny.obj in your fair zip
bunny <- readOBJ("bunny.obj")

library(igraph) # package for graph theory

# 3D coordinates of vertices in the shape of n x 3 matrix
V.xyz <- t(bunny[[1]][1:3,])

# Enumerate edges of triangle in n x 2 matrix shape
Edges <- rbind(t(bunny[[2]][1:2,]),t(bunny[[2]][2:3,]),t(bunny[[2]][c(3,1),]))
# Remove duplicates of edges
Edges <- t(apply(Edges,1,sort))
Edges <- unique(Edges)
# length of edges
Edge.len <- sqrt(apply((V.xyz[Edges[,1],] - V.xyz[Edges[,2],])^2,1,sum))
# make a graph object
g <- graph.edgelist(Edges,directed =FALSE)
#gt <- graph.edgelist(Edges)
# distance on the graph
#dt <- distances(gt,weights=Edge.len)
d <- distances(g,weights=Edge.len)


### Post-spherization
# Spherization maps all the vertices on the bunny on the S2 sphere
# Along the geodesics on the S2 sphere, we can "re-measure" the distance.
# The geodesics on the S2 sphere is drawn back on the bunny and
# along the geodesics, we can measure the "distance between vertices on the bunny".
# We can make the distance matrix of all vertices pairs along this back-to-the-bunny geodesics.
# Can you compare the two distance matrices; one is the distance matrix that the R codes above generates and the other is the back-to-the-bunny distance matrix.
# Showing distance from a point.

# Normalization of distance value ranging [0,1] so that you can use the values for color function rgb arguments
d. <- d/max(d)
# Coloring the vertices with distance from the n-th vertex
# Select one of (many vertices) arbitrarily
# The color of bunny should depend on the distance from the n-th vertex.
n <- 30
# rgb() function takes three [0,1] real values specifying red/green/blue
# to generate colors
# d.[n,] is the vector of normalized distance from the n-th vertex to all vertices
col <- rgb(d.[n,],1-d.[n,],1)
# Integer values >= 1 are generated depending on distance
col2 <- ceiling(d.[n,] * 15)+1
# plot them
plot3d(V.xyz)
spheres3d(V.xyz,radius=0.005,col=col2)

# When you get different distance matrix
# you can replace the object d to the new distance matrix
# and draw the similar contours.
# Then, you can visually compare two distance matrices.

n <- 30
col <- rgb(d.[n,],1-d.[n,],1)
col2 <- ceiling(d.[n,] * 15)+1
# color values are 1 or 2
col3 <- col2 %% 2 + 1
plot3d(V.xyz)
spheres3d(V.xyz,radius=0.005,col=col3)

# Change the vertex, distance form which determines colors

n2 <- 400
col_ <- rgb(d.[n2,],1-d.[n2,],1)
col2_ <- ceiling(d.[n2,] * 15)+1
# color values are 1 or 2
col3_ <- col2_ %% 2 + 1
plot3d(V.xyz)
spheres3d(V.xyz,radius=0.005,col=col3_)


# combining two distance vectors from two vertices and generate color values 1 or 2
# This should make square lattice type colorings
plot3d(V.xyz)
spheres3d(V.xyz,radius=0.005,col=as.numeric(col3 == col3_)+1)

my.tips <- function(g,v,d,id=FALSE){
	n <- length(d[1,])
	adj <- adjacent_vertices(g, 1:n)
	ret <- rep(0,n)
	for(i in 1:n){
		tmp <- c(i,adj[[i]])
		if(max(d[v,tmp]) == d[v,i]){
			ret[i] <- 1
		}
	}
	if(id){
		ret <- which(ret==1)
	}
	return(ret)
}

my.graph.voronoi <- function(d,tips.id){
	tmp.d <- d[,tips.id]
	col.region <- apply(tmp.d,1,which.min)
	return(col.region)
}

my.scale.rgb <- function(x){
	x. <- (x-min(x))/(max(x)-min(x))
	return(rgb(x.,1-x.,1))
}

my.heat.colors <- function(x){
	n <- max(x)
	return(heat.colors(n)[x])
}

my.rainbow.colors <- function(x){
	n <- max(x)
	return(rainbow(n)[x])
}

v <- 1000
tips.id <- my.tips(g,v,d,id=TRUE)

selected.tip <- tips.id[3]
tips.id <- my.tips(g,tips.id[3],d,id=TRUE)

tips.id2 <- c(selected.tip,tips.id)

tmp.d <- d[,tips.id]

col.region <- apply(tmp.d,1,which.min)

col.region <- my.graph.voronoi(d,tips.id)
col.region2 <- my.graph.voronoi(d,tips.id2)

#col.region <- my.scale.rgb(col.region)
#col.region2 <- my.scale.rgb(col.region2)
col.region <- my.rainbow.colors(col.region)
col.region2 <- my.rainbow.colors(col.region2)

plot3d(V.xyz)
spheres3d(V.xyz[v,],radius=0.01,col=1)
spheres3d(V.xyz[tips.id2,],radius=0.01,col=my.scale.rgb(1:length(tips.id2)))
spheres3d(V.xyz,radius=0.005,col=col.region2)


plot3d(rbind(V.xyz,V.xyz+max(V.xyz)))
spheres3d(V.xyz[v,]+max(V.xyz),radius=0.01,col=1)
spheres3d(V.xyz[tips.id,]+max(V.xyz),radius=0.01,col=2)
spheres3d(V.xyz[v,],radius=0.01,col=1)
spheres3d(V.xyz[tips.id,],radius=0.01,col=1:length(tips.id)+1)
spheres3d(V.xyz,radius=0.005,col=col.region)


my.paths2tips <- function(g,v,tips,edge.len){
	ret <- list()
	for(i in 1:length(tips)){
		ret[[i]] <- shortest_paths(g,v,tips[i],weights=edge.len)[[1]][[1]]
	}
	return(ret)
}

my.tipPathGraph <- function(g,v,d,edge.len){
	tips <- which(my.tips(g,v,d)==1)
	pt <- my.paths2tips(g,v,tips,edge.len)
	el <- matrix(0,0,2)
	for(i in 1:length(pt)){
		tmp <- pt[[i]]
		tmp2 <- tmp[1:(length(tmp)-1)]
		tmp3 <- tmp[2:length(tmp)]
		tmp4 <- cbind(tmp2,tmp3)
		el <- rbind(el,tmp4)
	}
	el <- t(apply(el,1,sort))
	el <- unique(el)
	#ret.g <- graph.edgelist(el,directed=FALSE)
	return(list(el=el,tips=tips))
}

v <- 100

tpg <- my.tipPathGraph(g,v,d,Edge.len)
el. <- matrix(paste("",tpg$el),ncol=2)
gg <- graph.edgelist(el.,directed=FALSE)
plot(gg)


plot(tpg$g)


v <- 100
col_ <- rgb(d.[v,],1-d.[v,],1)
col2_ <- ceiling(d.[v,] * 15)+1
# color values are 1 or 2
col3_ <- col2_ %% 2 + 1
plot3d(V.xyz)
spheres3d(V.xyz,radius=0.005,col=col3_)

tips <- my.tips(g,v,d)
col4 <- rep(1,length(d[1,]))
col4[v] <- 2
col4[which(tips==1)] <- 3
plot3d(V.xyz)
spheres3d(V.xyz[which(tips==0),],radius=0.00001,col=1)
spheres3d(V.xyz[which(tips==1),],radius=0.01,col=3)
spheres3d(V.xyz[v,], radius=0.01,col=2)

sh.path <- all_shortest_paths(g,v,weights=Edge.len,mode="out")
len.path <- sapply(sh.path[[1]],length)
n.path.edge <- sum(len.path)-length(len.path)

path.edge <- matrix(0,2,n.path.edge)
cnt <- 1
#for(i in 1:3){
for(i in 1:length(d[1,])){
	tmp <- sh.path[[1]][[i]]
	#tmp <- shortest_paths(g,v,i,weights=Edge.len)[[1]][[1]]
	if(length(tmp)>1){
		tmp2 <- tmp[1:(length(tmp)-1)]
		tmp3 <- tmp[2:length(tmp)]
		tmp4 <- rbind(tmp2,tmp3)
		#path.edge <- cbind(path.edge,tmp4)
		path.edge[,cnt:(cnt+length(tmp2)-1)] <- tmp4
		cnt <- cnt+length(tmp2)
		#path.edge <- unique(path.edge)
		#segments3d(V.xyz[c(tmp4),])
	}
}
path.edge2 <- unique(t(path.edge))
segments3d(V.xyz[c(t(path.edge2)),])

plot3d(V.xyz,size=0.01)
spheres3d(V.xyz[which(tips==0),],radius=0.0000001,col=gray(0.1))
spheres3d(V.xyz[which(tips==1),],radius=0.001,col=3)
spheres3d(V.xyz[v,], radius=0.001,col=2)

tips.id <- which(tips==1)
for(i in 1:length(tips.id)){
	tmp <- sh.path[[1]][[tips.id[i]]]
	#tmp <- shortest_paths(g,v,i,weights=Edge.len)[[1]][[1]]
	if(length(tmp)>1){
		tmp2 <- tmp[1:(length(tmp)-1)]
		tmp3 <- tmp[2:length(tmp)]
		tmp4 <- rbind(tmp2,tmp3)
		#path.edge <- cbind(path.edge,tmp4)
		#path.edge[,cnt:(cnt+length(tmp2)-1)] <- tmp4
		#cnt <- cnt+length(tmp2)
		#path.edge <- unique(path.edge)
		segments3d(V.xyz[c(tmp4),])
	}
}
#mst <- mst(g)

d.tips <- d.[tips.id,tips.id]

library(ape)
tip.tr <- nj(d.tips)
plot(tip.tr)

v <- tips.id[1]
col_ <- rgb(d.[v,],1-d.[v,],1)
col2_ <- ceiling(d.[v,] * 15)+1
# color values are 1 or 2
col3_ <- col2_ %% 2 + 1
plot3d(V.xyz)
spheres3d(V.xyz,radius=0.005,col=col3_)

tips <- my.tips(g,v,d)
col4 <- rep(1,length(d[1,]))
col4[v] <- 2
col4[which(tips==1)] <- 3
plot3d(V.xyz)
spheres3d(V.xyz[which(tips==0),],radius=0.00001,col=1)
spheres3d(V.xyz[which(tips==1),],radius=0.01,col=3)
spheres3d(V.xyz[v,], radius=0.01,col=2)

sh.path <- all_shortest_paths(g,v,weights=Edge.len,mode="out")
len.path <- sapply(sh.path[[1]],length)
n.path.edge <- sum(len.path)-length(len.path)

path.edge <- matrix(0,2,n.path.edge)
cnt <- 1
#for(i in 1:3){
for(i in 1:length(d[1,])){
	tmp <- sh.path[[1]][[i]]
	#tmp <- shortest_paths(g,v,i,weights=Edge.len)[[1]][[1]]
	if(length(tmp)>1){
		tmp2 <- tmp[1:(length(tmp)-1)]
		tmp3 <- tmp[2:length(tmp)]
		tmp4 <- rbind(tmp2,tmp3)
		#path.edge <- cbind(path.edge,tmp4)
		path.edge[,cnt:(cnt+length(tmp2)-1)] <- tmp4
		cnt <- cnt+length(tmp2)
		#path.edge <- unique(path.edge)
		#segments3d(V.xyz[c(tmp4),])
	}
}
path.edge2 <- unique(t(path.edge))
segments3d(V.xyz[c(t(path.edge2)),])

plot3d(V.xyz,size=0.01)
spheres3d(V.xyz[which(tips==0),],radius=0.0000001,col=gray(0.1))
spheres3d(V.xyz[which(tips==1),],radius=0.001,col=3)
spheres3d(V.xyz[v,], radius=0.001,col=2)

tips.id <- which(tips==1)
for(i in 1:length(tips.id)){
	tmp <- sh.path[[1]][[tips.id[i]]]
	#tmp <- shortest_paths(g,v,i,weights=Edge.len)[[1]][[1]]
	if(length(tmp)>1){
		tmp2 <- tmp[1:(length(tmp)-1)]
		tmp3 <- tmp[2:length(tmp)]
		tmp4 <- rbind(tmp2,tmp3)
		#path.edge <- cbind(path.edge,tmp4)
		#path.edge[,cnt:(cnt+length(tmp2)-1)] <- tmp4
		#cnt <- cnt+length(tmp2)
		#path.edge <- unique(path.edge)
		segments3d(V.xyz[c(tmp4),])
	}
}
#mst <- mst(g)

d.tips <- d.[tips.id,tips.id]

library(ape)
tip.tr <- nj(d.tips)
plot(tip.tr)


