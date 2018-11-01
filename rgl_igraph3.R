library(rgl) # package for 3d object handling
library(Ronlyryamada)

library(RFOC)
my.cell.sim.2 <- function (n = 5, k = 5, n.mesh = 32, n.step = 50, scale.shift = 0.1, 
    scale.rotation1 = 0.1, scale.rotation2 = 0.1, scale.shear = 0.1, 
    file.name = "hoge") 
{
    A. <- matrix(runif(n^2), n, n)
    A.[1, 1] <- k
    B <- matrix(rnorm(n^2), n, n)
    xxx <- my.spherical.harm.mesh(A = A., B = B, n = n.mesh)
    #plot3d(xxx$v)
    #segments3d(xxx$v[c(t(xxx$edge)), ])
    M <- diag(rep(1, 4))
    M[1:3, 4] <- runif(3) * scale.shift
    theta1 <- runif(1) * scale.rotation1
    theta2 <- runif(1) * scale.rotation2
    Mm <- diag(rep(1, 3))
    Mm[1:2, 1:2] <- my.2d.rot(theta1) %*% Mm[1:2, 1:2]
    Mm[2:3, 2:3] <- my.2d.rot(theta2) %*% Mm[2:3, 2:3]
    M[1:3, 1:3] <- Mm
    M[4, 1:3] <- rnorm(3) * scale.shear
    for (i in 1:n.step) {
        A. <- A. + rnorm(n^2, 0, 0.05)
        xxx <- my.spherical.harm.mesh(A = A., n = n.mesh)
        xxxx <- cbind(xxx$v, rep(1, length(xxx$v[, 1])))
        rot.xxxx <- t(M %*% t(xxxx))
        rot.xxxx. <- rot.xxxx[, 1:3]/rot.xxxx[, 4]
				
        #plot3d(rot.xxxx.)
        #segments3d(rot.xxxx.[c(t(xxx$edge)), ])
        file.out <- paste(file.name, i, ".obj", sep = "")
        if(i == n.step){
					#my.write.obj(rot.xxxx., xxx$f, file.out)
					#open3d()
					#plot3d(rot.xxxx.)
					#segments3d(rot.xxxx.[c(t(xxx$edge)),])
				#writeOBJ(filename)
					#objfilein <- readOBJ(file.out)
					ret <- list(vb = t(rot.xxxx.),it=t(xxx$f))
					#xxx$v <- t(rot.xxxx.)
					#xxx$f <- t(xxx$f)
				}
        
    }
		return(ret)
}
################
# reading the bunny.obj in your fair zip
#bunny <- readOBJ("bunny.obj")
bunny <- my.cell.sim.2()

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


my.path.edge <- function(g,v,Edge.len){
	sh.path <- all_shortest_paths(g,v,weights=Edge.len,mode="out")
	len.path <- sapply(sh.path[[1]],length)
	n.path.edge <- sum(len.path)-length(len.path)
	path.edge <- matrix(0,2,n.path.edge)
	cnt <-1 
	for(i in 1:vcount(g)){
		tmp <- sh.path[[1]][[i]]
		if(length(tmp)>1){
			tmp2 <- tmp[1:(length(tmp)-1)]
			tmp3 <- tmp[2:length(tmp)]
			tmp4 <- rbind(tmp2,tmp3)
			path.edge[,cnt:(cnt+length(tmp2)-1)] <- tmp4
			cnt <- cnt+length(tmp2)
		}
	}
	ret <- unique(t(path.edge))
	ret <- t(apply(ret,1,sort))
	return(ret)
}


v <- 1
path.edge2 <- my.path.edge(g,v,Edge.len)

# あるノードからの最短パス描図
plot3d(V.xyz)
spheres3d(V.xyz[which(tips==0),],radius=0.00001,col=1)
spheres3d(V.xyz[which(tips==1),],radius=0.1,col=3)
spheres3d(V.xyz[v,], radius=0.1,col=2)

segments3d(V.xyz[c(t(path.edge2)),])

# 三角形を構成する３頂点の
# vからの最短距離の差と
# 三角形の辺の長さの関係を調べる
tmp.ret <- matrix(0,length(bunny$it[1,]),3)
Z.ctr <- tmp.ret
cnts <- rep(0,4)
for(i in 1:length(bunny$it[1,])){
	vs <- bunny$it[,i]
	vs <- sort(vs)
	tmp.dist <- c(d[v,vs])
	#print("===")
	#print(c(dist(tmp.dist)))
	#print(dist(tmp.dist))
	#print("---")
	edge.id.1 <- which(Edges[,1] == vs[1] & Edges[,2] == vs[2])
	edge.id.2 <- which(Edges[,1] == vs[1] & Edges[,2] == vs[3])
	edge.id.3 <- which(Edges[,1] == vs[2] & Edges[,2] == vs[3])

	edge.lens <- Edge.len[c(edge.id.1,edge.id.2,edge.id.3)]
	tmp.ret[i,] <- (tmp.dist-edge.lens)
	
	Z.ctr[i,1] <- (d[v,vs[1]] + d[v,vs[2]])/2 - Edge.len[edge.id.1]
	Z.ctr[i,2] <- (d[v,vs[2]] + d[v,vs[3]])/2 - Edge.len[edge.id.2]
	Z.ctr[i,3] <- (d[v,vs[3]] + d[v,vs[1]])/2 - Edge.len[edge.id.3]
	
	path.edge.TF <- rep(NA,3)
	path.edge.id.1 <- which(path.edge2[,1]==vs[1] & path.edge2[,2]==vs[2])
	path.edge.id.2 <- which(path.edge2[,1]==vs[1] & path.edge2[,2]==vs[3])
	path.edge.id.3 <- which(path.edge2[,1]==vs[2] & path.edge2[,2]==vs[3])
	if(length(path.edge.id.1)>0){
		path.edge.TF[1] <- path.edge.id.1
	}
	if(length(path.edge.id.2)>0){
		path.edge.TF[2] <- path.edge.id.2
	}
	if(length(path.edge.id.3)>0){
		path.edge.TF[3] <- path.edge.id.3
	}

	if(length(which(is.na(path.edge.TF)))==0){ # 3辺がパスグラフに含まれる(ないはず)
		cnts[1] <- cnts[1] + 1
	}else if(length(which(is.na(path.edge.TF)))==1){ # 2辺がパスグラフに含まれる
		cnts[2] <- cnts[2] + 1
		# パスに含まれない辺の長さが、結ばれていない頂点間距離となるように
		# １点加えて、Y字にする
	}else if(length(which(is.na(path.edge.TF)))==2){ # 1辺がパスグラフに含まれる
		cnts[3] <- cnts[3] + 1
		# 
	}else{ # 0辺がパスグラフに含まれる
		cnts[4] <- cnts[4] + 1
	}
}
# 実質的に0以下であるかどうかの確認
range(tmp.ret)
head(Z.ctr)
cnts


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
spheres3d(V.xyz,radius=0.05,col=col3)

# Change the vertex, distance form which determines colors

n2 <- 400
col_ <- rgb(d.[n2,],1-d.[n2,],1)
col2_ <- ceiling(d.[n2,] * 15)+1
# color values are 1 or 2
col3_ <- col2_ %% 2 + 1
plot3d(V.xyz)
spheres3d(V.xyz,radius=0.5,col=col3_)


# combining two distance vectors from two vertices and generate color values 1 or 2
# This should make square lattice type colorings
plot3d(V.xyz)
spheres3d(V.xyz,radius=0.05,col=as.numeric(col3 == col3_)+1)

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
spheres3d(V.xyz,radius=0.05,col=col.region2)


plot3d(rbind(V.xyz,V.xyz+max(V.xyz)))
spheres3d(V.xyz[v,]+max(V.xyz),radius=0.01,col=1)
spheres3d(V.xyz[tips.id,]+max(V.xyz),radius=0.01,col=2)
spheres3d(V.xyz[v,],radius=0.01,col=1)
spheres3d(V.xyz[tips.id,],radius=0.01,col=1:length(tips.id)+1)
spheres3d(V.xyz,radius=0.05,col=col.region)


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
spheres3d(V.xyz,radius=0.05,col=col3_)

tips <- my.tips(g,v,d)
col4 <- rep(1,length(d[1,]))
col4[v] <- 2
col4[which(tips==1)] <- 3
plot3d(V.xyz)
spheres3d(V.xyz[which(tips==0),],radius=0.00001,col=1)
spheres3d(V.xyz[which(tips==1),],radius=0.01,col=3)
spheres3d(V.xyz[v,], radius=0.01,col=2)



plot3d(V.xyz,size=0.01)
spheres3d(V.xyz[which(tips==0),],radius=0.0000001,col=gray(0.1))
spheres3d(V.xyz[which(tips==1),],radius=0.01,col=3)
spheres3d(V.xyz[v,], radius=0.01,col=2)

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


