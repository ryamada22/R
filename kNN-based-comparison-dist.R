#https://stats.stackexchange.com/questions/165047/how-are-graphs-of-k-nearest-neighbors-built-for-clustering 
n.gr <- 2
N <- sample(100:150,n.gr)

p <- 3
q <- 2

X <- matrix(0,0,p)
Y <- matrix(0,0,q)
sh <- sample(1:sum(N))

for(i in 1:n.gr){
	tmpX <- matrix(rnorm(N[i]*p),ncol=p)
	tmpY <- matrix(rnorm(N[i]*q),ncol=q)
	mx <- rnorm(p)*5
	my <- rnorm(q)*5
	tmpX <- t(t(tmpX) + mx)
	tmpY <- t(t(tmpY) + my)
	X <- rbind(X,tmpX)
	Y <- rbind(Y,tmpY)
}

library(igraph)

adj <- matrix(1,sum(N),sum(N))
diag(adj) <- 0

g <- graph.adjacency(adj,mode="undirected")
d.X <- dist(X)
d.Y <- dist(Y)
d.shY <- dist(Y[sh,])

mst.X <- mst(g,c(d.X))
mst.Y <- mst(g,c(d.Y))

mst.XY <- mst(g,c(d.X) + c(d.Y))
mst.XY.sh <- mst(g,c(d.X) + c(d.shY))

ed.mst.X <- get.edgelist(mst.X)
ed.mst.Y <- get.edgelist(mst.Y)
ed.mst.XY <- get.edgelist(mst.XY)
ed.mst.XY.sh <- get.edgelist(mst.XY.sh)

w.X <- rep(0,length(ed.mst.X[,1]))
w.Y <- rep(0,length(ed.mst.Y[,1]))
w.XY <- rep(0,length(ed.mst.XY[,1]))
w.XY.sh <- rep(0,length(ed.mst.XY.sh[,1]))

for(i in 1:length(w.X)){
	w.X[i] <- as.matrix(d.X)[ed.mst.X[i,1],ed.mst.X[i,2]]
}

for(i in 1:length(w.Y)){
	w.Y[i] <- as.matrix(d.Y)[ed.mst.Y[i,1],ed.mst.Y[i,2]]
}

for(i in 1:length(w.XY)){
	w.XY[i] <- (as.matrix(d.X)+as.matrix(d.Y))[ed.mst.XY[i,1],ed.mst.XY[i,2]]
}

for(i in 1:length(w.XY.sh)){
	w.XY.sh[i] <- (as.matrix(d.X)+as.matrix(d.shY))[ed.mst.XY.sh[i,1],ed.mst.XY.sh[i,2]]
}

d.g.X <- distances(mst.X,weight=w.X)
d.g.Y <- distances(mst.Y,weight=w.Y)
d.g.XY <- distances(mst.XY,weight=w.XY)
d.g.XY.sh <- distances(mst.XY.sh,weight=w.XY,sh)

sum(d.g.X^2)
sum(d.g.Y^2)
sum(d.g.XY^2)
sum(d.g.XY.sh^2)

