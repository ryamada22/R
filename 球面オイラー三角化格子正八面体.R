

my.partition.tris <- function(tri.x){
	tris <- tri.x$tris
	X <- tri.x$X
	col <- tri.x$col
	tris <- t(apply(tris,1,sort))
	pairs <- rbind(tris[,1:2],tris[,2:3],tris[,c(1,3)])
	pairs2 <- t(apply(pairs,1,sort))
	pairs2 <- unique(pairs)
	new.v <- (1:length(pairs2[,1])) + max(pairs2)
	new.X <- (X[pairs2[,1],] + X[pairs2[,2],] ) /2
	ret.X <- rbind(X,new.X)
	new.tris <- matrix(0,0,3)
	new.col <- c()
	for(i in 1:length(tris[,1])){
		tri.now <- tris[i,]
		tri.add1 <- which(pairs2[,1]==tri.now[1] & pairs2[,2]==tri.now[2])+ max(pairs2)
		tri.add2 <- which(pairs2[,1]==tri.now[2] & pairs2[,2]==tri.now[3])+ max(pairs2)
		tri.add3 <- which(pairs2[,1]==tri.now[1] & pairs2[,2]==tri.now[3])+ max(pairs2)
		
		tmp <- rbind(c(tri.now[1],tri.add1,tri.add3),c(tri.now[2],tri.add1,tri.add2),c(tri.now[3],tri.add2,tri.add3),c(tri.add1,tri.add2,tri.add3))
		new.tris <- rbind(new.tris,tmp)
		if(col[i]==1){
			tmp.col <- c(2,2,2,1)
		}else{
			tmp.col <- c(1,1,1,2)
		}
		new.col <- c(new.col,tmp.col)
	}
	new.tris <- t(apply(new.tris,1,sort))
	
	ret.X <- ret.X/sqrt(apply(ret.X^2,1,sum))
	return(list(tris=new.tris,X=ret.X,col=new.col))
}

my.plot.eultri <- function(tri){
	plot3d(tri$X,cex=0.01)
	n.Face <- length(tri$tris[,1])
	col <- rep(2:3,n.Face/2)
	col <- c(rbind(tri$col,tri$col,tri$col))
	tmp.ord <- c(t(tri$tris))
	triangles3d(tri$X[tmp.ord,],col=col+1)
}

tris8=matrix(c(1,2,3,1,3,4,1,4,5,1,5,2,2,3,6,3,4,6,4,5,6,5,2,6),byrow=TRUE,ncol=3)
X8 = matrix(c(0,0,1,1,0,0,0,1,0,-1,0,0,0,-1,0,0,0,-1),byrow=TRUE,ncol=3)
col8 = c(1,2,1,2,2,1,2,1)
init.tri.x <- list(tris=tris8,X=X8,col=col8)

n.step <- 4
tri.x.list <- list()
tri.x.list[[1]] <- init.tri.x
for(i in 2:n.step){
	tri.x.list[[i]] <- my.partition.tris(tri.x.list[[i-1]])
}

my.plot.eultri(tri.x.list[[n.step]])

