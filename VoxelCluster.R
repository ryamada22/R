ns <- 100

na <- 50
nt <- 10
X <- matrix(0,ns,ns)

for(i in 1:na){
	tmp <- matrix(sample((-1):1,2*(nt-1),replace=TRUE,prob=c(0.2,0.6,0.2)),nt-1,2)
	tmp <- rbind(sample(ns,2),tmp)
	tmp2 <- apply(tmp,2,cumsum)
	for(j in 1:nt){
		rg <- range(tmp2[j,])
		if(rg[1]>=1 & rg[2] <=ns){
			X[tmp2[j,1],tmp2[j,2]] <- 1 #runif(1)*0.5+0.5
		}
		
	}
}

image(X)

X.add <- which(X> min(X)-1,arr.ind=TRUE)

XV <- cbind(X.add,c(X)*ns)

library(rgl)

plot3d(XV)

clout <- hclust(dist(XV))

plot(clout)

X.1 <- which(X>0)
X.add.1 <- which(X>0,arr.ind=TRUE)
XV1 <- cbind(X.add.1,X[X.1])

plot3d(XV1)

clout1 <- hclust(dist(XV1))
plot(clout1)

library(igraph)

dmat <- as.matrix(dist(XV1))

dmat.add <- which(dmat>(min(dmat)-1),arr.ind=TRUE)

ord <- order(dmat)

dmat.add.sort <- dmat.add[ord,]

checker <- rep(-1,max(dmat.add))
cnt <- 1
loop <- TRUE
while(loop){
	checker[dmat.add.sort[i,]] <- checker[dmat.add.sort[i,]] + 1
	cnt <- cnt+1
	print(cnt)
	if(sum(checker)==max(dmat.add)){
		loop <- FALSE
	}
}

plot(X.add.1,pch=20,cex=0.1)
for(i in 1:cnt){
	segments(X.add.1[dmat.add.sort[i,1],1],X.add.1[dmat.add.sort[i,1],2],X.add.1[dmat.add.sort[i,2],1],X.add.1[dmat.add.sort[i,2],2])
}

g <- graph.adjacency(as.matrix(dist(XV1)))



m.s.t <- mst(g)

plot(m.s.t)

