my.ip <- function(v1,v2,M){
	matrix(v1,nrow=1) %*% M %*% matrix(v2,ncol=1)
}
my.ip.mat <- function(x,M){
	n <- length(x[,1])
	ret <- matrix(0,n,n)
	for(i in 1:n){
		for(j in 1:n){
			ret[i,j] <- my.ip(x[i,],x[j,],M)
		}
	}
	return(ret)
}
library(devtools)
install_github("ryamada22/Ronlyryamada")
library(Ronlyryamada)
library(RFOC)

n <- 5
k <- 5
n.mesh <- 32
A. <- matrix(runif(n^2), n, n)
A.[1, 1] <- k
A. <- A. + rnorm(n^2, 0, 0.05)
xxx <- my.spherical.harm.mesh(A = A., n = n.mesh)

X <- xxx[[1]]
g <- graph.edgelist(xxx$edge,directed=FALSE)
#plot(g)

e.len <- rep(0,length(xxx$edge[,1]))
for(i in 1:length(e.len)){
	e.len[i] <- sqrt(sum((X[xxx$edge[i,1],]-X[xxx$edge[i,2],])^2))
}
D <- distances(g,weights=e.len)

IP <- diag(rep(1,length(X[,1])))
for(i in 1:(length(X[,1])-1)){
	for(j in (i+1):length(X[,1])){
		IP[i,j] <- IP[j,i] <- (D[i,j]^2 - IP[i,i]-IP[j,j])/2
	}
}

#IP

eout <- eigen(IP)

plot(eout[[1]])

s <- sign(eout[[1]])
M <- diag(s)

sq.ev <- sqrt(eout[[1]] * s)
V <- eout[[2]] %*% diag(sq.ev)
plot(sq.ev)

range(IP - my.ip.mat(V,M))

library(rgl)
#plot3d(V[,c(1,n-1,n)])

