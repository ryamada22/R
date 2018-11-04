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
# 1次元グラフ
library(igraph)
n <- 1000
el <- cbind(1:(n-1),2:n)
el <- rbind(el,c(n,1)) # ぐるり
g <- graph.edgelist(el,directed=FALSE)

# edgeに長さを規則的に入れる(何かしら意味のある１次元多様体のつもり)
# e.len <- sin((1:(n-1))*0.1) + 2
e.len <- runif(length(el[,1]))
e.len <- rep(1,length(el[,1]))

D <- distances(g,weights=e.len)

# 次を仮定する
# 各ノードは、原点から、1,2,...,nの距離にあるものとする
# ||x-y||^2 = ||x||^2+||y||^2 - 2(x,y) から
# (x,y)の値を定めることができる

IP <- diag(rep(1,n))
for(i in 1:(n-1)){
	for(j in (i+1):n){
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

#range(IP - my.ip.mat(V,M))

library(rgl)
plot3d(V[,c(1,n-1,n)])

