# g: igraph object
# w: edge length
# t: type of nodes; 0 means non-specified
# label all 0-labeled vertices with typed vertex id
my.graph.voronoi <- function(g,w,t){
	typed <- which(t != 0)
	nontyped <- which(t == 0)
	d <-distances(g,weights=w)
	new.t <- 1:length(t)
	for(i in 1:length(nontyped)){
		v <- nontyped[i]
		d.typed <- d[v,typed]
		nearest <- typed[which(d.typed==min(d.typed))]
		new.t[v] <- nearest
	}
	return(new.t)
}

# devtools::install_github("ryamada22/Ronlyryamada")
library(Ronlyryamada)
library(rgl)
library(RFOC)
library(igraph)
library(knitr)
library(tagcloud)
library(e1071)

# 形の凹凸・複雑さをコントロールするパラメタ、n,k
n <- 6
k <- 5
# メッシュのノード数をコントロールするパラメタ
n.mesh <- 64 # 色々試すなら、32くらいにしておくのが無難。送ったhtmlファイルはn.mesh=64
# 形を球面調和関数係数ベクトルで指定する
A. <- matrix(runif(n^2), n, n)
A.[1, 1] <- k
B <- matrix(rnorm(n^2), n, n)
# 閉曲面オブジェクトを作る
xxx <- my.spherical.harm.mesh(A = A., B = B, n = n.mesh)
plot3d(xxx$v)
segments3d(xxx$v[c(t(xxx$edge)), ])

g <- graph.edgelist(xxx$edge,directed=FALSE)
# edge lengths
w <- sqrt(apply((xxx$v[xxx$edge[,1],]-xxx$v[xxx$edge[,2],])^2,1,sum))

ad <- get.adjacency(g)

nv <- length(V(g))

plot3d(xxx$v)
segments3d(xxx$v[c(t(xxx$edge)), ])


# Types some nodes
typed <- rep(0,nv)
s <- sample(1:nv,round(nv*0.1))
typed[s] <- sample(1:7,length(s),replace=TRUE)

# Voronoi label
new.t <- my.graph.voronoi(g,w,typed)

# color only typed nodes
plot3d(xxx$v)
spheres3d(xxx$v,color=typed+1,radius=0.05)
segments3d(xxx$v[c(t(xxx$edge)), ])

# color all voronoi labeled nodes
plot3d(xxx$v)
spheres3d(xxx$v,color=typed[new.t]+1,radius=0.05)
segments3d(xxx$v[c(t(xxx$edge)), ])


