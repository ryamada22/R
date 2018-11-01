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

# グラフが平面グラフであるとする
# 球面と同相なので、単位球面に貼り付けることを考える
# 頂点間距離から、各点の対蹠ノードを決める
# 対蹠ノードまでの距離は、単位球面に張り付いたとするとpiである
# 対蹠ノード以外までの距離は、次のようにして定める
# Viの対蹠ノードがVjであるとき、ViからVkまでの距離は
# (Vi-Vkグラフ距離)/(Vi-Vkグラフ距離＋Vk-Vjグラフ距離) * pi
# このようにすると、Vi-Vk距離とVk-Vi距離とが異なるので
# 距離に対称性を持たせるために、両者の算術平均をVi-Vk間距離とする
# これにより全ノードのペアワイズ対象距離行列が得られるが
# この距離は単位球面上の測地線に擬せられているから
# ２点を単位ベクトルと見たときの、ベクトルがなす角度とも看做せる
# これを用いて、２つの単位ベクトルの内積(cos(theta))も計算できる
# 内積行列を固有値分解して、各ノードの座標を出し
# ３次元分単位球面部分だけを取り出し、残差は形のゆがみによるものとして捨てることにする


my.spherization <- function(g,w){
	d <- distances(g,weights=w)
	nv <- length(d[,1])
	taiseki <- apply(d,1,which.max)
	k <- d
	for(i in 1:nv){
		k[i,] <- d[i,]/(d[i,] + d[taiseki[i],])
	}
	k. <- (k + t(k))/2 * pi
	k.. <- cos(k.)
	out <- eigen(k..)
	x <- out[[2]][,1:3]
	x <- x/sqrt(apply(x^2,1,sum))
	return(list(x=x,dgr=d,dst.un=k,dst.sym=k.,dcos=k..,eigen.out = out))
}

my.spherization.d <- function(g,w){
	d <- distances(g,weights=w)
	nv <- length(d[,1])
	taiseki <- apply(d,1,which.max)
	k <- d
	for(i in 1:nv){
		k[i,] <- d[i,]/(d[i,] + d[taiseki[i],])
	}
	k. <- (k + t(k))/2 * pi
	k.. <- cos(k.)
	#out <- eigen(k..)
	#x <- out[[2]][,1:3]
	#x <- x/sqrt(apply(x^2,1,sum))
	return(list(dgr=d,dst.un=k,dst.sym=k.,dcos=k..,taiseki=taiseki))
}



library(RFOC)
library(igraph)
library(knitr)
library(tagcloud)
library(e1071)

# 形の凹凸・複雑さをコントロールするパラメタ、n,k
n <- 6
k <- 3
# メッシュのノード数をコントロールするパラメタ
n.mesh <- 20 # 色々試すなら、32くらいにしておくのが無難。送ったhtmlファイルはn.mesh=64
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
s <- sample(1:nv,round(nv*0.05))
typed[s] <- sample(1:7,length(s),replace=TRUE)

# Voronoi label
new.t <- my.graph.voronoi(g,w,typed)

# color only typed nodes
plot3d(xxx$v)
spheres3d(xxx$v,color=typed+1,radius=0.05)
segments3d(xxx$v[c(t(xxx$edge)), ])


# my.spherization
spout <- my.spherization(g,w)
x <- spout$x

# color all voronoi labeled nodes
#plot3d(xxx$v)
#spheres3d(xxx$v,color=typed[new.t]+1,radius=0.1)
#segments3d(xxx$v[c(t(xxx$edge)), ])

#open3d()
#plot3d(x)
#spheres3d(x,color=typed[new.t]+1,radius=0.1)
#segments3d(x[c(t(xxx$edge)), ])


# たとえば、第200番ノードからの測地線のうち、構成エッジ数が最大のものをプロットしてみる
# オリジナルの形表面ではのたくるが
# 強制球面上では、大円の半分に近くなっている
sh <- shortest_paths(g,200,weights=w)

len <- sapply(sh[[1]],length)

plot3d(xxx$v)
spheres3d(xxx$v[sh[[1]][[which.max(len)]],],color = "blue",radius=0.1)
segments3d(xxx$v[c(t(xxx$edge)), ])

open3d()
plot3d(x)
segments3d(x[c(t(xxx$edge)), ])
spheres3d(x[sh[[1]][[which.max(len)]],],color = i,radius=0.1)
# 各ノードごとに描図してみる
#for(i in 1:nv){
#	sh <- shortest_paths(g,i,weights=w)
#	len <- sapply(sh[[1]],length)
#	plot3d(x)
#	segments3d(x[c(t(xxx$edge)), ])
#	spheres3d(x[sh[[1]][[which.max(len)]],],color = i,radius=0.1)
#}

# 対称化測地距離で色塗りをすると、まずまずの同心円が現れる

plot3d(x)
segments3d(x[c(t(xxx$edge)), ])
col <- round(k.[200,]/max(k.[200,]) * 6) + 1
spheres3d(x,radius=0.07,color=col)

# オリジナルのグラフ距離でも、まずまずの同心円が現れる
plot3d(x)
segments3d(x[c(t(xxx$edge)), ])
col <- round(d[200,]/max(d[200,]) * 6) + 1
spheres3d(x,radius=0.07,color=col)

#for(i in 1:nv){
#	plot3d(x)
#	segments3d(x[c(t(xxx$edge)), ])	
#	#col <- round(k.[i,]/max(k.[i,]) * 6) + 1
#	col <- round(d[i,]/max(d[i,]) * 6) + 1
#	spheres3d(x,radius=0.07,color=col)
#	for(i in 1:10^3){
#		col <- sample(col)
#	}
#}
