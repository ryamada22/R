library(Ronlyryamada)
library(rgl)
library(RFOC)
library(igraph)
library(knitr)
library(tagcloud)
library(e1071)

my.WarshallFloyd <- function(g,w){
  v <- V(g)
  E(g)$weight <- w

  d <- as_adjacency_matrix(g,attr="weight")
  d <- as.matrix(d)
  d[which(d==0)] <- NA
  sh <- allShortestPaths(d)
  return(sh)
}

# 形の凹凸・複雑さをコントロールするパラメタ、n,k
n <- 6
k <- 5
# メッシュのノード数をコントロールするパラメタ
n.mesh <- 16 # 色々試すなら、32くらいにしておくのが無難。送ったhtmlファイルはn.mesh=64
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

sh <- my.WarshallFloyd(g,w)

nv <- length(V(g))

paths <- list()
for(i in 1:nv){
	paths[[i]] <- list()
}
for(i in 1:(nv-1)){
	for(j in (i+1):nv){
		paths[[i]][[j]] <- extractPath(sh,i,j)
		paths[[j]][[i]] <- paths[[i]][[j]][length(paths[[i]][[j]]):1]
	}
}

tips <- list()

for(i in 1:nv){
	tmp <- unlist(paths[[i]])
	tips[[i]] <- which(table(tmp) == 1)
}

cycles <- list()
cycles.half1 <- cycles.half2 <- list()
cnt <- 1
for(i in 1:nv){
	this.tips <- tips[[i]]
	for(j in 1:(length(this.tips)-1)){
		for(k in (j+1):length(this.tips)){
			if(ad[this.tips[j],this.tips[k]] == 1){
				path1 <- paths[[i]][[this.tips[j]]]
				path2 <- paths[[i]][[this.tips[k]]]
				two <- c(path1,path2[length(path2):1])
				un <- unique(two)
				if(length(two) == length(un)+1){
					cycles[[cnt]] <- two[-1]
					cycles.half1[[cnt]] <- path1
					cycles.half2[[cnt]] <- path2
					cnt <- cnt + 1
				}
			}
		}
	}
}

cycle.len <- sapply(cycles,length)

plot3d(xxx$v)
segments3d(xxx$v[c(t(xxx$edge)), ])

for(i in 1:length(cycles)){
	el <- cbind(cycles[[i]],c(cycles[[i]][-1],cycles[[i]][1]))
	segments3d(xxx$v[c(t(el)),],color=i,lw=4)
}

kubis <- list()
cnt <- 1
for(i in 1:length(cycles)){
	len1 <- length(cycles.half1[[i]])
	len2 <- length(cycles.half2[[i]])
	br <- FALSE
	for(j in 2:(len1)){
		if(br){
			break
		}
		for(k in 2:(len2)){
			tmp.path <- extractPath(sh,cycles.half1[[i]][j],cycles.half2[[i]][k])
			if(j==len1 & k==len2){
				br <- FALSE
				break
			}else{
				if(length(tmp.path)==2){
					br <- TRUE
					break
				}
				un <- unique(c(cycles[[i]],tmp.path))
				if(length(un) > length(cycles[[i]])){
					br <- TRUE
					break
				}
			}
			#print(tmp.path)
			#print(un)
		}
	}
	if(!br){
		kubis[[cnt]] <- cycles[[i]]
		cnt <- cnt + 1
	}
}

plot3d(xxx$v)
segments3d(xxx$v[c(t(xxx$edge)), ])

for(i in 1:length(kubis)){
	el <- cbind(kubis[[i]],c(kubis[[i]][-1],kubis[[i]][1]))
	segments3d(xxx$v[c(t(el)),],color=i+1,lw=4)
}


