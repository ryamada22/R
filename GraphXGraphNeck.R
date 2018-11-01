library(rgl) # package for 3d object handling
library(igraph) # package for graph theory
# reading the bunny_200.obj 
bunny <- readOBJ("bunny_200.obj")
# 3D coordinates of vertices in the shape of n x 3 matrix
V.xyz <- t(bunny[[1]][1:3,])
# Enumerate edges of triangle in n x 2 matrix shape
Edges <- rbind(t(bunny[[2]][1:2,]),t(bunny[[2]][2:3,]),t(bunny[[2]][c(3,1),]))
# Remove duplicates of edges
Edges <- unique(Edges)
# length of edges
Edge.len <- sqrt(apply((V.xyz[Edges[,1],] - V.xyz[Edges[,2],])^2,1,sum))
# make a graph object
g <- graph.edgelist(Edges,directed =FALSE)
# distance on the graph
d <- distances(g,weights=Edge.len)
# 隣接ノード
adj.v <- adjacent_vertices(g,1:length(V(g)))
# path
#path.list <- list()
#for(i in 1:length(V.xyz[,1])){
#    path.list[[i]] <- all_shortest_paths(g,i,mode="out")
#}
# Triangle list
tris <- list()
for(i in 1:length(V.xyz[,1])){
	tris[[i]] <- which(apply(bunny[[2]]-i,2,prod)==0)
}

angle.v <- rep(0,length(tris))

for(i in 1:length(angle.v)){
	for(j in 1:length(tris[[i]])){
		ids <- bunny[[2]][,tris[[i]][j]]
		this <- which(ids==i)
		others <- which(ids!=i)
		v1 <- V.xyz[others[1],] - V.xyz[this,]
		v2 <- V.xyz[others[2],] - V.xyz[this,]
		cos <- sum(v1*v2)/(sqrt(sum(v1^2))*sqrt(sum(v2^2)))
		angle.v[i] <- angle.v[i] + acos(cos)
	}
}

plot(sort(angle.v))


my.cycle <- function(g){
    n <- length(V(g))
    Edges <- get.edgelist(g)
    ret <- matrix(0,0,3)
    for(i in 1:length(Edges[,1])){
        v1 <- Edges[i,1]
        v2 <- Edges[i,2]
        for(j in 1:n){
            p1 <- shortest_paths(g,v1,j,mode="out")
            p2 <- shortest_paths(g,v2,j,mode="out")
            p12 <- c(unlist(p1[[1]]),unlist(p2[[1]]))
            p12. <- unique(p12)
            if(length(p12) == length(p12.)+1){
                ret <- rbind(ret,c(v1,v2,j))
            }
        }
    }
    ret
}

# 極遠点検出関数
my.Far <- function(d,adj.v){
    n <- length(adj.v)
    ret <- matrix(0,n,n)
    for(i in 1:n){
        tmp <- d[,adj.v[[i]]]
        tmp2 <- apply(tmp < d[,i],1,prod)
        ret[which(tmp2==1),i] <- 1
    }
    return(ret)
}

#lapply(adj.v,length)
fars <- my.Far(d,adj.v)
fars.face <- my.Far(d.face,adj.v.face)

image(fars)
hist(apply(fars,1,sum))
plot(sort(apply(fars,2,sum)))

my.minPair <- function(d){
    n <- length(d[,1])
    m <- length(d[1,])
    tmp <- t(apply(d,1,order))
    tmp2 <- apply(d,2,order)
    tmp.pair <- rbind(cbind(1:n,tmp[,1]),cbind(tmp2[1,],1:m))
    #tmp.pair <- t(apply(tmp.pair,1,sort))
    tmp.pair <- unique(tmp.pair)
    return(tmp.pair)
}
tmp.d <- matrix(c(1,2,3,4,5,6),ncol=3)
my.minPair(tmp.d)

# el: 元のグラフのエッジリスト
# x: 元のグラフの頂点座標
# fars: 極遠点情報行列

my.dist.pair <- function(X1,X2){
    L1 <- apply(X1^2,1,sum)
    L2 <- apply(X2^2,1,sum)
    IP <- X1 %*% t(X2)
    tmp <- outer(L1,L2,"+") - 2 * IP
    tmp[which(tmp < 0)] <- 0
    sqrt(tmp)
}

my.farManifold <- function(fars,el,x){
    # 頂点ペア(新たな多様体グラフの頂点になる)
    V <- which(fars==1,arr.ind=TRUE)
    n.el <- length(el[,1])
    ret <- matrix(0,0,4)
    ncol <- length(x[1,])
    for(i in 1:n.el){
        f1 <- which(fars[el[i,1],]==1)
        f2 <- which(fars[el[i,2],]==1)
        #F1 <- cbind(rep(el[i,1],length(f1)),f1)
        #F2 <- cbind(rep(el[i,2],length(f2)),f2)
        #X1 <- cbind(x[F1[,1],],x[F1[,2],])
        #X2 <- cbind(x[F2[,1],],x[F2[,2],])
        X1 <- matrix(x[f1,],ncol=ncol)
        X2 <- matrix(x[f2,],ncol=ncol)
        
        tmp <- my.dist.pair(X1,X2)
        tmp2 <- my.minPair(tmp)
        n.pair <- length(tmp2[,1])
        ret <- rbind(ret,cbind(rep(el[i,1],n.pair),tmp2[,1],rep(el[i,2],n.pair),tmp2[,2]))
    }
    return(list(V=V,E=ret))
}

el <- Edges
x <- V.xyz
out <- my.farManifold(fars,el,x)

# V は(点,極遠点)ペア
# E はV のペア((点,極遠点)ペア,(点,極遠点)ペア)
my.graphM2 <- function(V,E){
    nv <- max(V)
    new.V.id <- (V[,1]-1) * nv + V[,2]
    new.E.id <- cbind(paste("",(E[,1]-1) * nv + E[,2]), paste("",(E[,3]-1) * nv + E[,4]))
    new.g <- graph.edgelist(new.E.id,directed=FALSE)
    return(list(g=new.g,vid = new.V.id,eid=new.E.id))
}

graphM2 <- my.graphM2(out$V,out$E)

# V は(点,極遠点)ペア
# x はオリジナルグラフのノード座標
my.coordsM2 <- function(V,x){
    cbind(x[V[,1],],x[V[,2],])
}

is.connected(graphM2$g)

plot(graphM2$g,vertex.size=0.1,vertex.label=NA)

