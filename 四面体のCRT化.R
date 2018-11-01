#br.ex <- seq(from=0,to=6.5,by=0.5)
#br.ex <- c(br.ex,6,6.5,7,6.5,6,5.5,5,4.5,4,4.5,5,4.5,4,3.5,3,3.5,4,3.5,3,2.5,2,1.5,1,0.5,0)

br.ex <- c(0,6.5,6,7,4,5,3,4,0)

plot(br.ex,pch=20,type="b")
my.Rtree.dist <- function(g,s,t){
    g[s] + g[t] - 2 * min(g[s:t])
}

my.Rtree.dist.mat <- function(g){
    L <- length(g)
    D.mat <- matrix(0,L,L)
    
    for(i in 1:(L-1)){
        for(j in (i+1):L){
            D.mat[i,j] <- D.mat[j,i] <- my.Rtree.dist(g,i,j)
        }
    }
    return(D.mat)
}
D.mat <- my.Rtree.dist.mat(br.ex)
image(D.mat)

library(ape)
tr <- nj(D.mat)

plot(nj(my.Rtree.dist.mat(br.ex)),type="u",show.tip.label=FALSE)

library(igraph)
my.extree <- function(br.ex,minval = 10^(-15)){
    n <- length(br.ex)
    D.mat <- my.Rtree.dist.mat(br.ex)
    tr <- nj(D.mat)
    el <- tr$edge
    el[which(el==n)] <- 1
    L <- tr$edge.length
    shortL <- (L < minval)
    while(max(el) > n){
        el <- t(apply(el,1,sort))
        tmp <- el - (n+0.5)
        tmp2 <- tmp[,1] * tmp[,2]
        tmp3 <- which(shortL & (tmp2 < 0))
        for(i in 1:length(tmp3)){
            el[which(el==el[tmp3[i],2])] <- el[tmp3[i],1]
        }
    }
    el.true <- which(el[,1]!=el[,2])
    el2 <- el[el.true,]
    L2 <- L[el.true]
    return(list(el=el2,len=L2))
}


extree <- my.extree(br.ex)
g <- graph.edgelist(extree$el,directed=FALSE)
g$weight <- extree$len
plot(g)

my.circuit.extree <- function(br.ex){
    extree <- my.extree(br.ex)
    g <- graph.edgelist(extree$el,directed=FALSE)
    n <- length(br.ex)
    ret <- c(1)
    for(i in 1:(n-2)){
        tmp <- shortest_paths(g,from=i,to=i+1)[[1]][[1]]
        ret <- c(ret,tmp[-1])
    }
    tmp <- shortest_paths(g,from=n-1,to=1)[[1]][[1]]
    ret <- c(ret,tmp[-1])
    ret[length(ret)] <- n
    return(ret)
}

circuit <- my.circuit.extree(br.ex)
plot(br.ex[circuit])

g.circuit <- graph.edgelist(cbind(circuit[1:(length(circuit)-1)],circuit[2:length(circuit)]))
plot(g.circuit) # 第１ノードと最終ノードは別ノードとして扱っているが、同一地点に相当することに注意

circuit

my.extreeZ.2 <- function(extree,edgeZlen){
    #extree <- my.extree(br.ex)
    g <- graph.edgelist(extree$el,directed=FALSE)
    #edgeZlen <- my.edgeWiener(extree$len)
    path.from.1 <- shortest_paths(g,from=1,output="epath")
    Z <- rep(0,length(br.ex)-1)
    for(i in 2:length(Z)){
        tmp <- as_ids(path.from.1[[2]][[i]])
        Z[i] <- sum(edgeZlen[tmp])
    }
    Z <- c(Z,0)
    return(Z)
}

Z <- my.extreeZ.2(extree,c(1,1,-2,-0.5,1,1,3))

my.D.circuit <- function(v,Z,i,j){
    n <- length(v)
    Z. <- Z[v]
    vij <- i:j
    vji <- c(1:n,1:n)[j:(length(v)+i)]
    tmp <- max(min(Z.[vij]),min(c(Z.,Z.)[vji]))
    return(Z.[i]+Z.[j]-2*tmp)
}
# すべての
my.D.circuit.mat <- function(v,Z){
    n <- length(v)
    ret <- matrix(0,n,n)
    for(i in 1:(n-1)){
        for(j in (i+1):n){
            ret[i,j] <- ret[j,i] <- my.D.circuit(v,Z,i,j)
        }
    }
    return(ret)
}

Dcmat <- my.D.circuit.mat(circuit,Z)
image(Dcmat)

# d_e(s,t)=0のIDを取り出す
my.exUpDown <- function(br.ex){
    L <- length(br.ex)
    d <- diff(br.ex)
    #print(d)
    d.sign <- sign(d)
    #print(d.sign)
    d.d.sign <- diff(d.sign)
    dd.positive <- which(d.d.sign>0) + 1
    dd.negative <- which(d.d.sign<0) + 1
    #print(dd.positive)
    #print(dd.negative)
    col <- rep(1,L)
    col[dd.positive] <- 2
    col[dd.negative] <- 3
    #plot(br.ex,col=col,pch=20,type="b")
    return(col)
}
col <- my.exUpDown(br.ex)
plot(br.ex,col=col,pch=20,type="b")

plot(c(Z[circuit],Z[circuit]),type="b",pch=20,col=c(col[circuit],col[circuit]))
library(gtools)
valleyID <- which(my.exUpDown(br.ex)==2)
# d_e(s,t)に相当するcircuitノードペアを作る
my.valley.ids <- function(br.ex){
  valleyID <- which(my.exUpDown(br.ex)==2)
  el <- matrix(0,0,2)
  for(i in 1:length(valleyID)){
      tmp <- which(circuit==valleyID[i])
      tmp.comb <- combinations(length(tmp),2)
      el <- rbind(el,cbind(tmp[tmp.comb[,1]],tmp[tmp.comb[,2]]))
  }
  return(list(id=unique(c(el)),pairs=el))
}
my.valley.ids(br.ex)

my.Dc.circuit <- function(i,j,Dcmat,valleys){
  tmp <- Dcmat[c(i,j,valleys$id),c(i,j,valleys$id)]
  n <- length(tmp[,1])
  pairs <- as.matrix(expand.grid(1:n,1:n))
  g <- graph.edgelist(pairs,directed=FALSE)
  ret <- distances(g,1,2,weights=c(tmp))
  return(ret)
}
my.Dc.circuit.mat <- function(br.ex,Z){
  valleys <- my.valley.ids(br.ex)
  circuit <- my.circuit.extree(br.ex)
  n <- length(circuit)
  ret <- matrix(0,n,n)
  Dcmat <- my.D.circuit.mat(circuit,Z)
  Dcmat2 <- Dcmat
  Dcmat2[(valleys$pairs[,1]-1) * n + valleys$pairs[,2]] <- 0

  for(i in 1:(n-1)){
    for(j in (i+1):n){
      Dcmat3 <- Dcmat2
      Dcmat3[i,] <- Dcmat[i,]
      Dcmat3[j,] <- Dcmat[j,]
      Dcmat3[,i] <- Dcmat[,i]
      Dcmat3[,j] <- Dcmat[,j]
      ret[i,j] <- ret[j,i] <- my.Dc.circuit(i,j,Dcmat3,valleys)
    }
  }
  ret
}
D <- my.Dc.circuit.mat(br.ex,Z)
image(D)

