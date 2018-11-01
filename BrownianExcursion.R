Bt <- rnorm(500, 0, 1)
Bt <- cumsum(Bt)
Bt <- rnorm(500, 0, 1)
Bt <- Bt[Bt > 0 & Bt < 1]
Bt <- cumsum(Bt)
Bt <- rnorm(500, 0, 1/sqrt(500))
  Bt <- cumsum(Bt)
  Bt <- ts(Bt, start = 1/500, frequency = 500)
  ts(abs(Bt - time(Bt) * as.vector(Bt)[500]), start = 1/500, 
     frequency = 500)
     
library(e1071)
succeeded <- FALSE
while(!succeeded)
{
  bridge <- rbridge(end = 1, frequency = 500)
  succeeded=all(bridge>=0)
}
plot(bridge)


while(!succeeded)
{
  bridge <- rbridge(end = 1, frequency = 500)
  succeeded=all(bridge>=0)||all(bridge<=0)
}
bridge = abs(bridge)
plot(bridge)

# library(e1071)
my.rexcursion <- function(frequency = 1000, end = 1){
	succeeded <- FALSE
	while(!succeeded)
	{
		bridge <- rbridge(end = end, frequency = frequency)
		succeeded=all(bridge>=0)||all(bridge<=0)
	}
	return(c(0,c(abs(bridge))))
}
my.rwiener <- function(frequency = 1000, end = 1){
	c(0,cumsum(rnorm(end * frequency)/sqrt(frequency)))
}
br.ex <- my.rexcursion(frequency=100)
plot(br.ex)
length(br.ex)


fr <- 10
t <- (1:fr)/fr
br.ex <- my.rexcursion(frequency=fr)
br.z <- my.rwiener(frequency=fr)

library(rgl)
library(ape)
plot3d(t,br.ex,br.z,type="l")

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

tr <- nj(D.mat)

plot(nj(my.Rtree.dist.mat(br.ex)),type="u",show.tip.label=FALSE)

# ノードはbr.exにあるものだけなので、それに集約する。

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
		#print(el)
		#print(max(el))
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
	#ret <- c(ret,n)
	return(ret)
}

circuit <- my.circuit.extree(br.ex)

my.edgeWiener <- function(lens,frequency=1000){
	ret <- rep(0,length(lens))
	for(i in 1:length(ret)){
		tmp <- rwiener(lens[i],frequency=frequency/lens[i])
		ret[i] <- tmp[length(tmp)]
	}
	return(ret)
}
my.edgeWiener(extree$len)

my.extreeZ <- function(br.ex){
	extree <- my.extree(br.ex)
	g <- graph.edgelist(extree$el,directed=FALSE)
	edgeZlen <- my.edgeWiener(extree$len)
	path.from.1 <- shortest_paths(g,from=1,output="epath")
	Z <- rep(0,length(br.ex)-1)
	for(i in 2:length(Z)){
		tmp <- as_ids(path.from.1[[2]][[i]])
		Z[i] <- sum(edgeZlen[tmp])
	}
	Z <- c(Z,0)
	return(Z)
}

Z <- my.extreeZ(br.ex)

my.D.circuit <- function(v,Z,i,j){
	n <- length(v)
	Z. <- Z[v]
	vij <- i:j
	vji <- c(1:n,1:n)[j:(length(v)+i)]
	tmp <- max(min(Z.[vij]),min(c(Z.,Z.)[vji]))
	return(Z.[i]+Z.[j]-2*tmp)
}
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

my.D.circuit.cycle <- function(v,Z){
	n <- length(v)
	ret <- rep(0,n-1)
	for(i in 1:(n-1)){
		ret[i] <- my.D.circuit(v,Z,i,i+1)
	}
	#ret[n] <- my.D.circuit(v,Z,n,1)
	return(ret)
}

Dccycle <- my.D.circuit.cycle(circuit,Z)

valleyID <- which(my.exUpDown(br.ex)==2)

circuit.el <- cbind(1:(length(circuit)-1),2:length(circuit))
circuit.el.len <- Dccycle
circuit.el <- rbind(circuit.el,c(1,length(circuit)))
circuit.el.len <- c(circuit.el.len,0)

library(gtools)

for(i in 1:length(valleyID)){
	tmp <- which(circuit==valleyID[i])
	tmp.comb <- combinations(length(tmp),2)
	circuit.el <- rbind(circuit.el,cbind(tmp[tmp.comb[,1]],tmp[tmp.comb[,2]]))
	circuit.el.len <- c(circuit.el.len,rep(0,length(tmp.comb[,1])))
}
circuit.el
circuit.el.len

g2 <- graph.edgelist(circuit.el,directed=FALSE)
DD2 <- distances(g2,weights=circuit.el.len)

#tr2 <- nj(DD)


pre.circuit <- c(length(br.ex),circuit[-length(circuit)])
post.circuit <- c(circuit[-1],1)

g.cycle.a <- cbind(circuit,post.circuit)
g.cycle.a <- rbind(g.cycle.a,cbind(valleyID,valleyID))

g.cycle.b <- cbind(circuit,pre.circuit)
g.cycle.b <- rbind(g.cycle.b,cbind(valleyID,valleyID))

el <- get.edgelist(g)
el.weight <- rep(0,length(el[,1]))

for(i in 1:length(el.weight)){
	el.weight[i] <- Dcmat[el[i,1],el[i,2]]
}
el.weight

shortest_paths(g,from=1,weights=el.weight)

my.rwiener.tr <- function(tr,frequency=1000){
	len <- tr$edge.length
	v <- rep(0,length(len))
	#len[which(len<0)] <- 0
	for(i in 1:length(len)){
		if(len[i] > 0){
			tmp <- rwiener(len[i],frequency=frequency/len[i])
			v[i] <- c(tmp[length(tmp)])
		}
	}
	v
}
v.len <- my.rwiener.tr(tr)


my.Z.tr <- function(tr,frequency=1000){
	v.len <- my.rwiener.tr(tr)
	
}
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

my.edgeUpDown <- function(br.ex){
	col <- my.exUpDown(br.ex)
	Ups <- Downs <- list()
	dir <- 1
	Ups[[1]] <- c(1)
	for(i in 2:length(col)){
		if(col[i]==1){
			if(dir == 1){
				Ups[[length(Ups)]] <- c(Ups[[length(Ups)]],i)
			}else{
				Downs[[length(Downs)]] <- c(Downs[[length(Downs)]],i)
			}
		}else if(col[i]==2){
			dir <- dir * (-1)
			Ups[[length(Ups)+1]] <- c(i)
			Downs[[length(Downs)]] <- c(Downs[[length(Downs)]],i)
		}else{
			dir <- dir * (-1)
			Downs[[length(Downs)+1]] <- c(i)
			Ups[[length(Ups)]] <- c(Ups[[length(Ups)]],i)
		}
	}
	return(list(Ups=Ups,Downs=Downs))
}

my.edgeUpDown(br.ex)

my.make.circle <- function(br.ex){
	eUD <- my.edgeUpDown(br.ex)
	Ups <- eUD$Ups
	Downs <-eUD$Downs
	currentUp <- Ups[[1]]
	
	n <- length(br.ex)
	circ <- 1:n
	
	for(i in 1:length(Downs)){
		UpsInDown <- which(currentUp < Downs[[i]][1] & currentUp > Downs[[i]][length(Downs[[i]]))
		
		DownsInUp <- which(Downs[[i]] > currentUp[1] & Downs[[i]] < currentUp[length(currentUp))
		
		down.st <- which(circ == Downs[[i]][1])
		down.end <- which(cir == Downs[[i]][length(Downs[[i]]))
		circ <- c(circ[1:down.st],order(br.ex[UpsInDown],decreasing=TRUE))
		
	}
}

my.make.circle0 <- function(br.ex){
	eUD <- my.edgeUpDown(br.ex)

	#circ <- currentUp
	for(i in 1:length(Downs)){
		UpsInDown <- which(currentUp < Downs[[i]][1] & currentUp > Downs[[i]][length(Downs[[i]]))
		
		DownsInUp <- which(Downs[[i]] > currentUp[1] & Downs[[i]] < currentUp[length(currentUp))
		
		ids <- unique(c(currentUp, Downs[[i]]))
		ord <- order(br.ex[ids])
		down.st <- which(ord==Downs[[i]][1])
		down.end <- which(ord==Downs[[i]][length(Downs[[i]]))
		tmp <- ord[(down.st-1):down.end]
		#circ <- c(circ,tmp)
		currentUp <- ord[1:down.end]
		
		if(i < length(eUD$Downs)){
			circ <- c(circ,eUD$Ups[[i+1]][-1])
			currentUp <- unique(c(currentUp,eUD$Ups[[i+1]]))
		}
	}
}
#my.ex2tr(br.ex)
