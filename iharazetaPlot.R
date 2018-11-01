library(igraph)
library(complexplus)
my.Ihara.zeta.elem <- function(g){
	A <- as.matrix(as_adjacency_matrix(g))
	n.v <- length(V(g))
	n.e <- length(E(g))
	r <- n.e - n.v
	D <- diag(degree(g))
	H <- diag(degree(g)-1) 
	I <- diag(rep(1,n.v))
	return(list(r=r,A=A,H=H,I=I,D=D,n.v=n.v,n.e=n.e))
}
my.Ihara.zeta.ori <- function(g,u){
	elem <- my.Ihara.zeta.elem(g)
	1/((1-u^2)^elem$r * Det(elem$I-u*elem$A + u^2*elem$H))
}

my.Bartholdi.zeta <- function(g,u,t){
	elem <- my.Ihara.zeta.elem(g)
	tmp <- (1-(1-u)^2*t^2)^elem$r * Det(elem$I-t*elem$A+(1-u)*(elem$D-(1-u)*elem$I)*t^2)
	return(1/tmp)
}

my.Ihara.zeta <- function(g,u){
	my.Bartholdi.zeta(g,0,u)
}
# (エッジ長ｘ２）^2 行列でのゼータ関数
my.Ihara.zeta.e <- function(g,u){
	el <- as_edgelist(g)
	el2 <- rbind(el,cbind(el[,2],el[,1]))
	n.e <- length(el[,1])
	edge.mat <- matrix(0,n.e*2,n.e*2)
	for(i in 1:(n.e*2)){
		st <- el2[i,1]
		ed <- el2[i,2]
		tmp <- which(el2[,1]==ed & el2[,2]!=st)
		edge.mat[i,tmp] <- 1
	}
	I <- diag(rep(1,n.e*2))
	tmpmat <- I - edge.mat * u
	return(1/Det(tmpmat))
}
# エッジにウェイトがある場合
# Wマトリックスの要素がエッジウェイトによって値を変える
# エッジペアがつながっていたら、２つのエッジの長さの平均を重みにとる
my.Ihara.zeta.e.w <- function(g,u){
	w <- E(g)$weight
	w <- c(w,w)
	if(is.null(w)){
		w <- rep(1,length(E(g))*2)
	}
	el <- as_edgelist(g)
	el2 <- rbind(el,cbind(el[,2],el[,1]))
	n.e <- length(el[,1])
	edge.mat <- matrix(0,n.e*2,n.e*2)
	for(i in 1:(n.e*2)){
		st <- el2[i,1]
		ed <- el2[i,2]
		tmp <- which(el2[,1]==ed & el2[,2]!=st)
		edge.mat[i,tmp] <- u^((w[i]+w[tmp])/2) # u^w[i]にしてもよい？？？
	}
	I <- diag(rep(1,n.e*2))
	tmpmat <- I - edge.mat
	return(1/Det(tmpmat))
}
my.Ihara.zeta.e.w2 <- function(g,u){
	w <- E(g)$weight
	w <- c(w,w)
	if(is.null(w)){
		w <- rep(1,length(E(g))*2)
	}
	el <- as_edgelist(g)
	el2 <- rbind(el,cbind(el[,2],el[,1]))
	n.e <- length(el[,1])
	edge.mat <- matrix(0,n.e*2,n.e*2)
	for(i in 1:(n.e*2)){
		st <- el2[i,1]
		ed <- el2[i,2]
		tmp <- which(el2[,1]==ed & el2[,2]!=st)
		edge.mat[i,tmp] <- u^(w[i]) # u^w[i]にしてもよい->my.Ihara.zeta.e.w2()
	}
	I <- diag(rep(1,n.e*2))
	tmpmat <- I - edge.mat
	return(1/Det(tmpmat))
}


my.Ihara.zeta.path <- function(g,u){
	# まずはedge zetaと同じ処理
	# その後で全域木(minimum spaning tree)を取り出して、その補を取る
	w <- E(g)$weight
	if(is.null(w)){
		w <- rep(1,length(E(g)))
	}
	mst <- mst(g,weights=w)
	
	comst <- g- mst
	n.e.co <- length(E(comst))
	Z <- matrix(0,n.e.co*2,n.e.co*2)
	w.co <- E(comst)$weight
	w.co <- c(w.co,w.co)
	#diag(Z) <- u^w.co
	el <- as_edgelist(comst)
	el2 <- rbind(el,cbind(el[,2],el[,1]))

	for(i in 1:((n.e.co*2))){
		for(j in 1:(n.e.co*2)){
			if(abs(i-j) == n.e.co){
				Z[i,j] <- 0
			#}else if(i==j){
				
			}else{
				a <- el2[i,2]
				b <- el2[j,1]
				tmp <- w.co[i]/2 + w.co[j]/2 + shortest.paths(mst,a,b)
				#tmp <- shortest.paths(mst,a,b,weights=E(mst)$weight)
				Z[i,j] <- u^tmp
				#a <- el2[j,2]
				#b <- el2[i,1]
				#tmp <- w.co[i] + w.co[j] + shortest.paths(mst,a,b)
				#Z[j,i] <- u^tmp
			}
		}
	}
	I <- diag(rep(1,n.e.co*2))
	tmpmat <- I - Z
	return(1/Det(tmpmat))
}

my.Ihara.zeta.path(g,xy.[30])
my.Ihara.zeta.e.w(g,xy.[30])

# エッジに整数長さがあるときに、すべてのエッジ長が１となるように
# ノードを加えたグラフに変換する
# E(g)$weight <- sample(1:20,n.e)

my.graph.inflation.weight <- function(g){
	w <- E(g)$weight
	el <- as_edgelist(g)
	A <- as.matrix(as_adjacency_matrix(g))
	for(i in 1:length(w)){
		if(w[i] > 1){
			nv <- w[i]-1
			st <- el[i,1]
			ed <- el[i,2]
			A[st,ed] <- A[ed,st] <- 0
			n <- length(A[,1])
			newA <- matrix(0,n+nv,n+nv)
			newA[1:n,1:n] <- A
			newA[n+1,st] <- newA[st,n+1] <- 1
			newA[n+nv,ed] <- newA[ed,n+nv] <- 1
			if(nv > 1){
				for(j in 1:(nv-1)){
					newA[n+j,n+j+1] <- newA[n+j+1,n+j] <-1
				}
			}
			A <- newA
		}
	}
	newg <- graph.adjacency(A,mode="undirected")
	return(newg)
}

E(g)$weight <- rep(3,n.e)
newg <- my.graph.inflation.weight(g)


fu <- function(u){
	tmp <- (-1)*(3*u-1)*(u+1)^5*(u-1)^6*(3*u^2+u+1)^4
	return(1/tmp)
}

#g <- sample_gnp(10, 2/10)
n <- 5
adj <- matrix(1,n,n)
diag(adj) <- 0
g <- graph.adjacency(adj,mode="undirected")

x <- y <- seq(from=-0.8,to=0.8,length=100)

xy <- expand.grid(x,y)
xy. <- xy[,1] + 1i*xy[,2]

vs <- rep(0,length(xy.))
vs2 <- vs
for(i in 1:length(vs)){
	vs[i] <- my.Ihara.zeta(g,xy.[i])
	vs2[i] <- fu(xy.[i])
}
my.Ihara.zeta(g,vs[55])
fu(vs[55])
plot(xy,col=abs(Mod(vs)))


fu2 <- function(u){
	tmp <- (1-u^10)^5*(1-3*u^5)*(1-u^5)*(1+u^5+3*u^10)
	return(1/tmp)
}


fu(vs[515]^5)
fu2(vs[515])

