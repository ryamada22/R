# 三角形を(x,y,z)分割する
library(partitions)
library(combinat)
library(rgl)

my.divide.n <- function(n){
	if(n==1){
		return(list(composition=diag(rep(1,3)),tris=matrix(c(1,2,3),nrow=3)))
	}
	comp <- compositions(n,3)
	nv <- length(comp[1,])
	trios <- matrix(combn(1:nv,3),nrow=3)
	e1 <- comp[,trios[1,]] - comp[,trios[2,]]
	e2 <- comp[,trios[2,]] - comp[,trios[3,]]
	e3 <- comp[,trios[3,]] - comp[,trios[1,]]
	tris <- which((apply(abs(e1),2,sum) == 2) & (apply(abs(e2),2,sum) == 2) & (apply(abs(e3),2,sum) == 2))
	return(list(composition=comp,tris=trios[,tris])) # n^2個の三角形ができる
}

my.divide.n(2)

my.integer.division <- function(n,m){
	p <- min(n,m)
	q <- max(n,m)
	tmp1 <- q %% p
	tmp2 <- q %/% p
	ret <- rep(tmp2,p)
	tmp3 <- floor(((1:p) * tmp1)/p)
	tmp4 <- c(0,diff(tmp3))
	tmp5 <- which(tmp4==1)
	ret[tmp5] <- ret[tmp5] + 1
	return (ret)
}

my.integer.division(3,22)


my.divide.tri <- function(tri,d.n=NULL){
	if(sum(abs(tri[[1]]-c(1,1,1)))==0){
		ret <- list()
		#ret[[1]] <- list()
		#ret[[1]][[1]] <- tri[[1]]
		#ret[[1]][[2]] <- tri[[2]]
		ret[[1]] <- list(tri[[1]],tri[[2]])
		return(ret)
	}
	if(min(tri[[1]])==1){
		snd <- sort(tri[[1]])[2]
		if(snd == 1){
			ret <- my.divide.tri.11(tri)
			return(ret)
		}else{
			ret <- my.divide.tri.1(tri)
			return(ret)
		}
	}
	if(is.null(d.n)){
		
		d.n <- my.divide.n(min(tri[[1]]))
	}
	edge.div <- list()
	edge.frac <- list()
	for(i in 1:3){
		edge.div[[i]] <- my.integer.division(min(tri[[1]]),tri[[1]][i])
		edge.frac[[i]] <- cumsum(c(0,edge.div[[i]]))
		edge.frac[[i]] <- edge.frac[[i]]/sum(edge.div[[i]])
	}
	ret <- list()
	for(i in 1:length(d.n$tris[1,])){
		first <- d.n$composition[,d.n$tris[1,i]]+1
		second <- d.n$composition[,d.n$tris[2,i]]+1
		third <- d.n$composition[,d.n$tris[3,i]]+1


		frac1 <- c(edge.frac[[1]][first[1]],edge.frac[[2]][first[2]],edge.frac[[3]][first[3]])
		frac1 <- frac1/sum(frac1)
		frac2 <- c(edge.frac[[1]][second[1]],edge.frac[[2]][second[2]],edge.frac[[3]][second[3]])
		frac2 <- frac2/sum(frac2)
		frac3 <- c(edge.frac[[1]][third[1]],edge.frac[[2]][third[2]],edge.frac[[3]][third[3]])
		frac3 <- frac3/sum(frac3)
		
		x1 <- frac1[1] * tri[[2]][,1]+frac1[2] * tri[[2]][,2]+frac1[3] * tri[[2]][,3]
		x2 <- frac2[1] * tri[[2]][,1]+frac2[2] * tri[[2]][,2]+frac2[3] * tri[[2]][,3]
		x3 <- frac3[1] * tri[[2]][,1]+frac3[2] * tri[[2]][,2]+frac3[3] * tri[[2]][,3]
		#ret[[i]] <- list()
		
		e1 <- d.n$composition[,d.n$tris[,i][1]] - d.n$composition[,d.n$tris[,i][2]]
		e2 <- d.n$composition[,d.n$tris[,i][2]] - d.n$composition[,d.n$tris[,i][3]]
		e3 <- d.n$composition[,d.n$tris[,i][3]] - d.n$composition[,d.n$tris[,i][1]]
		dir1 <- which(e1==0)
		dir2 <- which(e2==0)
		dir3 <- which(e3==0)

		n1 <- edge.div[[dir1]][min(d.n$composition[dir1,d.n$tris[,i][1]],d.n$composition[dir1,d.n$tris[,i][2]])+1]
		n2 <- edge.div[[dir2]][min(d.n$composition[dir2,d.n$tris[,i][2]],d.n$composition[dir2,d.n$tris[,i][3]])+1]
		n3 <- edge.div[[dir3]][min(d.n$composition[dir3,d.n$tris[,i][3]],d.n$composition[dir3,d.n$tris[,i][1]])+1]
		ord <- order(c(dir1,dir2,dir3))

		#ret[[i]][[1]] <- c(n1,n2,n3)[ord]
		#ret[[i]][[2]] <- cbind(x1,x2,x3)
		ret[[i]] <- list(c(n1,n2,n3)[ord],cbind(x1,x2,x3))
	}
	return(ret)
}
my.divide.tri.11 <- function(tri){ # tri[[1]] = {1,1,n!=1}
	ord <- order(tri)
	snd <- tri[ord[2]]
	mx <- tri[ord[3]]
	
	ret <- list()
	
	n.triangle <- mx
	for(i in 1:n.triangle){
		
		ret[[i]] <- list(c(1,1,1),cbind(x1,x2,x3))
	}
	return(ret)
}

my.divide.tri.1 <- function(tri){ # tri[[1]] = {1,n!=1,m!=1}
	ord <- order(tri)
	snd <- tri[ord[2]]
	mx <- tri[ord[3]]
	
	n.d <- my.integer.division(snd,mx)
	
	
	
}

my.plot.divtri <- function(q,face=TRUE,edge=TRUE){
	mat <- matrix(0,0,3)
	for(i in 1:length(q)){
		mat <- rbind(mat,t(q[[i]][[2]]))
	}
	plot3d(mat)
	if(face){
		col <- matrix(0,0,3)
		for(i in 1:length(q)){
			col <- rbind(col,q[[i]][[1]])
		}
		col <- col
		col <- t(t(col)/apply(col,2,max))
		col <- (1-col)

		for(i in 1:length(q)){
			triangles3d(t(q[[i]][[2]]),col=rgb(col[i,1],col[i,2],col[i,3]))
		}
	}

	if(edge){
		ed <- matrix(0,0,3)
		for(i in 1:length(q)){
			#segments3d(rbind(q[[i]][[2]][,1],q[[i]][[2]][,2]))
			#segments3d(rbind(q[[i]][[2]][,2],q[[i]][[2]][,3]))
			#segments3d(rbind(q[[i]][[2]][,3],q[[i]][[2]][,1]))
			ed <- rbind(ed,rbind(q[[i]][[2]][,1],q[[i]][[2]][,2]))
			ed <- rbind(ed,rbind(q[[i]][[2]][,2],q[[i]][[2]][,3]))
			ed <- rbind(ed,rbind(q[[i]][[2]][,3],q[[i]][[2]][,1]))
		}
		segments3d(ed)
	}
}

my.divide.tri.recursive <- function(tri){
	loop <- TRUE
	cnt <- 1
	
	while(loop){
		loop <- FALSE
		ret <- list()
		print("cnt")
		print(cnt)
		cnt <- cnt+1
		for(i in 1:length(tri)){
			if(sum(abs(tri[[i]][[1]]-c(1,1,1)))==0){
				cnt2 <- length(ret)+1
				#ret[[cnt2]] <- list()
				#ret[[cnt2]][[1]] <- tri[[i]][[1]]
				#ret[[cnt2]][[2]] <- tri[[i]][[2]]
				#ret[[cnt2]] <- list(tri[[i]][[1]],tri[[i]][[2]])
				ret[[cnt2]] <- tri[[i]]
				print("in")
				print(abs(tri[[i]][[1]]-c(1,1,1)))

			}else{
				print(tri[[i]])
				out <- my.divide.tri(tri[[i]])
				print(out)
				for(j in 1:length(out)){
					cnt2 <- length(ret)+1
					#ret[[cnt2]] <- list()
					#ret[[cnt2]][[1]] <- out[[j]][[1]]
					#ret[[cnt2]][[2]] <- out[[j]][[2]]
					#ret[[cnt2]] <- list(out[[j]][[1]],out[[j]][[2]])
					ret[[cnt2]] <- out[[j]]
				}
				#ret <- my.divide.tri.recursive(tri=out,ret=ret)
				#print(j)
				print(length(ret))
				print("----")
				print(ret)
				print("===")
				#loop <- TRUE
			}
		}
		#print(tri)
		tri <- ret
		#print(tri)
	}

	return(tri)
}

tri.org <- list(le = c(2,2,3),x = diag(rep(3,3)))
tri <- list()
tri[[1]] <- list(tri.org[[1]],tri.org[[2]])

#tri[[1]] <- list()
#tri[[1]][[1]] <- tri.org[[1]]
#tri[[1]][[2]] <- tri.org[[2]]


tmpout <- my.divide.tri(tri[[1]])

my.plot.divtri(tmpout)
tmpout2 <- my.divide.tri.recursive(tri)

my.plot.divtri(tmpout2)


