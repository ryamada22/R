# 三角形を(x,y,z)分割する
library(partitions)
library(combinat)
library(rgl)
# 同一平面上にある４点が作る直線の交点座標を出す
my.2d.intersect <- function(v1,v2,v3,v4){
	M <- matrix(c(v2[1]-v1[1],v3[1]-v4[1],v2[2]-v1[2],v3[2]-v4[2]),byrow=TRUE,2,2)
	b <- matrix(c(-v1[1]+v3[1],-v1[2]+v3[2]),ncol=1)
	a <- solve(M,b)
	return(a)
}

my.3d.intersect.tri <- function(v1,v2,v3,v4){
	a <- my.2d.intersect(v1[1:2],v2[1:2],v3[1:2],v4[1:2])
	t <- a[1]
	v1 + (v2-v1) * t
}

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

my.div.series <- function(d){
	ret <- list()
	ret[[1]] <- d
	for(i in 2:length(d)){
		n <- length(ret[[i-1]])
		tmp1 <- ret[[i-1]][1:(n-1)]
		tmp2 <- ret[[i-1]][2:n]
		ret[[i]] <- apply(cbind(tmp1,tmp2),1,min)
	}
	return(ret)
}
id <- my.integer.division(3,22)
id
my.div.series(id)

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
	
	# 等分割の場合の、三角形のリストと、そのBarycentric coordinates相当整数値
	
	if(is.null(d.n)){
		
		d.n <- my.divide.n(min(tri[[1]]))
	}
	
	# 整数分割
	edge.div <- list()
	edge.frac <- list()
	div.series <- list()
	for(i in 1:3){
		edge.div[[i]] <- my.integer.division(min(tri[[1]]),tri[[1]][i])
		div.series[[i]] <- my.div.series(edge.div[[i]])
		edge.frac[[i]] <- cumsum(c(0,edge.div[[i]]))
		edge.frac[[i]] <- edge.frac[[i]]/sum(edge.div[[i]])
	}
	

	
	# edge上の点座標の計算
	nv <- max(d.n$tris)
	divnum <- max(d.n$composition)
	v <- c(1:3,1:3,1)
	edge.x <- list()
	for(i in 1:3){
		edge.x[[i]] <- matrix(0,3,divnum+1)
		for(j in 1:(divnum+1)){
			edge.x[[i]][,j] <- (1-edge.frac[[i]][j]) * tri[[2]][,v[i+1]] + edge.frac[[i]][j] * tri[[2]][,v[i+2]]
			print(edge.x[[i]][,j])
		}
	}
	# 分割後の頂点座標の計算

	x <- matrix(0,3,nv)
	v.x <- list()
	
	vv <- matrix(0,0,3)
	for(i in 1:nv){
		tmp <- d.n$composition[,i]
		num.zero <- length(which(tmp==0))
		if(num.zero==2){
			nonzeroid <- which(d.n$composition[,i]!=0)
			x[,i] <- tri[[2]][,nonzeroid]
						print("$$$")
			print(x[,i])
		}else if(num.zero==1){
			zeroid <- which(d.n$composition[,i]==0)
			val <- tmp[v[zeroid+2]]
			x[,i] <- edge.x[[zeroid]][,val+1]
			print("&&&")
			print(x[,i])
		}else{
			v.x[[i]] <- list()
			three.lines <- list()
			for(j in 1:3){
				one <- edge.x[[v[j+1]]][,tmp[j]+1]
				another <- edge.x[[v[j+2]]][,divnum-tmp[j]+1]
				#print(one)
				#print(another)
				three.lines[[j]] <- cbind(one,another)
			}
			tmp <- cbind(three.lines[[1]],three.lines[[2]],three.lines[[3]])
			cross1 <- my.3d.intersect.tri(three.lines[[1]][,1],three.lines[[1]][,2],three.lines[[2]][,1],three.lines[[2]][,2])
			cross2 <- my.3d.intersect.tri(three.lines[[2]][,1],three.lines[[2]][,2],three.lines[[3]][,1],three.lines[[3]][,2])
			cross3 <- my.3d.intersect.tri(three.lines[[3]][,1],three.lines[[3]][,2],three.lines[[1]][,1],three.lines[[1]][,2])

			x[,i] <- apply(cbind(cross1,cross2,cross3),1,mean)

		}


		
		
	}
	x
	plot3d(t(x))
	spheres3d(t(x),radius=0.1)

	ret <- list()
	for(i in 1:length(d.n$tris[1,])){
		# ３頂点のID
		vid <- d.n$tris[,i]
		# ３頂点の座標
		vx <- x[,vid]
		
		# ３辺の長さ
		# 第１辺(第１頂点の対辺)
		tmp <- d.n$composition[,vid]
		tmptmp1 <- tmp[,1]-tmp[,2]
		tmptmp2 <- tmp[,2]-tmp[,3]
		tmptmp3 <- tmp[,3]-tmp[,1]
		
		dir1 <- which(tmptmp1==0)
		level1 <- tmptmp1[dir1]
		mm1 <- tmptmp1[v[dir1+2]]
		n1 <- div.series[[dir1]][[level1+1]][mm1+1]
		dir2 <- which(tmptmp2==0)
		level2 <- tmptmp2[dir2]
		mm2 <- tmptmp2[v[dir2+2]]
		n2 <- div.series[[dir2]][[level2+1]][mm2+1]
		dir3 <- which(tmptmp3==0)
		level3 <- tmptmp3[dir3]
		mm3 <- tmptmp3[v[dir3+2]]
		n3 <- div.series[[dir3]][[level3+1]][mm3+1]

		ret[[i]] <- list(c(n1,n2,n3),vx)
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
			#triangles3d(t(q[[i]][[2]]),col=rgb(col[i,1],col[i,2],col[i,3]))
			triangles3d(t(q[[i]][[2]]),col=i)
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

tri.org <- list(le = c(5,7,8),x = diag(rep(3,3)))
tri <- list()
tri[[1]] <- list(tri.org[[1]],tri.org[[2]])

#tri[[1]] <- list()
#tri[[1]][[1]] <- tri.org[[1]]
#tri[[1]][[2]] <- tri.org[[2]]


tmpout <- my.divide.tri(tri[[1]])

my.plot.divtri(tmpout)
for(i in 1:length(tmpout)){
	ctr <- apply(tmpout[[i]][[2]],1,mean)
	txt = paste("",i)
	text3d(ctr+0.1,text=txt)
}



tmpout2 <- my.divide.tri.recursive(tri)

my.plot.divtri(tmpout2)

