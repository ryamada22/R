library(devtools)
library(rgl)
library(knitr)
library(e1071)
library(ape)
library(igraph)
library(gtools)
install_github("ryamada22/Ronlyryamada")
library(Ronlyryamada)
library(RFOC)
library(onion)


my.partition.tris <- function(tri.x){
	tris <- tri.x$tris
	X <- tri.x$X
	col <- tri.x$col
	tris <- t(apply(tris,1,sort))
	pairs <- rbind(tris[,1:2],tris[,2:3],tris[,c(1,3)])
	pairs2 <- t(apply(pairs,1,sort))
	pairs2 <- unique(pairs)
	new.v <- (1:length(pairs2[,1])) + max(pairs2)
	new.X <- (X[pairs2[,1],] + X[pairs2[,2],] ) /2
	ret.X <- rbind(X,new.X)
	new.tris <- matrix(0,0,3)
	new.col <- c()
	for(i in 1:length(tris[,1])){
		tri.now <- tris[i,]
		tri.add1 <- which(pairs2[,1]==tri.now[1] & pairs2[,2]==tri.now[2])+ max(pairs2)
		tri.add2 <- which(pairs2[,1]==tri.now[2] & pairs2[,2]==tri.now[3])+ max(pairs2)
		tri.add3 <- which(pairs2[,1]==tri.now[1] & pairs2[,2]==tri.now[3])+ max(pairs2)
		
		tmp <- rbind(c(tri.now[1],tri.add1,tri.add3),c(tri.now[2],tri.add1,tri.add2),c(tri.now[3],tri.add2,tri.add3),c(tri.add1,tri.add2,tri.add3))
		new.tris <- rbind(new.tris,tmp)
		if(col[i]==1){
			tmp.col <- c(2,2,2,1)
		}else{
			tmp.col <- c(1,1,1,2)
		}
		new.col <- c(new.col,tmp.col)
	}
	new.tris <- t(apply(new.tris,1,sort))
	
	ret.X <- ret.X/sqrt(apply(ret.X^2,1,sum))
	return(list(tris=new.tris,X=ret.X,col=new.col))
}

my.plot.eultri <- function(tri,alpha=0.5){
	plot3d(tri$X,cex=0.01)
	n.Face <- length(tri$tris[,1])
	col <- rep(2:3,n.Face/2)
	col <- c(rbind(tri$col,tri$col,tri$col))
	tmp.ord <- c(t(tri$tris))
	triangles3d(tri$X[tmp.ord,],col=col+1,alpha=alpha)
}


tris8=matrix(c(1,2,3,1,3,4,1,4,5,1,5,2,2,3,6,3,4,6,4,5,6,5,2,6),byrow=TRUE,ncol=3)
X8 = matrix(c(0,0,1,1,0,0,0,1,0,-1,0,0,0,-1,0,0,0,-1),byrow=TRUE,ncol=3)
col8 = c(1,2,1,2,2,1,2,1)
init.tri.x <- list(tris=tris8,X=X8,col=col8)

n.step <- 4
tri.x.list <- list()
tri.x.list[[1]] <- init.tri.x
for(i in 2:n.step){
	tri.x.list[[i]] <- my.partition.tris(tri.x.list[[i-1]])
}


my.plot.eultri(tri.x.list[[n.step]])


my.Tri2Ms <- function(tri){
  xyz <- tri$X
  edge <- black.tri.id <- which(tri$col==1)
  el <- rbind(tri$tris[black.tri.id,1:2],tri$tris[black.tri.id,2:3],tri$tris[black.tri.id,c(3,1)])
  face <- tri$tris
  return(list(xyz=xyz,edge=el,face=face))
}
# The faceid-th face triangle is the outer most and the planar graph is drawn on the plane. Ms is trianglular mesh on the sphere. faceid is an integer. k is magnification factor of the outer most triangle.
my.planegraph <- function(Ms,faceid,k=1.1){
  face.v <- Ms$face[faceid,]
  
  face.xyz <- Ms$xyz[face.v,]
  direc <- apply(face.xyz,2,mean)
  direc <- direc/sqrt(sum(direc^2))
  xyz <- my.3drotation(Ms$xyz,direc,c(0,0,1))
  rang.complex <- my.Sp2Plane(xyz)
  max.R <- max(rang.complex$r)
  inside.max.R <- max(rang.complex$r[-face.v])
  
  face.v.ang <- rang.complex$ang[face.v]
  tmp1 <- cos((face.v.ang[1]-face.v.ang[2])/2)
  tmp2 <- cos((face.v.ang[1]-face.v.ang[3])/2)
  tmp3 <- cos((face.v.ang[2]-face.v.ang[3])/2)
  k2 <- min(abs(c(tmp1,tmp2,tmp3)))
  rang.complex$r[Ms$face[faceid,]] <- inside.max.R/k2 * k
  tmp <- rang.complex$r * exp(1i * rang.complex$ang)
  return(cbind(Re(tmp),Im(tmp)))
}

# Plot the planar graph
# When edge=NULL, all edges are drawn
# When edge is given as 2-column matrix, the given edges only are drawn (overlap=FALSE) or are highlighted (overlap=TRUE).
my.plot.plane.graph <- function(Mr,faceid,edge=NULL,overlap=TRUE){
  vs <- sort(Mr$face[faceid,])
  if(is.null(edge)){
    edge <- Mr$edge
    plot(Mr$xyz,pch=20,cex=0.1)
    edge1 <- edge
    segments(Mr$xyz[edge1[,1],1],Mr$xyz[edge1[,1],2],Mr$xyz[edge1[,2],1],Mr$xyz[edge1[,2],2])
  }else{
    if(overlap){
      edge.ori <- Mr$edge
      plot(Mr$xyz,pch=20,cex=0.1)
      edge1 <- edge.ori
      segments(Mr$xyz[edge1[,1],1],Mr$xyz[edge1[,1],2],Mr$xyz[edge1[,2],1],Mr$xyz[edge1[,2],2],col="gray")
      edge1 <- edge
      segments(Mr$xyz[edge1[,1],1],Mr$xyz[edge1[,1],2],Mr$xyz[edge1[,2],1],Mr$xyz[edge1[,2],2],col="red")
    }else{
      plot(Mr$xyz,pch=20,cex=0.1)
      edge1 <- edge
      segments(Mr$xyz[edge1[,1],1],Mr$xyz[edge1[,1],2],Mr$xyz[edge1[,2],1],Mr$xyz[edge1[,2],2])
    }
  }
}


# v1,v2,v3 are points on a unit sphere
# return 0 means counter-clockwise when looking the triangle from outside of sphere
my.tri.direction <- function(v1,v2,v3){
	d12 <- v2-v1
	d13 <- v3-v1
	axis <- c(d12[2]*d13[3]-d12[3]*d13[2],d12[3]*d13[1]-d12[1]*d13[3],d12[1]*d13[2]-d12[2]*d13[1])
	
	#ctr.direction <- apply(cbind(v1,v2,v3),1,mean)
	tmp <- sum(v1 * axis)
	if(tmp>=0){
		ret <- 1
	}else{
		ret <- 0
	}
	return(ret)
}

my.direction.eultri.xxx <- function(tri,black=1){
  n.tri <- length(tri$tris[,1])
  n.edge <- n.tri*3/2
  black.tri.id <- which(tri$col==black)
  el <- rbind(tri$tris[black.tri.id,1:2],tri$tris[black.tri.id,2:3],tri$tris[black.tri.id,c(3,1)])
  
  el.dir <- rep(0,length(el[,1]))
  
  first.tri.id <- which(tri$col==black)[1]
  first.tri <- tri$tris[first.tri.id,]
  cnt <- 1
  el[cnt,] <- c(fisrt.tri[1],first.tri[2])
  el[cnt+1,] <- c(first.tri[2],first.tri[3])
  el[cnt+2,] <- c(first.tri[3],first.tri[1])
  cnt <- cnt+3
  
}
my.direction.eultri <- function(tri,black=1){
  n.tri <- length(tri$tris[,1])
  n.edge <- n.tri*3/2
  el <- matrix(0,n.edge,2)
  black.tris <- tri$tris[which(tri$col==black),]
  cnt <- 1
  for(i in 1:length(black.tris[,1])){
    v <- black.tris[i,]
    print("v")
    print(v)
    M <- tri$X[v,]
    d2 <- det(M)
    d <- my.tri.direction(M[1,],M[2,],M[3,])
    print("d and d2")
    print(c(d,d2))
    print("d")
    print(d)
    if(d == 1){
    #if(d2 >0){
			print("èá")
      el[cnt,] <- c(v[1],v[2])
      el[cnt+1,] <- c(v[2],v[3])
      el[cnt+2,] <- c(v[3],v[1])
      print(el[cnt,])
      print(el[cnt+1,])
      print(el[cnt+2,])
    }else{
			print("ãt")
      el[cnt,] <- c(v[2],v[1])
      el[cnt+1,] <- c(v[3],v[2])
      el[cnt+2,] <- c(v[1],v[3])
      print(el[cnt,])
      print(el[cnt+1,])
      print(el[cnt+2,])
    }
    cnt <- cnt+3

  }
  return(el)
}


tri <- tri.x.list[[n.step]]
#tri <- tri.x.list[[1]]


el <- my.direction.eultri(tri,black=1)


g <- graph.edgelist(el)
plot(g)


rho <- sample(1:length(tri$X[,1]),1)
Z <- shortest.paths(g,rho,mode="out")


my.plot.eultri(tri)
spheres3d(tri$X,radius=Z*0.01)
spheres3d(tri$X[rho,],radius=0.1,color="red")


my.eultree <- function(tri,rho=1,black=1){
  el <- my.direction.eultri(tri,black)
  g <- graph.edgelist(el)
  Z <- shortest.paths(g,rho,mode="out")
  n.v <- length(tri$X[,1])
  n.f <- length(tri$tris[,1])
  col.lab <- unique(tri$col)
  white <- col.lab[which(col.lab!=black)]
  
  tmpX <- matrix(tri$X[tri$tris,1],ncol=3)
  tmpY <- matrix(tri$X[tri$tris,2],ncol=3)
  tmpZ <- matrix(tri$X[tri$tris,3],ncol=3)
  new.X <- cbind(apply(tmpX,1,mean),apply(tmpY,1,mean),apply(tmpZ,1,mean))
  white.tri <- which(tri$col==white)
  black.tri <- which(tri$col==black)
  v.id.black.tri <- tri$tris[black.tri,]
  
  unfl.v <- which(tri$col==white)
  fl.v <- which(tri$col==black)
  #unfl.el <- matrix(0,length(unfl.v)*2,2)
  #fl.el <- matrix(0,length(fl.v),2)
  unfl.el <- matrix(0,0,2)
  fl.el <- matrix(0,0,2)
  unfl.cnt <- 1
  fl.cnt <- 1
  #flgs <- matrix(0,length(fl.v),2)
  flgs <- matrix(0,0,2)
  for(i in 1:length(white.tri)){
    v <- tri$tris[white.tri[i],]
    v. <- c(v,v)
    for(j in 1:3){
      this.el <- v.[c(j,j+1)]
      # direction check of this.el
      tmp1 <- el[,1] - this.el[1]
      tmp2 <- el[,2] - this.el[2]
      tmp3 <- abs(tmp1) + abs(tmp2)
      tmp4 <- which(tmp3==0)
      if(length(tmp4)>0){
				
			}else{
				this.el <- this.el[c(2,1)]
			}
      this.Z <- Z[this.el]
      m <- this.Z[1]
      n <- this.Z[2]
      if(n==(m+1)){
        # v.id of ctr of white tri is i + n.v 
        #unfl.el[unfl.cnt,] <- c(i+n.v,this.el[2])
        unfl.el <- rbind(unfl.el,c(i+n.v,this.el[2]))
        unfl.cnt <- unfl.cnt + 1
      }else if(n<=m){
        # v.id of ctr of black tri is i + 2 * n.vou
        #fl.el[fl.cnt,] <- c(i+n.v,i+2*n.v)
        tmp1 <- v.id.black.tri - this.el[1]
        tmp2 <- v.id.black.tri - this.el[2]
        tmp11 <- apply(tmp1,1,prod)
        tmp22 <- apply(tmp2,1,prod)
        black.face <- which(tmp11==0 & tmp22==0)
        #print(v)
        #print(c(i,this.el))
        #print(black.face)
        fl.el <- rbind(fl.el,c(i+n.v,black.face+n.v+n.f/2))
        #flgs[fl.cnt,] <- this.el
        flgs <- rbind(flgs,this.el)
        fl.cnt <- fl.cnt + 1
      }
    }
  }
  ret.X <- rbind(tri$X,new.X[white.tri,],new.X[black.tri,])
  ret.el <- rbind(el,unfl.el,fl.el)
  ret.el.2 <- rbind(unfl.el,fl.el)
  ret.g <- graph.edgelist(ret.el)
  # flgs
  return(list(g=ret.g,el=ret.el,el.2 = ret.el.2,X=ret.X,flgs=flgs,flg.v=(2*n.v+1):length(ret.X[,1]),flg.e = (length(ret.el[,1])-length(fl.v)+1):length(ret.el[,1])))

}


out <- my.eultree(tri,rho=rho,black=1)


my.plot.eultri(tri,alpha=0.03)

spheres3d(out$X,radius=0.005,color="blue")
segments3d(out$X[t(out$el.2),],lwd=1)
spheres3d(tri$X,radius=0.05)


