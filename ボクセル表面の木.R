# ボクセルはその座標最小頂点(x,y,z)座標で登録する

# 各ボクセルは６面を持つ
# 面(正方形)はその座標最小頂点(x,y,z)と、面の方向で区別する
# xy平面、yz平面、zx平面に広がるものを、1,2,3と呼び分ける

# 各頂点は４頂点を持つ

Vox.list <- as.matrix(expand.grid(0:2,0:2,0:2))

my.vox.faces <- function(xyz){
	rbind(c(xyz,1),c(xyz,2),c(xyz,3),c(xyz+c(1,0,0),2),c(xyz+c(0,1,0),3),c(xyz+c(0,0,1),1))
}

my.face.nodes <- function(xyzw){
	xyz <- xyzw[1:3]
	w <- xyzw[4]
	if(w==1){
		ret <- rbind(xyz,xyz+c(1,0,0),xyz+c(0,1,0),xyz+c(1,1,0))
	}else if(w==2){
		ret <- rbind(xyz,xyz+c(0,1,0),xyz+c(0,0,1),xyz+c(0,1,1))
	}else{
		ret <- rbind(xyz,xyz+c(0,0,1),xyz+c(1,0,0),xyz+c(1,0,1))
	}
	return(ret)
}

my.onlyone <- function(x){
	n <- length(x[,1])
	a <- !duplicated(x)
	b <- !duplicated(x[n:1,])[n:1]
	
	return(list(a & b,x[a & b,]))
}

my.vox.surface <- function(Vlist){
	faces <- matrix(0,0,4)
	
	for(i in 1:length(Vlist[,1])){
		faces <- rbind(faces,my.vox.faces(Vlist[i,]))
		
	}
	omoteura <- rep(c(0,0,0,1,1,1),length(faces[,1]))
	onlyones <- my.onlyone(faces)
	#print(onlyones[[1]])
	return(list(onlyones[[2]],omoteura[which(onlyones[[1]])]))
}

facesout <- my.vox.surface(Vox.list)
faces <- facesout[[1]]
faces.omoteura <- facesout[[2]]

my.surface.nodes <- function(faces){
	vs <- matrix(0,0,3)
	for(i in 1:length(faces[,1])){
		vs <- rbind(vs,my.face.nodes(faces[i,]))
	}
	return(unique(vs))
}

nodes <- my.surface.nodes(faces)

library(rgl)
plot3d(nodes)

library(igraph)
my.surface.graph <- function(nodes){
	d <- as.matrix(dist(nodes))
	d <- matrix(as.numeric(d==1),ncol=length(nodes[,1]))
	d
}

g <- graph.adjacency(my.surface.graph(nodes),mode="undirected")
plot(g)

edges <- t(apply(get.edgelist(g),1,sort))

my.face.nodesid <- function(faces,faces.omoteura,nodes){
	ret <- matrix(0,length(faces[,1]),4)
	for(i in 1:length(faces[,1])){
		nodes.xyz <- my.face.nodes(faces[i,])
		nodes.xyz <- nodes.xyz[c(1,2,4,3),] # 反時計回り
		if(faces[i,4]==1){
			if(faces.omoteura[i]==0){
				#nodes.xyz <- nodes.xyz[c(1,4,3,2),]
			}else{
				nodes.xyz <- nodes.xyz[c(1,4,3,2),]
			}
		}else if(faces[i,4]==2){
			if(faces.omoteura[i]==0){
				#nodes.xyz <- nodes.xyz[c(1,4,3,2),]
			}else{
				nodes.xyz <- nodes.xyz[c(1,4,3,2),]
			}
		}else{
			if(faces.omoteura[i]==0){
				#nodes.xyz <- nodes.xyz[c(1,4,3,2),]
			}else{
				nodes.xyz <- nodes.xyz[c(1,4,3,2),]
			}
		}
		nodes.id <- rep(0,4)
		for(j in 1:4){
			tmp <- t(nodes) - nodes.xyz[j,]
			tmp2 <- apply(tmp^2,2,sum)
			nodes.id[j] <- which(tmp2==0)
			ret[i,] <- nodes.id
		}
	}
	return(ret)
}

face.nodes <- my.face.nodesid(faces,faces.omoteura,nodes)

rootid <- 10

shdist <- distances(g)

rootdist <- shdist[rootid,]

my.tree.edge <- function(vals){ # 反時計回り
	tmp <- which(vals == max(vals))
	if(length(tmp)==2){
		ret <- tmp
		print("diagonal")
	}else{
		ret <- c(tmp,tmp+1)
		if(ret[2] > 4){
			ret[2] <- 1
		}
	}
	return(ret)
}

my.tree.edge(c(1,2,3,2))

my.select.face.edge <- function(nodeid,rootdist){
	vals <- rootdist[nodeid]
	tmp <- my.tree.edge(vals)
	return(nodeid[tmp])
}

my.tree.edges <- function(face.nodes,rootdist){
	tmp <- apply(face.nodes,1,my.select.face.edge,rootdist)
	ret <- matrix(tmp,nrow=2)
	return(t(ret))
}

tree.edges <- my.tree.edges(face.nodes,rootdist)

my.draw.surface <- function(nodes,root,edges,tree.edges=Null,axes = FALSE,radius.root = 0.1){
	plot3d(nodes,axes=axes,xlab="x",ylab="y",zlab="z")
	tmp <- c(t(edges))
	#segments3d(nodes[tmp,])
	if(!is.null(tree.edges)){
		tmp <- c(t(tree.edges))
		segments3d(nodes[tmp,],col=2)
	}
	spheres3d(nodes[root,],col=3,radius=radius.root)
	
}

my.draw.surface(nodes,rootid,edges,tree.edges=tree.edges)

one.nodes <- which(rootdist==1)
for(i in 1:length(one.nodes)){
	segments3d(rbind(nodes[rootid,],nodes[one.nodes[i],]),col="purple")
}


Vox.list1 <- as.matrix(expand.grid(0:4,0:4,0:4))
Vox.list2 <- as.matrix(expand.grid(2:5,2:5,2:5))

Vox.list3 <- as.matrix(expand.grid(1:10,3:4,4:5))
Vox.list4 <- as.matrix(expand.grid(2:5,3:20,4:5))
Vox.list <- unique(rbind(Vox.list1,Vox.list2,Vox.list3,Vox.list4))

n.step <- 1000
xxx <- matrix(rep(0,3),ncol=3)

for(i in 1:n.step){
	r <- sample(1:length(xxx[,1]),1)
	p <- sample(1:3,1)
	tmp <- xxx[r,]
	tmp[p] <- tmp[p] + 1
	xxx <- rbind(xxx,tmp)
	xxx <- unique(xxx)
}
Vox.list <- unique(xxx)
