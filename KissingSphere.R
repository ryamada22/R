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
	return(list(x=x,eigen.out = out))
}

library(RFOC)
library(igraph)
library(knitr)
library(tagcloud)
library(e1071)


my.rsphere <- function(n,d=3){
	X <- matrix(rnorm(n*d),ncol=d)
	X <- X/sqrt(apply(X^2,1,sum))
	X
}
n1 <- 300
n2 <- 300
r1 <- 1
r2 <- 0.7
S1 <- my.rsphere(n1) 
S2 <- my.rsphere(n2) 

SS1 <- S1 %*% t(S1)
SS2 <- S2 %*% t(S2)

SS1[which(SS1>1)] <- 1
SS1[which(SS1< -1)] <- -1
SS2[which(SS2>1)] <- 1
SS2[which(SS2 < -1)] <- -1

D1 <- acos(SS1)* r1
D2 <- acos(SS2)* r2

D <- matrix(0,n1+n2,n1+n2)
D[1:n1,1:n1] <- D1
D[(n1+1):(n1+n2),(n1+1):(n1+n2)] <- D2

S <- matrix(0,n1+n2,n1+n2)
S[1:n1,1:n1] <- SS1*r1
S[(n1+1):(n1+n2),(n1+1):(n1+n2)] <- SS2*r2


A <- sample(1:n1,1)
B <- sample(1:n2,1)

for(i in 1:n1){
	for(j in 1:n2){
		d1 <- D1[i,A]
		d2 <- D2[j,B]
		D[i,n1+j] <- D[n1+j,i] <- d1 + d2
		s1 <- SS1[i,A]*r1
		s2 <- SS2[j,B]*r2
		S[i,n1+j] <- S[n1+j,i] <- s1+s2
	}
}

image(D)

eoutD <- eigen(D)
eoutD1 <- eigen(D1)
eoutS <- eigen(S)
eoutS1 <- eigen(SS1)

plot3d(eoutS[[2]][,c(1,3,4)])
open3d()
plot3d(eoutS[[2]][,c(2,5,6)])


col1 <- eoutS[[2]][,599]
col2 <- eoutS[[2]][,600]

plot3d(eoutS[[2]][,c(1,3,4)])
spheres3d(eoutS[[2]][,c(1,3,4)],color=round((col2-min(col2))/(max(col2)-min(col2)) * 5)+1,radius=0.05)

spheres3d(eoutS[[2]][A,c(1,3,4)],color=2,radius=0.07)


plot3d(eoutS[[2]][,c(2,5,6)])
spheres3d(eoutS[[2]][,c(2,5,6)],color=round((col1-min(col1))/(max(col1)-min(col1)) * 5)+1,radius=0.05)
spheres3d(eoutS[[2]][B+n1,c(2,5,6)],color=2,radius=0.07)


# ‚¾‚é‚Ü

n1 <- 300
n2 <- 300
r1 <- 1
r2 <- 1
S1 <- my.rsphere(n1) 
S2 <- my.rsphere(n2) 


# x = p‚ÅŒğ‚í‚é
p <- 0.7
S1. <- S1[which(S1[,1] < p),]
S2. <- S2[which(S2[,1] > -p),]

# x = p‚É‰~‚ğ”¼Œasqrt(1-p^2)‚Ì‰~‚ğì‚é
n3 <- 100
t <- seq(from=0,to=n3-1,length=n3)/n3 * 2*pi
R <- sqrt(1-p^2)
L <- cbind(rep(p,n3),R*cos(t),R*sin(t))
S2.. <- S2.
S2..[,1] <- S2..[,1] + 2*p
plot3d(rbind(S1.,S2..,L))

SS1. <- S1. %*% t(S1.)
SS2. <- S2. %*% t(S2.)

LS1 <- L %*% t(S1.)
LS2 <- L %*% t(S2.)

SS1.[which(SS1.>1)] <- 1
SS1.[which(SS1.< -1)] <- -1
SS2.[which(SS2.>1)] <- 1
SS2.[which(SS2. < -1)] <- -1

n1. <- length(SS1.[,1])
n2. <- length(SS2.[,1])

D1. <- acos(SS1.)* r1
D2. <- acos(SS2.)* r2

D. <- matrix(0,n1.+n2.,n1.+n2.)
D.[1:n1.,1:n1.] <- D1.
D.[(n1.+1):(n1.+n2.),(n1.+1):(n1.+n2.)] <- D2.

S. <- matrix(0,n1.+n2.,n1.+n2.)
S.[1:n1.,1:n1.] <- SS1.*r1
S.[(n1.+1):(n1.+n2.),(n1.+1):(n1.+n2.)] <- SS2.*r2


#A <- sample(1:n1,1)
#B <- sample(1:n2,1)

for(i in 1:n1.){
	for(j in 1:n2.){
		s1 <- LS1[,i]
		s2 <- LS2[,j]
		d1 <- acos(s1)
		d2 <- acos(s2)
		maxid <- which.max(d1+d2)
		s12 <- (s1+s2)[maxid]
		d12 <- (d1+d2)[maxid]
		D[i,n1.+j] <- D[n1.+j,i] <- d12
		
		S[i,n1.+j] <- S[n1.+j,i] <- s12
	}
}

image(D.)

eoutD. <- eigen(D.)
eoutD1. <- eigen(D1.)
eoutS. <- eigen(S.)
eoutS1. <- eigen(SS1.)

plot3d(eoutS.[[2]][,c(1,4,6)])
open3d()
plot3d(eoutS.[[2]][,c(2,3,5)])


col1. <- eoutS.[[2]][,n1.+n2.-1]
col2. <- eoutS.[[2]][,n1.+n2.]

plot3d(eoutS.[[2]][,c(1,4,6)])
spheres3d(eoutS.[[2]][,c(1,4,6)],color=round((col2.-min(col2.))/(max(col2.)-min(col2.)) * 5)+1,radius=0.01)



plot3d(eoutS.[[2]][,c(2,3,5)])
spheres3d(eoutS.[[2]][,c(2,3,5)],color=round((col1.-min(col1.))/(max(col1.)-min(col1.)) * 5)+1,radius=0.01)

