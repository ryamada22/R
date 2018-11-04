library(onion)

q <- rquat(1)
p <- rquat(1)
w <- rquat(1)

Re(p) <- 0
Re(w) <- 0

q <- q/Mod(q)
Mod(q)
sum(as.vector(q * p * Conj(q)) * as.vector(w))

sum(as.vector(q * p) * as.vector(w * q))
q.v <- as.vector(q)

matrix(q.v,nrow=1) %*% (my.mat1(p) %*% my.mat2(w)) %*% matrix(q.v,ncol=1)
matrix(q.v,nrow=1) %*% ((my.mat1(p) %*% my.mat2(w))) %*% matrix(q.v,ncol=1)

theta <- runif(1)
psi < -runif(1)

p <- Hi * sin(theta) * cos(psi) + Hj * sin(theta) * sin(psi) + Hk * cos(theta)
w <- p
sum(as.vector(q * p * Conj(q)) * as.vector(w))

sum(as.vector(q * p) * as.vector(w * q))
q.v <- as.vector(q)

matrix(q.v,nrow=1) %*% (my.mat1(p) %*% my.mat2(w)) %*% matrix(q.v,ncol=1)
(my.mat1(p) %*% my.mat2(w))

u1 <- matrix(q.v,nrow=1) %*% my.mat1(p)

u2 <- my.mat2(w) %*% matrix(q.v,ncol=1)

u1 %*% u2


my.mat1 <- function(p){
	ret <- matrix(Re(p),4,4)
	ret[1,2] <- (-1) * i(p)
	ret[1,3] <- (-1) * j(p)
	ret[1,4] <- (-1) * k(p)
	ret[2,1] <- i(p)
	ret[2,3] <- k(p)
	ret[2,4] <- (-1) * j(p)
	ret[3,1] <- j(p)
	ret[3,2] <- (-1) * k(p)
	ret[3,4] <- i(p)
	ret[4,1] <- k(p)
	ret[4,2] <- j(p)
	ret[4,3] <- (-1) * i(p)
	return(ret)
}
my.mat2 <- function(p){
	ret <- matrix(Re(p),4,4)
	ret[1,2] <- (-1) * i(p)
	ret[1,3] <- (-1) * j(p)
	ret[1,4] <- (-1) * k(p)
	ret[2,1] <- i(p)
	ret[2,3] <- (-1) * k(p)
	ret[2,4] <-  j(p)
	ret[3,1] <- j(p)
	ret[3,2] <-  k(p)
	ret[3,4] <- (-1) *i(p)
	ret[4,1] <- k(p)
	ret[4,2] <- (-1) *j(p)
	ret[4,3] <-  i(p)
	return(ret)
}

my.mat12 <- function(theta,psi){
	ret <- matrix(0,4,4)
	ret[1,1] <- -1
	ret[2,2] <- cos(theta)^2-sin(theta)^2*cos(2*psi)
	ret[3,3] <- sin(theta)^2*cos(2*psi) + cos(theta)^2
	ret[4,4] <- -cos(2*theta)
	ret[2,3] <- ret[3,2] <- -sin(theta)^2*sin(2*psi)
	ret[2,4] <- ret[4,2] <- -cos(psi) * sin(2*theta)
	ret[3,4] <- ret[4,3] <- -sin(psi) * sin(2*theta)
	return(ret)
}

theta <- runif(1)
psi <- runif(1)

p <- Hi * sin(theta) * cos(psi) + Hj * sin(theta) * sin(psi) + Hk * cos(theta)

my.mat1(p) %*% my.mat2(p)
my.mat12(theta,psi)
round(my.mat12(theta,psi) - my.mat1(p) %*% my.mat2(p),5)
eigen(my.mat12(theta,psi))

my.mat1(p)
p
my.mat2(p)
