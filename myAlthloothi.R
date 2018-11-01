npt <- 20
d <- 3

X <- matrix(rnorm(npt*d),ncol=d)

library(GPArotation)

R <- Random.Start(3)

X. <- X %*% R

X. <- X. + rnorm(length(X.),0,0.0001)

M <- t(X.) %*% X

N <- matrix(0,4,4)

N[1,1] <- M[1,1] + M[2,2] + M[3,3]
N[1,2] <- M[2,3] - M[3,2]
N[1,3] <- M[3,1] - M[1,3]
N[1,4] <- M[1,2] - M[2,1]
N[2,1] <- N[1,2]
N[2,2] <- M[1,1] - M[2,2] - M[3,3]
N[2,3] <- M[1,2] + M[2,1]
N[2,4] <- M[3,1] + M[1,3]
N[3,1] <- N[1,3]
N[3,2] <- N[2,3]
N[3,3] <- -M[1,1] + M[2,2] - M[3,3]
N[3,4] <- M[2,3] + M[3,2]
N[4,1] <- N[1,4]
N[4,2] <- N[2,4]
N[4,3] <- N[3,4]
N[4,4] <- -M[1,1] - M[2,2] + M[3,3]

eigen.out <- eigen(N)

q <- Re(eigen.out[[2]][,1])

library(onion)

qh <- q[1] + Hi * q[2] + Hj * q[3] + Hk * q[4]

Xh <- Hi * X[,1] + Hj * X[,2] + Hk * X[,3]



Conj(qh) * Xh * qh

X.

my.althlooti <- function(X1,X2){
	M <- t(X2) %*% X1
	N <- matrix(0,4,4)

	N[1,1] <- M[1,1] + M[2,2] + M[3,3]
	N[1,2] <- M[2,3] - M[3,2]
	N[1,3] <- M[3,1] - M[1,3]
	N[1,4] <- M[1,2] - M[2,1]
	N[2,1] <- N[1,2]
	N[2,2] <- M[1,1] - M[2,2] - M[3,3]
	N[2,3] <- M[1,2] + M[2,1]
	N[2,4] <- M[3,1] + M[1,3]
	N[3,1] <- N[1,3]
	N[3,2] <- N[2,3]
	N[3,3] <- -M[1,1] + M[2,2] - M[3,3]
	N[3,4] <- M[2,3] + M[3,2]
	N[4,1] <- N[1,4]
	N[4,2] <- N[2,4]
	N[4,3] <- N[3,4]
	N[4,4] <- -M[1,1] - M[2,2] + M[3,3]
	
	eigen.out <- eigen(N)
	q <- Re(eigen.out[[2]][,1])
	qh <- q[1] + Hi * q[2] + Hj * q[3] + Hk * q[4]
	Xh <- Hi * X1[,1] + Hj * X1[,2] + Hk * X1[,3]

	RotX1 <- Conj(qh) * Xh * qh
	RotX1.mat <- cbind(i(RotX1),j(RotX1),k(RotX1))
	
	D0 <- sqrt(sum((X1-X2)^2))
	Dal <- sqrt(sum((RotX1.mat-X2)^2))
	return(list(X1=X1,X2=X2,M=M, N=N, eigen.out=eigen.out,q =q,RotX1 = RotX1.mat,D0=D0,Dal=Dal))
}
npt <- 1
d <- 3

n.obj <- 5
k <- 5
Obj <- list()
cnt <- 1
for(i in 1:n.obj){
	#Obj[[i]] <- matrix(rnorm(npt*d),ncol=d)
	tmp <- matrix(0,npt,d)
	tmp[1,] <- rnorm(d)
	if(npt>1){
	for(j in 2:npt){
		tmp[j,] <- tmp[j-1,] + rnorm(d)
	}
	}

	#av <- apply(tmp,2,mean)
	#Obj[[cnt]] <- t(t(tmp)-av)
	Obj[[cnt]] <- tmp
	#Obj[[cnt]] <- Obj[[cnt]]/sqrt(sum(Obj[[cnt]]^2))
	cnt <- cnt+1
	for(j in 1:k){
		R <- Random.Start(d)
		Obj[[cnt]] <- Obj[[cnt-1]] %*% R
		Obj[[cnt]] <- Obj[[cnt]] + rnorm(npt*d,0,0.0005)
		#Obj[[cnt]] <- Obj[[cnt]]/sqrt(sum(Obj[[cnt]]^2))
		cnt <- cnt+1
	}
}

nn <- length(Obj)
D0 <- Dal <- matrix(0,nn,nn)

for(i in 1:nn){
	for(j in 1:nn){
		tmp <- my.althlooti(Obj[[i]],Obj[[j]])
		D0[i,j] <- tmp$D0
		Dal[i,j] <- tmp$Dal
	}
}
range(Dal)


S <- Dal %*% solve(D0)

plot(eigen(S)[[1]])


# Sinv <- solve(S)

