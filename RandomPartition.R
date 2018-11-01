

n.iter <- 10^4
K <- 10^3
X <- rep(0,n.iter)
Ns <- X
N <- 10
for(i in 1:n.iter){
	N <- sample(5:15,1)
	Ns[i] <- N
	R <- runif(N)
	R <- R/sum(R)

	tmp <- sample(1:N,K,replace=TRUE,prob=R)
	X[i] <- R[tmp[1]]
}

Nrange <- range(Ns)
meanX <- rep(0,Nrange[2]-Nrange[1]+1)
for(i in Nrange[1]:Nrange[2]){
	tmp2 <- which(Ns==i)
	meanX[i-Nrange[1]+1] <- mean(1/X[tmp2])
}

plot((Nrange[1]:Nrange[2]),meanX)



plot(Ns,1/X)



mean(X)
1/mean(X)
N
