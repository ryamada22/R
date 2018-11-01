# おまじない・・・
# 
exp.m <- function(A,n){
	# 固有値分解
	eigen.out<-eigen(A)
	# P=V,P^{-1}=U
	V<-eigen.out[[2]]
	U<-solve(V)
	B<-diag(exp(eigen.out[[1]]*n))
	X <- V%*%B%*%U
	return(Re(X))
}

# 初期状態
x0 <- c(1,1)
# "線形代数からの準備" 『演習問題』1.
A <- matrix(c(3,-2,5,-2),2,2)
A

t <- seq(from=0,to=10,length=1000)

X <- matrix(0,length(t),2)

for(i in 1:length(t)){
	X[i,] <- exp.m(A,t[i]) %*% x0
}
par(mfcol=c(1,2))

plot(X,type="l")
matplot(X,type="l")
eigen(A)

A <- matrix(c(3,5,-2,-2),2,2)
t <- seq(from=0,to=10,length=1000)

X <- matrix(0,length(t),2)

for(i in 1:length(t)){
	X[i,] <- exp.m(A,t[i]) %*% x0
}

plot(X,type="l")
matplot(X,type="l")
eigen(A)


A <- matrix(c(-3,-2,5,3),2,2)
t <- seq(from=0,to=10,length=1000)

X <- matrix(0,length(t),2)

for(i in 1:length(t)){
	X[i,] <- exp.m(A,t[i]) %*% x0
}

plot(X,type="l")
matplot(X,type="l")
eigen(A)

A <- matrix(c(1,-1,1,3),2,2)
A

t <- seq(from=0,to=10,length=1000)

X <- matrix(0,length(t),2)

for(i in 1:length(t)){
	X[i,] <- exp.m(A,t[i]) %*% x0
}
plot(X,type="l")
matplot(X,type="l")
eigen(A)

