library(GPArotation)
d <- 5
I <- diag(rep(1,d))
l <- rep(-1,d)*1
u <- rep(1,d)*1
n.iter <- 1000
tn <- vol <- rep(0,n.iter)
for(i in 1:n.iter){
	R <- Random.Start(d)
	m <- (runif(d)) * 2
	P <- diag(1/m) %*% R
	A <- diag(m) %*% R
	A. <- cbind(rep(0,d),A)
	A.. <- rbind(A.,rep(1,d+1))
	vol[i] <- abs(det(A..))
	
	#P <- P/sqrt(apply(P^2,2,sum))
	tn[i] <- mvNcdf(l,u,P %*% t(P),1000)[[1]]
}


plot(tn,vol)

############

d <- 3
chi <- 1
l <- u <- rep(1,d) * chi
l <- l*(-1)

n.iter <- 1000
y <- vol <- rep(0,n.iter)
for(i in 1:n.iter){
	A <- matrix(rnorm(d^2),ncol=d)
	A <- A/sqrt(apply(A^2,1,sum))
	Sigma <- A %*% t(A)
	y[i] <- mvNcdf(l,u,Sigma,1000)[[1]]
	
	Ainv <- solve(A)
	
	V0 <- Ainv %*% l
	tmp <- matrix(rep(l,d),ncol=d)
	diag(tmp) <- u
	V <- Ainv %*% tmp
	K <- cbind(V0,V)
	K <- rbind(K,rep(1,d+1))
	vol[i] <- abs(det(K))
}

plot(y,vol)

###########
d <- 100
chi <- 1
l <- u <- rep(1,d) * chi
l <- l*(-1)

n.iter <- 100
k <- seq(from=0,to=5,length=n.iter+1)
k <- k[-1]
y <- vol <- rep(0,n.iter)
	A <- matrix(rnorm(d^2),ncol=d)
	A <- A/sqrt(apply(A^2,1,sum))
	Sigma <- A %*% t(A)

for(i in 1:n.iter){
	l. <- l*k[i]
	u. <- u*k[i]
	y[i] <- mvNcdf(l.,u.,Sigma,10)[[1]]
	
	Ainv <- solve(A)
	
	V0 <- Ainv %*% l.
	tmp <- matrix(rep(l.,d),ncol=d)
	diag(tmp) <- u.
	V <- Ainv %*% tmp
	K <- cbind(V0,V)
	K <- rbind(K,rep(1,d+1))
	vol[i] <- abs(det(K))
}

plot(y,vol)
plot(k,y)

