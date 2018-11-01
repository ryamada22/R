n <- 4
a <- runif(n)
b <- runif(n)

a <- a/sum(a)
b <- b/sum(b)

a <- rep(1/n,n)
b <- rep(1/n,n)

M <- matrix(runif(n^2),n,n)


lambda <- 0.1

K <- exp(-lambda*M)
ainvK <- diag(1/a) %*% K

u <- rep(1/n,n)
n.iter <- 1000

u.hx <- matrix(0,n.iter+1,n)
u.hx[1,] <- u
for(i in 1:n.iter){
	tmp2 <- t(K) %*% u
	tmp1 <- ainvK %*% (c(b)/c(tmp2))
	tmp <- 1/tmp1
	u.hx[i+1,] <- tmp
	u <- tmp
}

matplot(u.hx,type="l")

alpha <- 1/lambda * log(u) - sum(log(u))/(lambda*n)

alpha


P <- diag(c(alpha)) %*% M %*% diag(c(alpha))

