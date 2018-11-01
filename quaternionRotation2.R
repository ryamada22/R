library(onion)
# quaternionの４ｘ４行列表現は複数ある
# 実際４８種類ある(和と積を満足するものが)
# そのうち、共役quaternionが転置行列であって
# normの４乗がdeterminantになるものが２種類ある
# それが以下の２つ

my.H2mat1 <- function(q){
	r <- Re(q)
	x <- i(q)
	y <- j(q)
	z <- k(q)
	ret <- matrix(c(r,x,y,z,-x,r,-z,y,-y,z,r,-x,-z,-y,x,r),4,4)
	return(ret)
}

my.H2mat2 <- function(q){
	r <- Re(q)
	x <- i(q)
	y <- j(q)
	z <- k(q)
	ret <- matrix(c(r,x,y,z,-x,r,z,-y,-y,-z,r,x,-z,y,-x,r),4,4)
	return(ret)
}

# m=3次元、k個の点の座標行列を二つとって、その最近接回転を出す

my.OptimRot <- function(x1,x2){
	n <- length(x1[,1])
	x1 <- x1/sqrt(sum(x1^2))
	x2 <- x2/sqrt(sum(x2^2))
	x1H <- x1[,1] * Hi + x1[,2] * Hj + x1[,3] * Hk
	x2H <- x2[,1] * Hi + x2[,2] * Hj + x2[,3] * Hk
	retM <- matrix(0,4,4)
	for(i in 1:n){
		pH <- x1H[i]
		wH <- x2H[i]
		P <- my.H2mat1(pH)
		W <- my.H2mat2(wH)
		M <- t(P) %*% W
		retM <- retM + M
	}
	eigen.out <- eigen(retM)
	return(eigen.out)
}

n <- 3

m <- 3
k <- 5

X <- list()
for(i in 1:n){
	X[[i]] <- matrix(rnorm(m*k),ncol=m)
	X[[i]] <- X[[i]]/sqrt(sum(X[[i]]^2))
}

optim.rot <- array(0,c(n,n,4))

for(i in 1:n){
	x1 <- X[[i]]
	for(j in 1:n){
		x2 <- X[[j]]
		tmp <- my.OptimRot(x1,x2)
		optim.rot[i,j,] <- tmp[[2]][,1]
	}
}



theta <- runif(1) * 2 * pi
q.v <- rnorm(3)
q.v <- q.v/sqrt(sum(q.v^2))


q <- cos(theta/2) + sin(theta/2) * (q.v[1]*Hi + q.v[2]*Hj + q.v[3]*Hk)
q. <- Conj(q)



p <- rnorm(3)
p <- p/sqrt(sum(p^2))
w <- rnorm(3)
w <- w/sqrt(sum(w^2))

pH <- p[1] * Hi + p[2] * Hj + p[3] * Hk
wH <- w[1] * Hi + w[2] * Hj + w[3] * Hk





P <- my.H2mat1(pH)
W <- my.H2mat2(wH)

M <- t(P) %*% W


t(as.matrix(q * pH * q.)) %*% as.matrix(wH)

t(as.matrix(q * pH)) %*% as.matrix(wH * q)
t(as.matrix(q)) %*% M %*% as.matrix(q)
#q * pH * q. * wH

det(M)
round(M %*% t(M),5)

sum(diag(M))

eigen(M)
pH * wH

my.H2mat1(pH) %*% my.H2mat1(wH)

my.H2mat2(pH) %*% my.H2mat2(wH)

########
my.H2mat2(pH) + t(my.H2mat1(pH))

pH * 2

# 足し合わせ
n <- 10
library(MCMCpack)

Xp <- rdirichlet(1,rep(1,n))
Xw <- rdirichlet(1,rep(1,n))


Ms <- list()
for(i in 1:n){
p <- rnorm(3)
p <- p/sqrt(sum(p^2))
w <- rnorm(3)
w <- w/sqrt(sum(w^2))
p <- p * sqrt(Xp[i])
w <- w * sqrt(Xp[i])

pH <- p[1] * Hi + p[2] * Hj + p[3] * Hk
wH <- w[1] * Hi + w[2] * Hj + w[3] * Hk






P <- my.H2mat1(pH)
W <- my.H2mat2(wH)

M <- t(P) %*% W
Ms[[i]] <- M

}

sumM <- matrix(0,4,4)
for(i in 1:n){
	sumM <- sumM + Ms[[i]]
}

det(sumM)
round(sumM %*% t(sumM),5)

sum(diag(sumM))

eigen(sumM)

