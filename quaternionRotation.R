library(onion)

theta <- runif(1) * 2 * pi
q.v <- rnorm(3)
q.v <- q.v/sqrt(sum(q.v^2))


q <- cos(theta/2) + sin(theta/2) * (q.v[1]*Hi + q.v[2]*Hj + q.v[3]*Hk)
q. <- Conj(q)


# quaternion‚Ì‚S‚˜‚Ss—ñ•\Œ»‚Í•¡”‚ ‚é
# ŽÀÛ‚S‚WŽí—Þ‚ ‚é(˜a‚ÆÏ‚ð–ž‘«‚·‚é‚à‚Ì‚ª)
# ‚»‚Ì‚¤‚¿A‹¤–ðquaternion‚ª“]’us—ñ‚Å‚ ‚Á‚Ä
# norm‚Ì‚Sæ‚ªdeterminant‚É‚È‚é‚à‚Ì‚ª‚QŽí—Þ‚ ‚é
# ‚»‚ê‚ªˆÈ‰º‚Ì‚Q‚Â

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

# ‘«‚µ‡‚í‚¹
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

wH <- pH




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

