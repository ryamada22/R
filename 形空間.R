my.standardSh <- function(Sh){
  M <- t(Sh[1:3,])
  qrout <- qr(M)
  Q <- qr.Q(qrout)
  
  Sh_s <- t(t(Q) %*% t(Sh))
  return(list(Sh_s=Sh_s,P=Q))
}

####
# R“c”Å
####
library(onion)
# X1‚ğ‰ñ‚µ‚ÄX2‚É‹ß‚Ã‚¯‚é
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
	IP0 <- sum(X1*X2)
	IPal <- sum(RotX1.mat * X2)
	return(list(X1=X1,X2=X2,M=M, N=N, eigen.out=eigen.out,q =qh,qmat=my.q2rotmat(qh),RotX1 = RotX1.mat,D0=D0,Dal=Dal,IP0=IP0,IPal=IPal))
}

# lŒ³”‚©‚ç‘Î‰‚·‚é‚R‚˜‚R‰ñ“]s—ñ‚É•ÏŠ·
my.q2rotmat <- function(q){
  x <- Re(q)
  y <- i(q)
  z <- j(q)
  w <- k(q)
  R <- matrix(c(x^2+y^2-z^2-w^2,2*(y*z-x*w),2*(x*z+y*w),
                2*(x*w+y*z),x^2-y^2+z^2-w^2,2*(-x*y+z*w),
                2*(y*w-x*z),2*(z*w+x*y),x^2-y^2-z^2+w^2),
              byrow=TRUE,3,3)
  return(t(R))
  
}
# ‚R‚c‰ñ“]²’PˆÊƒxƒNƒgƒ‹‚Æ‰ñ“]Šp‚ğw’è‚µ‚ÄA‰ñ“]lŒ³”‚ğì‚é
my.rotq <- function(v,theta){
	v <- v/sqrt(sum(v^2))
	cos(theta/2) + sin(theta/2) * (Hi*v[1]+Hj*v[2]+Hk*v[3])
}
# ’PˆÊ’´‹…–Êã‚Å‚ÌüŒ`˜a
# x,y‚Í’PˆÊƒxƒNƒgƒ‹
my.vector.sum.sp <- function(x,y,p){
  ip <- sum(x*y)
  theta <- acos(ip)
  newangle <- p * theta
  cmp <- cos(newangle) + 1i * sin(newangle)
  arg <- Arg(cmp)
  perpen.v <- y-ip*x
  perpen.v.st <- perpen.v/sqrt(sum(perpen.v^2))
  #new.v.direction <- ip*x + perpen.v*tan(arg)/tan(theta)
  new.v.direction <- cos(arg) * x + sin(arg) * perpen.v.st
  #print(sum(new.v.direction^2))
  z <- new.v.direction/sqrt(sum(new.v.direction^2))
  return(list(z=z,V1=x,V2=perpen.v.st,comp1=cos(arg),comp2=sin(arg),angle=theta,newangle=newangle,ip=ip))
}
my.vector.sum.sp2 <- function(x,y,p){
  ip <- sum(x*y)
  theta <- acos(ip)
  #newangle <- p * theta
  newangle <- p
  cmp <- cos(newangle) + 1i * sin(newangle)
  arg <- Arg(cmp)
  perpen.v <- y-ip*x
  perpen.v.st <- perpen.v/sqrt(sum(perpen.v^2))
  #new.v.direction <- ip*x + perpen.v*tan(arg)/tan(theta)
  new.v.direction <- cos(arg) * x + sin(arg) * perpen.v.st
  #print(sum(new.v.direction^2))
  z <- new.v.direction/sqrt(sum(new.v.direction^2))
  return(list(z=z,V1=x,V2=perpen.v.st,comp1=cos(arg),comp2=sin(arg),angle=theta,newangle=newangle,ip=ip))
}
###########
library(GPArotation)
library(onion)
k <- 5
d <- 3
X <- matrix(rnorm(k*d),ncol=d)
Y <- matrix(rnorm(k*d),ncol=d)
X <- X/sqrt(sum(X^2))
Y <- Y/sqrt(sum(Y^2))

al.out <- my.althlooti(Y,X)
Y.rot <- al.out$RotX1

rot.vec <- rnorm(d)
rot.vec <- rot.vec/sqrt(sum(rot.vec^2))
# theta <- runif(1) * 2*pi # ”CˆÓ‚ÌŠp“x‚Å‘Î‰‚ªæ‚ê‚é‚±‚Æ‚ÍŠm”FÏ‚İ
theta <- pi
rot.q <- my.rotq(rot.vec,theta)
R <- my.q2rotmat(rot.q)

X. <- t(R %*% t(X))
Y.rot. <- t(R %*% t(Y.rot))
Y.rot.2 <- my.althlooti(Y.rot.,X.)$RotX1
Y.rot.22 <- my.althlooti(Y,X.)$RotX1

range(Y.rot.-Y.rot.2)
range(Y.rot.-Y.rot.22)

# 3d‰ñ“]‚Å‚Í‘å‰~‚Íì‚Á‚½‚ç‚¾‚ß(”½“]‚³‚¹‚é‚±‚Æ‚¾‚©‚ç)
# X.. <- my.vector.sum.sp2(X,X.,pi)$z
# Y.rot.. <- my.vector.sum.sp2(Y.rot,Y.rot.,pi)$z
# Y.rot..2 <- my.althlooti(Y.rot..,X..)$RotX1
# Y.rot..22 <- my.althlooti(Y,X..)$RotX1

# range(Y.rot..-Y.rot..2)
# range(Y.rot..-Y.rot..22)
