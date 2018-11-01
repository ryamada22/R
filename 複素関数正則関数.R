# zは複素数のベクトル
# int0,int1はintensityの上下限、sat0,sat1はSaturation(彩度)の上下限
my.hsv <- function(z,int0=0.6,sat0=0.3,int1=1,sat1=1){
# 複素数の偏角
	arg <- Arg(z)
	s <- which(arg<0)
	arg[s] <- arg[s]+2*pi
# 複素数の絶対値
	r <- Mod(z)
# 絶対値が非常に大きくてもそこそこの色になるように対数変換
	s <- which(r>1)
	r[s] <- log(r[s])
# 絶対値で周期性が出るように4のmod
	r. <- 4*(r%%1)
	k <- floor(r.)
	r. <- r.-k
# 明度が上限、明度が下限、彩度が上限、彩度が下限の４パターンを
# 4のmodに対応づける
# 明度・彩度を動かすときは、複素数の絶対値で１次線形変化
	inten <- sat <- rep(0,length(r))
	s <- which(k==0)
	inten[s] <- int1
	sat[s] <- sat1-(sat1-sat0)*r.[s]
	s <- which(k==1)
	inten[s] <- int1-(int1-int0)*r.[s]
	sat[s] <- sat0
	s <- which(k==2)
	inten[s] <- int0
	sat[s] <- sat1-(sat1-sat0)*(1-r.[s])
	s <- which(k==3)
	inten[s] <- int1-(int1-int0)*(1-r.[s])
	sat[s] <- sat1

	return(cbind(arg,inten,sat))
}

my.hsv2rgb <- function(h,s,v){
# 色相の6 のmodでぐるりの情報を作る
	hi <- floor(h/(2*pi)*6)
	hi[which(hi==6)] <- 0
# 色相のぐるりの余りをfに入れ、それと明度・彩度とでp,q,tという３変数を決める
# ３変数を色相からの値を取らせる１つの原色を除いた２原色の値を定めるために使う
# 使い方は巡回させることでうまいことやる
	f <- (h/(2*pi)*6) %%1
	p <- v*(1-s)
	q <- v *(1-f*s)
	t <- v *(1-(1-f)*s)
	r <- g <- b <- rep(0,length(h))
	s <- which(hi==0)
		r[s] <- v[s];g[s] <- t[s]; b[s] = p[s];
	s <- which(hi==1)
		r[s] <- q[s];g[s] <- v[s]; b[s] = p[s];
	s <- which(hi==2)
		r[s] <- p[s];g[s] <- v[s]; b[s] = t[s];
	s <- which(hi==3)

		r[s] <- p[s];g[s] <- q[s]; b[s] = v[s];
	s <- which(hi==4)
		r[s] <- t[s];g[s] <- p[s]; b[s] = v[s];
	s <- which(hi==5)
		r[s] <- v[s];g[s] <- p[s]; b[s] = q[s];
	return(cbind(r,g,b))
}

my.z.color <- function(z){
	a <- Arg(z)
	col <- a %/% (pi/2) + 3
	return(col)
}
x <- seq(from=-10,to=10,len=100)
xx <- expand.grid(x,x)
z <- xx[,1]+1i * xx[,2]

#my.f <- function(z){
#	(z^2-1)*(z-2-1i)^2/((1+1i)*z^4+z^2+2+2*1i)
#}
my.f <- function(z){
	(z^2-1)*(z-2-1i)^2/((1+1i)*z^4+z^2+1+1i)
}

w <- my.f(z)
hsv <- my.hsv(w,int0=0.1,sat0=0.1,int1=1,sat1=1)
col <- my.hsv2rgb(hsv[,1],hsv[,3],hsv[,2])
col2 <- rgb(col[,1],col[,2],col[,3])
plot(xx,pch=20,col=col2)

n.pt <- 1000
X <- matrix(rnorm(n.pt*3),ncol=3)
X <- X/sqrt(apply(X^2,1,sum))


library(rgl)
plot3d(X)

my.RiemannSphere <- function(z){
	a <- Arg(z)
	r <- Mod(z)
	theta <- atan(2/r)
	x3 <- cos(pi-2*theta)
	x1 <- sin(pi-2*theta) * cos(a)
	x2 <- sin(pi-2*theta) * sin(a)
	return(cbind(x1,x2,x3))
}

z.Riemann <- my.RiemannSphere(z)

plot3d(z.Riemann)
spheres3d(z.Riemann,radius=0.07,color=col2)

plot3d(z.Riemann)
spheres3d(z.Riemann,radius=0.07,color=my.z.color(w))

plot(xx,pch=20,col=my.z.color(w))

plot3d(xx[,1],xx[,2],Mod(w))
