# 考え方
# d次元空間に一般性を持つn個の点を配置する
# その配置の仕方はd * n 次元・自由度
# ここでd次元空間での原点中心の回転を同一視したい
# d次元空間でk-1回、回転させてやるとk個の配置が作成できる
# このk個の配置をd * n 次元空間のk個の点と見ることができる
# k>=d * nとし
# これをQR分解すると
# (d * n)^2の回転行列Qと、内積を保った新しい座標を表すRとに分解できる
# このときRのランクはd^2+1になる(実験的にそうなる)
# Rのランクはkを増やしても大きくならない
# これはd次元空間での回転によってできるd * n次元の点の集合が
# d^2 + 1次元という(狭いかもしれない)空間に含まれていることを意味する
# Qは最初に使うd * n個のベクトルに依存して決まる
# Rで確かめてみよう

# Randomな回転行列を作成する関数を有するパッケージ
library(GPArotation) 

# 入力：X はn行d列の行列
# 出力：Q,RはQR分解の結果,Yは回転後のd*n座標を束ねた行列
my.dimreductionWithRotation <- function(X,k=prod(dim(X))+1){
	n <- length(X[,1])
	d <- length(X[1,])

	Y <- c(X)
	tmp <- X
	for(i in 1:k){
		r <- Random.Start(d)
		tmp <- t(r %*% t(matrix(c(tmp),ncol=d)))
		#tmp <- t(r %*% t(X))
		Y <- cbind(Y,c(tmp))
	}
	qrout <- qr(Y)
	Q <- qr.Q(qrout)
	R <- qr.R(qrout)
	return(list(Q=Q,R=R,Y=Y))
}
# number of points
n <- 20
# space dimension
d <- 2

X1 <- matrix(rnorm(n*d),ncol=d)
X1 <- X1/sqrt(sum(X1^2))
X2 <- matrix(rnorm(n*d),ncol=d)
X2 <- X2/sqrt(sum(X2^2))
qrout1 <- my.dimreductionWithRotation(X1)
# 回転の数を増やしてやる
qrout1 <- my.dimreductionWithRotation(X1,k=500)
qrout2 <- my.dimreductionWithRotation(X2,k=500)
qr(qrout1$R)$rank
qr(qrout2$R)$rank
par(mfcol=c(1,2))
image(qrout1$R)
image(qrout2$R)
par(mfcol=c(1,1))
library(rgl)
plot3d(t(qrout2$R[sample(1:qr(qrout2$R)$rank,3),]))

plot(as.data.frame(t(qrout1$R[1:qr(qrout1$R)$rank,])))
plot(as.data.frame(t(qrout2$R[1:qr(qrout2$R)$rank,])))

# ランクd^2+1とはどういうことだろうか
# d次元球面を張るのに必要な次元
# 低次元に降りているのだが、低次元球面に降りているので、もう少し次元を下げたい
# d次元で行っていることも回転だし、Qも回転なので、座標をsqrt()してやると、曲面が平面になる
# このとき、d次元回転して得られる点(の座標をsqrt()した点)の集合はd-1次元超平面であって、原点からある距離のところにある
# そのような超平面上の点集合はランクとしてd


my.dimreductionWithRotationSimplex <- function(X,k=prod(dim(X))+1){
	n <- length(X[,1])
	d <- length(X[1,])
	Y <- c(X)
	tmp <- X
	for(i in 1:k){
		r <- Random.Start(d)
		#tmp <- t(r %*% t(matrix(tmp,ncol=d)))
		tmp <- t(r %*% t(X))
		Y <- cbind(Y,c(tmp))
	}
	Y. <- sqrt(abs(Y))
	qrout <- qr(Y.)
	Q <- qr.Q(qrout)
	R <- qr.R(qrout)
	return(list(Q=Q,R=R,Y=Y.))
}


qrout2s <- my.dimreductionWithRotationSimplex(X1,k=500)
qr(qrout2s$R)$rank
d
# 得られた座標を二乗する
Q <- qrout2s$Q
R <- qrout2s$R
Rsq <- R^2
library(rgl)
plot3d(t(Rsq[1:3,]))


Q12 <- (qrout1$Q) %*% t(qrout2$Q)

tmp <- Q12 %*% c(1,rep(0,n*d-1))
tmp[1]
tmp1 <- rnorm(n*d)
tmp1 <- tmp1/sqrt(sum(tmp1^2))
tmp2 <- Q12 %*% tmp1
sum(tmp1*tmp2)

t1 <- qrout1$Q %*% tmp1
t2 <- qrout2$Q %*% tmp1

sum(t1*t2)


X <- matrix(rnorm(d*n),ncol=d)
X2 <- matrix(rnorm(d*n),ncol=d)
#X <- X/sqrt(apply(X^2,1,sum))
Y <- c(X)
Y2 <- c(X2)
Y. <- Y
Y2. <- Y2
X. <- X
X2. <- X2
n.iter <- 500
for(i in 1:n.iter){
	#tmpX <- Random.Start(d*n) %*% Y
	rr <- Random.Start(d)
	tmpX <- rr %*% t(X)
	Y. <- cbind(Y.,c(tmpX))
	tmpX <- rr %*% t(X2)
	Y2. <- cbind(Y2.,c(tmpX))
}

A <- Y.
A2 <- Y2.
# qr 分解
qrout <- qr(A)
qrout. <- qr(A[,1:61])

qrout2 <- qr(A2)
Q <- qr.Q(qrout)
Q. <- qr.Q(qrout.)

range(Q-Q.)

R <- qr.R(qrout)
Q2 <- qr.Q(qrout2)
R2 <- qr.R(qrout2)


# t(R) %*% R はAのペアワイズ距離二乗行列
t(R) %*% R
t(A) %*% A

eigen(R)[[1]] - eigen(R2)[[1]]


qr(R)$rank

D.ori <- t(Y.) %*% Y.
D.post <- t(R) %*% R

range(D.ori-D.post)

D.ori2 <- as.matrix(dist(t(Y.)))
D.post2 <- as.matrix(dist(t(R)))

range(D.ori2-D.post2)

