
m <- 2
k <- 3

z <- matrix(rnorm(m*k),nrow=m)
z[1,] <- z[1,]*0.3 + 1
z[2,] <- z[2,]*0.2 + 2

ang <- Arg(z[1,] + i * z[2,])
z <- z[,order(ang)]


plot(t(z),pch=20,asp=1)

z.normed <- z - apply(z,1,mean)
plot(t(z.normed),pch=20,asp=1)

z.normed. <- z.normed/sqrt(sum(z.normed^2))

plot(t(cbind(z,z.normed,z.normed.)),pch=20,asp=1,type="p")
points(t(z),pch=20,col=1,type="b")
points(t(z.normed),pch=20,col=2,type="b")
points(t(z.normed.),pch=20,col=3,type="b")


library(MASS)
library(Matrix)
library(rgl)

n <- 1000 # データセットの数
# そのmxn次元座標を納める行列
Z <- matrix(0,n,m*k)

for(i in 1:n){
	z <- matrix(rnorm(m*k),nrow=m)
	z[1,] <- z[1,]*0.3 + 1
	z[2,] <- z[2,]*0.2 + 2
	z.normed <- z - apply(z,1,mean)
	z.normed. <- z.normed/sqrt(sum(z.normed^2))
	Z[i,] <- z.normed.
}

# ランクがmxk-m
# mxk-m次元空間に収まっていることがわかる
rankMatrix(Z)
# QR分解を用いて、mxk-m次元空間に移す回転行列を計算する
# 発生させたmxkベクトルのうち、ランク本のベクトルを用いる
df <- m*k - (m)
topZ <- Z[1:df,]
# QR分解
qr.out <- qr(t(topZ))
Q <- qr.Q(qr.out)
R <- qr.R(qr.out)

# QR分解により、t(topZ)のペアワイズなベクトル内積が
# R(減次元ベクトルの束)のペアワイズなベクトル内積と一致することを確認する
# 以下はessentially zero
range(t(R) %*% R - topZ %*% t(topZ))

# t(topZ) = QR と出来たから
# Q.inv t(topZ) = R のQ.invを作る
Q.inv <- ginv(Q)
# 検算
round(Q.inv %*% t(topZ) - R)
# Q.invを使って、すべてのmxkベクトルを減次元写像する
rotZ <- Q.inv %*% t(Z)
# ペアワイズ内積が保たれていることを確認する
d1 <- Z %*% t(Z)
d2 <- t(rotZ) %*% rotZ
range(d1-d2)

# 減次元座標のノルムの確認
range(apply(rotZ^2,2,sum))

# mxk -m = 4次元空間の球面はなので、３次元空間にプロットすると
# 球面と球の内部に点が分布する
plot3d(t(rotZ[1:3,]))
plot3d(t(rotZ[2:4,]))
