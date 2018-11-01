


my.theta.decomposition <- function(P){
	# log ( 函数内積) / 2 が\theta座標内積であるので、それをペアワイズに求める
	H <- log(P %*% t(P)) /2

	# 固有値分解する
	eigen.out <- eigen(H)
	# 正定値ではない
	#plot(eigen.out[[1]])
	# このままではうまくないが、むりやり、H = Theta %*% t(Theta) と分解することにする
	# 固有ベクトルを並べた行列を複素数行列オブジェクトにする
	V <- eigen.out[[2]] + 0*1i
	# 固有値の平方根を対角成分にする行列
	S <- diag((eigen.out[[1]] + 0*1i)^0.5)
	# theta座標を各標本分布に与える
	Theta <- V %*% S
	# H = Theta %*% t(Theta)の検算
	#range(Re(Theta%*% t(Theta) -H))
	#range(Im(Theta%*% t(Theta) -H))

	# 分布型標本にTheta という座標を与えることができた
	# この座標は一部の軸が実数で
	# 一部の軸が虚数である
	# 実数軸の成分は、値の差が大きいほど、分布間距離は大きく
	# 虚数軸の成分は、値の差が小さいほど、分布間距離は大きい

	# log(P) = C(x) + Theta F - \Psi(Theta)
	# のC,Fを解きたい

	logP <- log(P)
	Psi <- diag(H)

	# C(x):定数項の分をFに組み込むべく、Thetaに1列加える
	ThetaC <- cbind(Theta,rep(1,n))
	# log(P) = ThetaC %*% F. - \Psi(Theta)
	# Fの最下行がC(x)に相当する行

	# logP. = log(P) + \Psi(Theta) 
	logP. <- logP + Psi

	# logP. = ThetaC F.
	# 一般化逆行列を使って解く

	F. <- ginv(ThetaC) %*% logP.

	# 検算
	#range(Mod(logP. - ThetaC %*% F.))

	F <- F.[1:(length(F.[,1])-1),]
	Cx <- F.[length(F.[,1]),]

	# 検算２

	logP.calc <- matrix(rep(Cx,n),nrow=n,byrow=TRUE) + Theta %*% F - Psi
	#range(Re(logP.calc))
	#range(Im(logP.calc)) # 虚部はない
	# 一致する
	#range(Mod(logP - logP.calc))


	ThetaC.real <- matrix(0,n,n+1)
	ThetaC.real[,n+1] <- 1
	ThetaC.real[,which(eigen.out[[1]] >= 0)] <- Re(ThetaC[,which(eigen.out[[1]] >= 0)])
	ThetaC.real[,which(eigen.out[[1]] < 0)] <- Im(ThetaC[,which(eigen.out[[1]] < 0)])

	#matplot(ThetaC.real,type="l")

	F.real <- F
	F.real[which(eigen.out[[1]] >= 0),] <- Re(F.real[which(eigen.out[[1]] >= 0),])
	F.real[which(eigen.out[[1]] < 0),] <- -Im(F.real[which(eigen.out[[1]] < 0),])

	Theta.real <- ThetaC.real[,1:n]
	logP.calc.real <- matrix(rep(Cx,n),nrow=n,byrow=TRUE) + Theta.real %*% F.real - Psi
	#range(Re(logP.calc.real))
	#range(Im(logP.calc.real)) # 虚部はない
	# 一致する
	#range(Mod(logP - logP.calc.real))

	# Potential functionを、固有値　正の軸は+\theta^2、
	# 固有値　負の軸は -\theta^2とすることで
	# すべてのF,thetaが実数になる

	Psi.real <- apply(Theta.real[,which(eigen.out[[1]] >= 0)]^2,1,sum)-apply(Theta.real[,which(eigen.out[[1]] < 0)]^2,1,sum)

	#plot(Psi,Psi.real)

	return(list(Theta.real=Theta.real,F.real=F.real,Cx=Cx,Psi.real=Psi.real,P=P,logP=logP,H=H,eigen.out=eigen.out))

}

library(MASS)

# 単変量正規分布

my.discrete.normal <- function(m,s,bin=seq(from=-10,to=10,length=100)){
	tmp <- pnorm(bin,m,s)
	return(diff(c(0,tmp,1)))
}
n <- 500
ms <- rnorm(n)
#ms <- rep(0,n)
ss <- runif(n)*3
#ss <- rep(1,n)
bin <- seq(from=-2,to=2,length=500)
P <- matrix(0,n,length(bin)+1)

for(i in 1:n){
	P[i,] <- my.discrete.normal(ms[i],ss[i],bin=bin)
}
# 正規化されていないＰの対策をする
# 確率質量が0の要素の対策をする
my.st.P <- function(P){
	ret <- P/apply(P,1,sum)
	ret <- ret + 10^(-15)
	return(ret/apply(ret,1,sum))
}

out1 <- my.theta.decomposition(my.st.P(P))
library(rgl)

plot3d(out1$Theta.real[,1],out1$Theta.real[,n],ss)
plot3d(out1$Theta.real[,1],ms,ss)
plot3d(out1$Theta.real[,n],ms,ss)

# ポアッソン分布
n <- 100
lambdas <- runif(n)*10

N <- 20
P <- matrix(0,n,N+2)

for(i in 1:n){
	P[i,1:(N+1)] <- dpois(0:N,lambdas[i])
	P[i,N+2] <- 1-sum(P[i,])
}


out2 <- my.theta.decomposition(my.st.P(P))

plot3d(out2$Theta.real[,1],out2$Theta.real[,n],lambdas)

# 50:50混合ポアッソン分布

n <- 800
lambdas1 <- runif(n) * 10
lambdas2 <- runif(n) * 5

N <- 20
P <- matrix(0,n,N+2)

for(i in 1:n){
	P[i,1:(N+1)] <- dpois(0:N,lambdas1[i])*0.5 + dpois(0:N,lambdas2[i])*0.5
	P[i,N+2] <- 1-sum(P[i,])
}


out3 <- my.theta.decomposition(my.st.P(P))

plot3d(out3$Theta.real[,1],lambdas1,lambdas2)
plot3d(out3$Theta.real[,n],lambdas1,lambdas2)

plot3d(out3$Theta.real[,1],out3$Theta.real[,n],lambdas1)
spheres3d(out3$Theta.real[,1],out3$Theta.real[,n],lambdas1,color=gray(lambdas2/max(lambdas2)),radius=0.15)

