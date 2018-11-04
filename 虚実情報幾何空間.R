n <- 1000 # 分布の数
d <- 100 # 項数
library(MCMCpack)
P <- rdirichlet(n,runif(d)*3) 
ord <- order(P[,1])
P <- P[ord,]

# log ( 函数内積) / 2 が\theta座標内積であるので、それをペアワイズに求める
H <- log(P %*% t(P)) /2

# 固有値分解する
eigen.out <- eigen(H)
# 正定値ではない
plot(eigen.out[[1]])
# このままではうまくないが、むりやり、H = Theta %*% t(Theta) と分解することにする
# 固有ベクトルを並べた行列を複素数行列オブジェクトにする
V <- eigen.out[[2]] + 0*1i
# 固有値の平方根を対角成分にする行列
S <- diag((eigen.out[[1]] + 0*1i)^0.5)
# theta座標を各標本分布に与える
Theta <- V %*% S
# H = Theta %*% t(Theta)の検算
range(Re(Theta%*% t(Theta) -H))
range(Im(Theta%*% t(Theta) -H))

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
range(Mod(logP. - ThetaC %*% F.))

F <- F.[1:(length(F.[,1])-1),]
C <- F.[length(F.[,1]),]

# 検算２

logP.calc <- matrix(rep(C,n),nrow=n,byrow=TRUE) + Theta %*% F - Psi
range(Re(logP.calc))
range(Im(logP.calc)) # 虚部はない
# 一致する
range(Mod(logP - logP.calc))


ThetaC.real <- matrix(0,n,n+1)
ThetaC.real[,n+1] <- 1
ThetaC.real[,which(eigen.out[[1]] >= 0)] <- Re(ThetaC[,which(eigen.out[[1]] >= 0)])
ThetaC.real[,which(eigen.out[[1]] < 0)] <- Im(ThetaC[,which(eigen.out[[1]] < 0)])

matplot(ThetaC.real,type="l")

F.real <- F
F.real[which(eigen.out[[1]] >= 0),] <- Re(F.real[which(eigen.out[[1]] >= 0),])
F.real[which(eigen.out[[1]] < 0),] <- -Im(F.real[which(eigen.out[[1]] < 0),])

Theta.real <- ThetaC.real[,1:n]
logP.calc.real <- matrix(rep(C,n),nrow=n,byrow=TRUE) + Theta.real %*% F.real - Psi
range(Re(logP.calc.real))
range(Im(logP.calc.real)) # 虚部はない
# 一致する
range(Mod(logP - logP.calc.real))

# Potential functionを、固有値　正の軸は+\theta^2、
# 固有値　負の軸は -\theta^2とすることで
# すべてのF,thetaが実数になる

Psi.real <- apply(Theta.real[,which(eigen.out[[1]] >= 0)]^2,1,sum)-apply(Theta.real[,which(eigen.out[[1]] < 0)]^2,1,sum)

plot(Psi,Psi.real)

C

