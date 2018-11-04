library(MCMCpack)
library(mvtnorm)
# 多変量正規分布の分散共分散行列
Sigma <- matrix(c(3,2,2,3),2,2)
Sigma
# 変量の個数
p <- dim(Sigma)[1]
# サンプル数
n <- 10
# p次元多変量正規分布から
# n標本を取って
# n標本間の内積行列t(X)Xを作ると
# それはもとの分散共分散行列と同じ形をしているが
# その個々の成分の期待値は
# もとの分散共分散の成分のn倍で
# その個々の成分の分散は
# 元の分散共分散行列の成分を使って、n(v_{ij}^2 + v_{ii}v_{jj})となる
# 以下、確かめる

n.iter <- 1000 # nサンプルをサンプリングする試行の繰り返し回数
series.XtX <- matrix(0,n.iter,p^2) # 出来るt(X)Xの成分を格納する箱
for(i in 1:n.iter){
  X <- rmvnorm(n,sigma=Sigma)
  tmp <- t(X) %*% X
  series.XtX[i,] <- c(tmp)
}
apply(series.XtX,2,mean) # 各成分の平均
c(n * Sigma)
apply(series.XtX,2,var) # 各成分の分散
c(n * (Sigma^2 + matrix(diag(Sigma),ncol=1) %*% matrix(diag(Sigma),nrow=1)))

# p = 1、Sigma = 1のとき
# n個のサンプルをサンプリングすると
# t(X)Xは、n個の１次元標準正規乱数の距離の二乗の和
# それは自由度nのカイ二乗分布に従う
# 自由度nのカイ二乗乱数とのcoplotをすることでそれを確かめる
Sigma <- matrix(1,1,1)
# 変量の個数
p <- dim(Sigma)[1]
# サンプル数
ns <- 1:9
par(mfrow=c(3,3))
n.iter <- 1000 # nサンプルをサンプリングする試行の繰り返し回数
for(ii in ns){
  n <- ii
  series.XtX <- matrix(0,n.iter,p^2) # 出来るt(X)Xの成分を格納する箱
  for(i in 1:n.iter){
    X <- rmvnorm(n,sigma=Sigma)
    tmp <- t(X) %*% X
    series.XtX[i,] <- c(tmp)
  }
  plot(sort(series.XtX),sort(rchisq(n.iter,df=ii)),pch=20)
  abline(0,1,col=2)
}
