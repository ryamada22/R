# 内積行列を作って、svd分解する
# 分布ペアの内積を計算している
# 内積を要素とする行列ができる
# 指数型分布表現のtheta座標のベクトル間のペアワイズ内積を取ったときに
# 分布ペア内積行列となってくれたらよい

# やってみる
n <- 10 # 分布の数
d <- 2 # 項数
library(MCMCpack)
R <- rdirichlet(n,rep(0.6,d)) 

H <- R %*% t(R)


H <- matrix(0,n,n)

for(i in 1:n){
	for(j in 1:n){
		H[i,j] <- (sum(R[i,]*R[j,]))
	}
}

svdout <- svd(H)
# svd分解のvとuは本当は違うが、この場合は意味のある部分では一緒

plot(svdout$d)
# 分解を確認
range(H - svdout$v %*% diag(svdout$d) %*% t(svdout$u))
# 実は、svdout$vとsvdout$uは実質的な部分は同じ行列になっている

# 内積行列 = Q %*% t(Q)にして、Qを座標として取り出したい
Q <- svdout$v %*% diag(sqrt(svdout$d))

range(Q %*% t(Q) - H)
# 元のRが再現できる…(符号の反転とかはあるらしい)

plot(R[,1],Q[,1])

plot(R[,2],Q[,2])


# さて。問題は、これが、正規分布とか、任意の分布でできるだろうか、というのが問題か・・・