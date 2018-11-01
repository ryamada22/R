library(devtools)
library(mwaytable)

library(MCMCpack)

# 適当にnxm行列型確率テーブルを作る
n <- 2
m <- 3
# 周辺分布
margn <- rdirichlet(1,rep(1,n))
margm <- rdirichlet(1,rep(1,m))

# 期待値テーブル
e.table <- matrix(c(outer(margn,margm,"*")),n,m)
# 検算
apply(e.table,1,sum)
apply(e.table,2,sum)
margn
margm

# 観測テーブルを作る
o.table <- e.table
o.table[1,1] <- o.table[1,1] * runif(1)
o.table[1,2] <- o.table[1,2] * runif(1)
o.table[1,3] <- margn[1]-o.table[1,1]-o.table[1,2]
o.table[2,] <- margm - o.table[1,]

# バサクの回転行列
r <- c(n,m)
KM <- make.simplex.multi(r)
dim(KM)

# 末尾に１を加えたテーブルベクトル
e.vec <- c(e.table,1)
o.vec <- c(o.table,1)

# ベクトルの長さ。１長いものとテーブルの要素数と。
N1 <- length(e.vec)
N <- N1-1

# 回転せずに、期待値テーブルを原点にするための行列
#Q : Shifting matrix

Q <- diag(N1)
Q[1:N,N1] <- -e.vec[-N1]
# 確かに、e.vecは(0,...,0,1)に移っている
Q %*% o.vec
Q %*% e.vec

# バサクの回転行列をN1 x N1行列に拡大する
R <- diag(N1)
R[1:N,1:N] <- KM
# 期待テーブルが回転した先の座標を出す
fe <- R %*% e.vec
# バサクの回転行列を拡大したものの、右端列を、期待テーブル座標を引くために使う
R[1:N,N1] <- -fe[-N1]

# このRで座標変換する
R %*% o.vec # 0がたくさんある
R %*% e.vec # 原点

# やっぱ、dfに応じて行を選ばないとだめ！
Zeroing <- diag(c(1,0,1,0,0,0,0))
# t(o-e)t(Q) E^{-1} Q (o-e) = t(o-e) t(R) t(Zeroing) t(X) X Zeroing R (o-e)を解く
XtX <- Zeroing %*% solve(t(R)) %*% t(Q) %*% diag(1/e.vec) %*% Q %*% solve(R) %*% Zeroing

eigenout <- eigen(XtX)

X <- diag(sqrt(eigenout[[1]])) %*% t(eigenout[[2]])

# これが、テーブルベクトルに1を継ぎ足したベクトルから
# 期待値テーブルを原点に移しつつ、観察テーブルを自由度次元座標で表し
# かつ、距離^2がchisqになるような行列
# 第7列が『期待値を原点にするために必要』
W <- X %*% R

# 回してみる
We <- W %*% e.vec
Wo <- W %*% o.vec

# ちゃんと自由度次元
round(We,5)
round(Wo,5)

# chi-sq計算
sum((o.table-e.table)^2/e.table)
# t(v)t(Q)E^{-1} Q vを使って計算

t(matrix(o.vec,ncol=1)) %*% t(Q) %*% diag(1/e.vec) %*% Q %*% o.vec -1 
# 作った変換行列で内積がchisqになっている
sum(Wo^2)


