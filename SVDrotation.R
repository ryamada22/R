library(svd)

# x,y,z座標球面スカラー場に対して
# k個のスペクトルが得られたとする
# それが形情報
# x,y なるm * k 行列

# 3次元ベクトルがk個
m <- 3
k <- 5

# x,y,z座標を平均0にして
# 全要素の二乗和を1に標準化する
x <- matrix(rnorm(m*k),ncol=k)
y <- matrix(rnorm(m*k),ncol=k)

x <- t(t(x)-apply(x,2,mean))
y <- t(t(y)-apply(y,2,mean))
x <- x/sqrt(sum(x^2))
y <- y/sqrt(sum(y^2))

# 検算
apply(x,1,sum)
apply(y,1,sum)
sum(x^2)
sum(y^2)



# 藤井さんの四元数法での最適回転の確認

# Althloothi
QSigma <- function(clm1, clm2, scale=FALSE){
  if(scale){
    clm1_scale <- scale(clm1, scale=FALSE)
    clm2_scale <- scale(clm2, scale=FALSE)
    dist1 <- sum((clm1_scale - clm2_scale)^2) # ???????????L2 norm
    L2 <- sum(clm1_scale^2) + sum(clm2_scale^2)
    Sigma <- t(clm1_scale) %*% clm2_scale
  } else {
    dist1 <- sum((clm1 - clm2)^2) # ???????????L2 norm
    L2 <- sum(clm1^2) + sum(clm2^2)
    Sigma <- t(clm1) %*% clm2
  }
  Aij <- Sigma - t(Sigma)
  Delta <- c(Aij[2, 3], Aij[3, 1], Aij[1, 2])
  res <- diag(0, 4)
  res[1, 1] <- sum(diag(Sigma))
  res[1, 2:4] <- Delta
  res[2:4, 1] <- Delta
  res[-1, -1] <- Sigma + t(Sigma) - res[1,1]*diag(1, 3)
  eig <- eigen(res)
  theta <- acos(eig$vector[1,1])*2
  u <- eig$vector[-1,1]/sin(theta/2)
  optdist <- L2 - 2*eig$value[1]
  return(list(mat=res, v=eig$value, e=eig$vector, theta=theta, u=u, distance=dist1, optdist=optdist))
}


#########
# 藤井さんの関数QSigmaの返り値のe の第１固有ベクトルが、回転四元数なので
# それをベクトルで返すか、四元数で返すかの関数に改変させてもらいました
my.OptimRot.QSigma <- function(x1,x2,scale=FALSE,quaternion=FALSE){
	tmp <- QSigma(x1,x2,scale=scale)
	q <- tmp$e[,1]
	if(quaternion){
		return(as.quaternion(matrix(q,ncol=1)))
	}else{
		return(q)
	}
}

# 以下のqが最適回転に対応する四元数


q <- my.OptimRot.QSigma(t(y),t(x),quaternion=TRUE)


# このqを使ってyを回転する
out <- matrix(0,4,length(y[1,]))
for(i in 1:length(y[1,])){
	out[,i] <- q * as.quaternion(matrix(c(0,y[,i]),ncol=1)) * Conj(q)
}
# 実部はゼロ
round(out,10)
# 3次元座標
out[2:4,]


# svd分解法
# http://www2.stat.duke.edu/~ab216/monograph.pdf の
# p84 の式(6.4)の辺りを使って
# yを回転して、xと近づける最適化を行う

yxt <- y %*% t(x)

svdout <- svd(yxt)
U <- svdout$u
V <- svdout$v

# svd 分解の確認
# テキストの式表現とsvd()関数の式表現のずれに注意する
U %*% diag(svdout$d) %*% t(V)
yxt

# テキストで言うところの最適回転行列 TをWとする
W <- t(t(V)) %*% t(U)

# Trace(Ty t(x)) = sum(eigenvalue)の検算

sum(diag(W %*% y %*% t(x)))
sum(svdout$d)
W %*% y


# ２つの方法の一致を確認
round(out[2:4,] - W %*% y,10)
