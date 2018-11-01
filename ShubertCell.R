library(gtools)
library(complexplus)

n <- 7
k <- 3
cmb <- combinations(n,k)

my.shubertcell <- function(M){
	k <- length(M[,1])
	n <- length(M[1,])
	cmb <- combinations(n,k)
	for(i in length(cmb[,1]):1){
		d <- Det(M[,cmb[i,]])
		if(d != 0){
			U <- solve(M[,cmb[i,]]) %*% M
			break
		}
	}
	return(list(U=U,position=cmb[i,]))
}

n.iter <- 10

for(i in 1:n.iter){
	# canonical formのposition列を適当に選ぶ
	s <- cmb[sample(1:length(cmb[,1]),1),]
	M <- matrix(sample((1):9,k*n,replace=TRUE),k,n) + 1i * matrix(sample((1):9,k*n,replace=TRUE),k,n)
	# position 列を単位行列にする
	M[,s] <- diag(rep(1,k))
	# position列から右側を0にして canonical formにする
	for(j in 1:k){
		if(s[j]<n){
			M[j,(s[j]+1):n] <- 0
		}
	}
	# ポジション列
	print("s")
	print(s)
	# canonical formの確認
	print("M")
	print(M)
	# 線形変換
	R <- matrix(rnorm(k^2)+rnorm(k^2)*1i,k,k)
	M. <- R %*% M
	# さらに行をシャッフル(これも一種のkxk行列による線形変換だから大丈夫に決まっているが)
	M. <- M.[sample(1:k),]
	out <- my.shubertcell(M.)
	print("M.")
	print(M.)
	# 再計算したcanonical form
	print("U")
	print(round(out$U,5))
	print(round(M - out$U,5))
}


#M <- cbind(diag(rep(1,k)),matrix(sample((1:9),k*(n-k),replace=TRUE,prob=c(rep(1,9),10,rep(1,9))),k,n-k))



for(i in 1:n.iter){
	s <- cmb[sample(1:length(cmb[,1]),1),]
	M <- matrix(sample((1):9,k*n,replace=TRUE),k,n) + 1i * matrix(sample((1):9,k*n,replace=TRUE),k,n)
	M[,s] <- diag(rep(1,k))
	for(j in 1:k){
		if(s[j]<n){
			M[j,(s[j]+1):n] <- 0
		}
	}
	print("s")
	print(s)
	print("M")
	print(M)
	#R <- matrix(rnorm(n^2),n,n)
	#M. <- t(R %*% t(M))
	R <- matrix(rnorm(k^2),k,k)
	M. <- R %*% M
	print("M.")
	print(M.)
	for(ii in length(cmb[,1]):1){
		d <- Det(M.[,cmb[ii,]])
		if(d != 0){
			U <- solve(M.[,cmb[ii,]]) %*% M.
			print("U")
			print(round(U,5))
			print("M-U")
			print(round(M-U,5))
			break
		}
	}
}

for(i in 1:length(cmb[,1])){
	d <- Det(M[,cmb[i,]])
	if(d != 0){
		U <- solve(M[,cmb[i,]]) %*% M
		print(U)
	}
}

##########
n <- 3
V <- matrix(rnorm(n^2),n,n)
det(V)
Vinv <- solve(V)
Vinv %*% V
v <- V[,1] * rnorm(1)

Vinv %*% v

v <- V[,1] * rnorm(1) + V[,2] * rnorm(1)
Vinv %*% v

v <-V[,1] * rnorm(1) + V[,2] * rnorm(1) + V[,3] * rnorm(1)

Vinv %*% v

