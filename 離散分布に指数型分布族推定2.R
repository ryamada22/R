# 離散分布(p,1-p)に指数型分布族推定を試してみた！
# (1,0)と(0,1)とは無限に遠いことがわかる
# こんな例を載せるとよいかも知れません

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
dd <- 2
H - svdout$v[,1:dd] %*% diag(svdout$d[1:dd]) %*% t(svdout$u[,1:dd])

svdout$v[,1:dd] -svdout$u[,1:dd]




cholout <- chol(H)

out <- eigen(H)

V <- out[[2]]
eig <- out[[1]]
eig <- abs(eig) # 実質0の固有値に負がある

X <- V %*% diag(sqrt(eig))

H - X%*%t(X)


M <- matrix(0,n,n)

for(i in 1:n){
	for(j in 1:n){
		if(i==j){
			M[i,j] <- sum(R[i,]^2)
			M[i,j] <- 0
		}else{
			M[i,j] <- sum(R[i,]^2) + sum(R[j,]^2) - log(sum(R[i,]*R[j,]))
		}
	}
}


out <- cmdscale(M,k=n-1)


plot(R[,1],out[,1],xlab="p[1]",ylab="theta_1")
plot(apply(out,2,var))

plot(as.data.frame(out[,1:5]))

# こういうわかりやすい分布を使いつつ
# 均等に分布を発生させるのではなく、偏りを持たせると、よい例ができそうです

n <- 100 # 分布の数
n1 <- 40
n2 <- n-n1
d <- 2 # 項数
library(MCMCpack)
R1 <- rdirichlet(n1,c(1,2))
R2 <- rdirichlet(n2,c(0.4,3))
R <- rbind(R1,R2)

M <- matrix(0,n,n)

for(i in 1:n){
	for(j in 1:n){
		if(i==j){
			M[i,j] <- sum(R[i,]^2)
			M[i,j] <- 0
		}else{
			M[i,j] <- sum(R[i,]^2) + sum(R[j,]^2) - log(sum(R[i,]*R[j,]))
		}
	}
}


out <- cmdscale(M,k=n-1)


plot(R[,1],out[,1],xlab="p[1]",ylab="theta_1")
plot(apply(out,2,var))

plot(as.data.frame(out[,1:5]))
