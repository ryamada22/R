# クラスごとの標本数
n1 <- 10
n2 <- 15
n3 <- 18

# 変数の数
d <- 4
# クラスごとの標本データ
X1 <- matrix(rnorm(n1*d),ncol=d)
X2 <- matrix(rnorm(n2*d),ncol=d)
X3 <- matrix(rnorm(n3*d),ncol=d)

Xall <- rbind(X1,X2,X3)
# クラスごと、全体の変数ごとの平均
m1 <- apply(X1,2,mean)
m2 <- apply(X2,2,mean)
m3 <- apply(X3,2,mean)
mall <- apply(Xall,2,mean)

# クラスの全標本が、クラス重心にあるものとして
# 作った標本データ
X1. <- matrix(rep(m1,each=n1),ncol=d)
X2. <- matrix(rep(m2,each=n2),ncol=d)
X3. <- matrix(rep(m3,each=n3),ncol=d)

Xall. <- rbind(X1.,X2.,X3.)
ss1 <- ss2 <- ss3 <- ssbetween <- ssall <- matrix(0,d,d)
for(i in 1:n1){
	ss1 <- ss1 + matrix(X1[i,]-m1,ncol=1) %*% matrix(X1[i,]-m1,nrow=1)
}
for(i in 1:n2){
	ss2 <- ss2 + matrix(X2[i,]-m2,ncol=1) %*% matrix(X2[i,]-m2,nrow=1)
}

for(i in 1:n3){
	ss3 <- ss3 + matrix(X3[i,]-m3,ncol=1) %*% matrix(X3[i,]-m3,nrow=1)
}
for(i in 1:(n1+n2+n3)){
	ssbetween <- ssbetween + matrix(Xall.[i,]-mall,ncol=1) %*% matrix(Xall.[i,]-mall,nrow=1)
	ssall <- ssall + matrix(Xall[i,]-mall,ncol=1) %*% matrix(Xall[i,]-mall,nrow=1)
}

# 実質的にゼロ行列
ssall - ssbetween - (ss1+ss2+ss3)

ssratio <- solve(ssall) %*% ssbetween

eigen.out <- eigen(ssratio)


#cv1 <- cov(X1)
#cv2 <- cov(X2)
#cv3 <- cov(X3)

#cvwithin <- cv1 + cv2 + cv3

cvbetween <- cov(rbind(X1.,X2.,X3.))

cvall <- cov(Xall)



tmp



M <- rbind(m1,m2,m3)
cvm <- cov(M)


cvall <- cov(Xall)

