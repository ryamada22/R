my.ip <- function(v1,v2,M){
	matrix(v1,nrow=1) %*% M %*% matrix(v2,ncol=1)
}
my.ip.mat <- function(x,M){
	n <- length(x[,1])
	ret <- matrix(0,n,n)
	for(i in 1:n){
		for(j in 1:n){
			ret[i,j] <- my.ip(x[i,],x[j,],M)
		}
	}
	return(ret)
}
x <- matrix(c(0,0,2,1,5,3,8,4,4,-4),byrow=TRUE,ncol=2)

M <- diag(c(1,-1))

IP <- my.ip.mat(x,M)


IP

eout <- eigen(IP)

eout

newx <- eout[[2]][,1] * sqrt(eout[[1]][1])
newy <- eout[[2]][,5] * sqrt(-eout[[1]][5])

newX <- cbind(newx,newy)
newIP <- my.ip.mat(newX,M)
par(mfcol=c(1,2))
plot(x)

plot(newx,newy)


par(mfcol=c(1,1))

plot(rbind(x,newX),col=c(rep(1,5),rep(2,5)))
