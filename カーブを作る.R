n.step <- 100

x <- rnorm(n.step,10,1) * sample(c(-1,1),n.step,replace=TRUE)
x. <- cumsum(x)

plot(x.)
max.x <- max(abs(x))

x.st <- x/max.x * 0.95
x.st. <- cumsum(x.st)

plot(x.st.,type="l")

D <- matrix(0,n.step,n.step)

for(i in 1:(n.step-1)){
	for(j in i:n.step){
		D[i,j] <- D[j,i] <- (i-j)^2 - (x.st.[i]-x.st.[j])^2
	}
}

image(D)

IP <- matrix(0,n.step,n.step)
M <- matrix(c(1,0,0,-1),2,2)
for(i in 1:n.step){
	for(j in 1:n.step){
		IP[i,j] <- matrix(c(i,x.st.[i]),ncol=2) %*% M %*% matrix(c(j,x.st.[j]),ncol=1)
	}
}


e.out <- eigen(IP)

my.curve.update <- function(x){
	n <- length(x[,1])
	
	s <- sample(1:(n-1),1)
	v1 <- x[s,]
	v2 <- x[s+1,]
	
	v3 <- c((v1[1]+v2[1])/2+(v1[2]-v2[2])/2,(v1[1]-v2[1])/2 + (v1[2]+v2[2])/2)
	v4 <- c((v1[1]+v2[1])/2-(v1[2]-v2[2])/2,-(v1[1]-v2[1])/2 + (v1[2]+v2[2])/2)


	r <- runif(1)

	newv <- r*v3 + (1-r)*v4
	ret <- rbind(x[1:s,],newv,x[(s+1):n,])
	print(ret)
	return(ret)
}
my.curve <- function(n){
	s <- runif(1)  - 0.5
	x <- rbind(c(0,0),c(1,s))
	
	for(i in 1:(n-2)){
		x <- my.curve.update(x)
	}
	return(x)
	
}

n.step <- 100
cv <- my.curve(n.step)
plot(cv,type="b",ylim=c(-0.5,0.5))




D <- matrix(0,n.step,n.step)

for(i in 1:(n.step-1)){
	for(j in i:n.step){
		D[i,j] <- D[j,i] <- sqrt((cv[i,1]-cv[j,1])^2 - (cv[i,2]-cv[j,2])^2)
	}
}

image(D)

IP <- matrix(0,n.step,n.step)
M <- matrix(c(1,0,0,-1),2,2)
for(i in 1:n.step){
	for(j in 1:n.step){
		IP[i,j] <- matrix(cv[i,],ncol=2) %*% M %*% matrix(cv[j,],ncol=1)
	}
}

image(IP)

e.out <- eigen(IP)
plot(e.out[[1]])

par(mfcol=c(1,2))
plot(cv)
plot(e.out[[2]][,c(1,n.step)])

