my.Ewens.prob <- function(ms,theta,log=FALSE){
	n <- sum(1:length(ms)*ms)
	if(length(ms) < n){
		ms <- c(ms,rep(0,n-length(ms)))
	}
	tmp <- lfactorial(n) -sum(log((0:(n-1))+theta)) + sum(ms * log(theta))-sum(ms*log(1:n))-sum(lfactorial(ms))
	if(log){
		return(tmp)
	}else{
		return(exp(tmp))
	}
}

my.Ewens.prob(c(0,1),1)
my.Ewens.prob(c(2,0),1)

library(partitions)
n <- 6
theta <- 0.1
prts <- parts(n)
prbs <- rep(0,length(prts[1,]))
for(i in 1:length(prts[1,])){
	prt <- prts[,i]
	print(prt)
	tab <- tabulate(prt)
	prbs[i] <- my.Ewens.prob(tab,theta)
}

sum(prbs)

my.Prob.numallele <- function(n,theta){
	prts <- parts(n)
	prbs <- rep(0,length(prts[1,]))
	for(i in 1:length(prts[1,])){
		prt <- prts[,i]
		print(prt)
		tab <- tabulate(prt)
		prbs[i] <- my.Ewens.prob(tab,theta)
	}
	ks <- apply(prts,2,function(x){length(which(x!=0))})
	K <- sum(prbs*ks)
	ks2 <- 1:n
	prbs2 <- rep(0,n)
	for(i in 1:n){
		prbs2[i] <- sum(prbs[which(ks==i)])
	}
	return(list(K=K,ks=ks,prbs=prbs,ks2=ks2,prbs2=prbs2,prts=prts))
}
library(gmp) # Stirling1()
my.Prob.numallele2 <- function(n,theta){
	ks <- 1:n
	prbs <- rep(0,n)
	for(i in 1:n){
		tmp <- ks[i] * log(theta) - sum(log((0:(n-1)) + theta))
		prbs[i] <- abs(as.integer(Stirling1(n,ks[i]))) * exp(tmp)
	}
	K <- sum(ks*prbs)
	return(list(K=K,ks=ks,prbs=prbs))
}

# singleton ++
n <- 3
theta <- 0.3

tmp.prbs <- my.Prob.numallele(n,theta)$prbs2

tmp.prbs2 <- my.Prob.numallele(n+1,theta)$prbs2

a <- tmp.prbs- tmp.prbs2[2:(n+1)]
b <- tmp.prbs - a

a+b - tmp.prbs
a/tmp.prbs
b/tmp.prbs

my.singleton.pp <- function(n,theta){
	tmp <- my.Prob.numallele2(n,theta)
	tmp2 <- my.Prob.numallele2(n+1,theta)
	prbs <- tmp$prbs
	prbs2 <- tmp2$prbs

	a <- b <- rep(0,n)
	for(i in n:1){
		if(i == n){
			a[i] <- prbs2[i+1]
			b[i] <- prbs[i]-a[i]
		}else{
			a[i] <- prbs2[i+1]-b[i+1]
			b[i] <- prbs[i]-a[i]
		}
	}
	return(list(a=a,b=b))
	
}
n <- 5
theta <- runif(1)
out <- my.singleton.pp(n,theta)
out[[1]]/(out[[1]]+out[[2]])
theta/(theta+n)

n <- 3
st <- rep(0,n)
for(i in 1:n){
	st[i] <- as.integer(Stirling2(n,i))
}



my.Ewens.dist <- function(bs,theta,log=FALSE){
	k <- length(bs)
	n <- sum(bs)
	tmp <- k * log(theta) -sum(log((0:(n-1))+theta)) + sum(lfactorial(bs-1))
	if(log){
		return(tmp)
	}else{
		return(exp(tmp))
	}
}


my.Ewens.dist(c(3),theta)
p1 <- my.Ewens.dist(c(3),theta)
p2 <- my.Ewens.dist(c(1,2),theta)
p3 <- my.Ewens.dist(c(1,1,1),theta)
p1 + p2*3 + p3

my.Ewens.dist2 <- function(bs,theta,log=FALSE){
	n <- sum(bs)
	
	tab <- tabulate(bs)
	tmp <- my.Ewens.prob(tab,theta,log=TRUE)
	
	tmp2 <- tmp -lfactorial(n) + sum(tab*lfactorial(1:length(tab))) + sum(lfactorial(tab))
	if(log){
		return(tmp2)
	}else{
		return(exp(tmp2))
	}

}
p1. <- my.Ewens.dist2(c(3),theta)
p2. <- my.Ewens.dist2(c(1,2),theta)
p3. <- my.Ewens.dist2(c(1,1,1),theta)
p1. + p2.*3 + p3.

my.list.parts <- function(n,theta){
	lst <- listParts(n)
	ret <- rep(0,length(lst))
	for(i in 1:length(ret)){
		ret[i] <- my.Ewens.dist2(sapply(lst[[i]],length),theta)
	}
	return(list(lst=lst,prob=ret))
}

n <- 3
theta <- 0.3
my.list.parts(n,theta)

# Ewens-Pitman sampling formula

