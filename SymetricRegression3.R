library(partitions)
library(gtools)
library(TruncatedNormal)
# calculate monovariate polynomial function value
# xs : values in support
# as : coeffs of polynomial formula, with as[1] = a_0
my.poly <- function(xs,as){
	ret <- rep(0,length(xs))
	for(i in 1:length(as)){
		ret <- ret + as[i] * xs^(i-1)
	}
	return(ret)
}
# calculate symmetric multivariate polynomial function value
# symmetric multivariate polynomial function is multiplication of
# monovariate polynomial functions with the same coefficients
my.symmetric.poly <- function(Xs,as){
	if(!is.matrix(Xs)){
		Xs <- matrix(Xs,ncol=1)
	}
	ret <- rep(1,length(Xs[,1]))
	for(i in 1:length(Xs[1,])){
		ret <- ret * my.poly(Xs[,i],as)
	}
	return(ret)
}

# n : number of variables
# k : the degree of monovariate polynomial
# returns the partitions of 0,1,...,n*k whose element values <= k
# and returns integer list that stands for duplication due to the exchangeability of variables
# whose index in terms are identical
my.serial.partitions <- function(n,k){
	p <- n*k
	rp <- matrix(0,n,1)
	for(i in 1:p){
		tmp <- restrictedparts(i,n)
		rp <- cbind(rp,tmp[,which(apply(tmp,2,max)<=k)])
	}
	dupl <- sapply(apply(rp,2,table),function(x)prod(factorial(x)))
	return(list(rp=rp,dupl=dupl))
}

# Calculates coefficients of monovariate polynomial from
# coefficients of symmetric multivariate polynomial
# as : coefficients of all symmetric terms that should be estimated by linear regression
# rp : patterns of interger partitions that are used in the multivariate symmetric polynomial
my.polybeta.from.multpoly <- function(as,rp){
	p <- max(rp)
	n <- length(rp[,1])
	ret <- rep(0,p+1)
	for(i in 0:p){
		tmp <- which(apply(abs(rp-i),2,sum)==0)
		if(i>0){
		  onlyi <- which(apply(rp==0,2,sum)==(n-1) & apply(rp,2,sum)==i)
		  ret[i+1] <- abs(as[tmp])^(1/n)*sign(as[onlyi])
		}else{
		  ret[i+1] <- abs(as[tmp])^(1/n)* sign(as[1])
		}
		
		#ret[i+1] <- abs(as[tmp])^(1/n)*sign(as[tmp])
		
	}
	return(ret)
}
# standardize
my.standardize <- function(x){
	m <- mean(x)
	s <- sum((x-m)^2)/length(x)
	(x-m)/sqrt(s)
}
# Linear regression of symmetric polynomial
# args : Xs, y, and k (max degree of monovariate polynomial)
# values :
# lm.out : output of lm() with all symmetric terms
# pred.as : estimated coefficients of MONOvariate poly
# rp : All partitions considered
# dupl: coefficents for all partitions to take care exchangeability of variables

my.sym.poly.lm <- function(Xs,y,k){
	n <- length(Xs[1,])
	rp.dupl <- my.serial.partitions(n,k)
	rp <- rp.dupl[[1]]
	dupl <- rp.dupl[[2]]

	pms <- permutations(n,n)
	ret <- matrix(0,length(Xs[,1]),length(rp[1,]))
	for(i in 1:length(rp[1,])){
		for(j in 1:length(pms[,1])){
			tmp <- rep(1,length(Xs[,1]))
			for(j2 in 1:n){
				tmp <- tmp * Xs[,pms[j,j2]]^rp[j2,i]
			}
			ret[,i] <- ret[,i] + tmp
		}
	}
	Xs.poly.terms <- t(t(ret)/dupl)
	lm.out <- lm(y~Xs.poly.terms-1)
	pred.as <- my.polybeta.from.multpoly(lm.out[[1]],rp)
	
	return(list(lm.out=lm.out,pred.as=pred.as,rp=rp,dupl=dupl,Xs.poly.terms=Xs.poly.terms))

}
my.sym.poly.lm.st <- function(Xs,y,k){
	n <- length(Xs[1,])
	rp.dupl <- my.serial.partitions(n,k)
	rp <- rp.dupl[[1]]
	dupl <- rp.dupl[[2]]

	pms <- permutations(n,n)
	ret <- matrix(0,length(Xs[,1]),length(rp[1,]))
	for(i in 1:length(rp[1,])){
		for(j in 1:length(pms[,1])){
			tmp <- rep(1,length(Xs[,1]))
			for(j2 in 1:n){
				tmp <- tmp * Xs[,pms[j,j2]]^rp[j2,i]
			}
			ret[,i] <- ret[,i] + tmp
		}
	}
	Xs.poly.terms <- t(t(ret)/dupl)
	Xs.poly.terms.st <- apply(Xs.poly.terms,2,my.standardize)
	Xs.poly.terms.st[,1] <- Xs.poly.terms[,1]
	lm.out <- lm(y~Xs.poly.terms.st-1)
	pred.as <- my.polybeta.from.multpoly(lm.out[[1]],rp)
	
	return(list(lm.out=lm.out,pred.as=pred.as,rp=rp,dupl=dupl,Xs.poly.terms=Xs.poly.terms,Xs.poly.terms.st=Xs.poly.terms.st))

}

d <- 2

A <- matrix(rnorm(d^2),d,d)
A <- A/sqrt(apply(A^2,1,sum))

Sigma <- A %*% t(A)
l <- rep(-1,d)
u <- rep(1,d)
mvNcdf(l,u,Sigma,1000)

d <- 6
l <- rep(-1,d)
u <- rep(1,d)

n.a <- 1000
ss <- matrix(0,n.a,d*(d-1)/2)
y <- rep(0,n.a)
volp <- rep(0,n.a)
for(i in 1:n.a){
	A <- matrix(rnorm(d^2),d,d)
	A <- A/sqrt(apply(A^2,1,sum))
	
	Sigma <- A %*% t(A)
	ss[i,] <- Sigma[upper.tri(Sigma)]
	y[i] <- mvNcdf(l,u,Sigma,1000)[[1]]
	A. <- cbind(rep(0,d),A)
	volp[i] <- abs(det(rbind(A.,rep(1,d+1))))
}

plot(volp,y)


library(rgl)
y.st <- (y-min(y))/(max(y)-min(y))
plot3d(ss)

spheres3d(ss,col=rgb(y.st,1-y.st,1),radius=0.01)
k <- 2
logy <- exp(y)
out <-my.sym.poly.lm.st(ss,y,k)
plot(y,predict(out$lm.out))

out.lars <- lars(y,out$Xs.poly.terms.st)

plot(out.lars)


plot(out[[1]][[1]])


out$rp[,which(abs(out[[1]][[1]]) > mean(abs(out[[1]][[1]])))]

ord <- order(abs(out[[1]][[1]]))
out$rp[,ord]

plot(y,out$Xs.poly.terms.st[,ord[84]])

##########

q <- 4
plot(y,apply(ss^q,1,prod))

library(lars)

out.lars <- lars(out$Xs.poly.terms.st,y)
