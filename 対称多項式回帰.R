# xs : values in support
# as : coeffs of polynomial formula, with as[1] = a_0
my.poly <- function(xs,as){
	ret <- rep(0,length(xs))
	for(i in 1:length(as)){
		ret <- ret + as[i] * xs^(i-1)
	}
	return(ret)
}

my.symmetric.poly <- function(Xs,as){
	if(!is.matrix(Xs)){
		Xs <- matrix(Xs,ncol=1)
	}
	ret <- rep(01,length(Xs[,1]))
	for(i in 1:length(Xs[1,])){
		ret <- ret * my.poly(Xs[,i],as)
	}
	return(ret)
}

# partitions
library(partitions)
n <- 4 # number of variables
k <- 3 # the max degree of monovariable polynomial
p <- n * k # the max degree of multivariable polynomial

restrictedparts(p,n)

perms <- permutations(n,n)

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

my.polybeta.from.multpoly <- function(as,rp){
	p <- max(rp)
	n <- length(rp[,1])
	ret <- rep(0,p+1)
	for(i in 0:p){
		tmp <- which(apply(abs(rp-i),2,sum)==0)
		ret[i+1] <- abs(as[tmp])^(1/n)*sign(as[tmp])
	}
	return(ret)
}

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
	
	return(list(lm.out=lm.out,pred.as=pred.as,rp=rp,dupl=dupl))

}


n <- 4
n.pt <- 100
Xs <- matrix(rnorm(n.pt*n),ncol=n)

k <- 5
as <- rnorm(k+1)

y <- my.symmetric.poly(Xs,as)

out <- my.sym.poly.lm(Xs,y,k)

print(as)
print(out$pred.as)

plot(y,my.symmetric.poly(Xs,out$pred.as))

plot(y,predict(out$lm.out))


my.sym.poly.terms <- function(Xs,k){
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
	return(t(t(ret)/dupl))
}

Xs <- matrix(c(1,2,3,4,5,6),ncol=3)

n <- 3
n.pt <- 100
Xs <- matrix(rnorm(n.pt*n),ncol=n)

k <- 2
as <- rnorm(k+1)

y <- my.symmestric.poly(Xs,as)
#y <- y + apply(sin(Xs),1,sum) * sd(y)*10
#y <- apply(sin(Xs),1,sum)
Xs.poly.terms <- my.sym.poly.terms(Xs,k)

rp.dupl <- my.serial.partitions(n,k)

lm.out <- lm(y~Xs.poly.terms-1)
pre.as <- my.polybeta.from.multpoly(lm.out[[1]],rp.dupl[[1]])

preds <- predict(lm.out)

plot(y,preds)
abline(0,1,col=2)
