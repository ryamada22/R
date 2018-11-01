my.CGcoef <- function(L1,m,L2,j,k){
	if(k<0){
		return(my.CGcoef(L1,-m,L2,j,-k) * (-1)^(j-L1-L2))
	}
	if(L1 < L2){
		return(my.CGcoef(L2,m,L1,j,k) * (-1)^(j-L1-L2))
	}
	m1 <- m
	m2 <- k-m
	tmp1 <- sqrt((2*j+1)*factorial(j+L1-L2)*factorial(j-L1+L2)*factorial(L1+L2-j)/factorial(L1+L2+j+1))
	tmp2 <- sqrt(factorial(j+k)*factorial(j-k)*factorial(L1+m1)*factorial(L1-m1)*factorial(L2-m2)*factorial(L2+m2))
	min.k <- max(c(0,-j+L2-m1,-j+L1+m2))
	max.k <- min(c(L1+L2-j,L1-m1,L2+m2))
	tmp3 <- 0
	if(min.k <= max.k){
		for(i in min.k:max.k){
			tmp3 <- tmp3 + (-1)^i/(factorial(i)*factorial(L1+L2-j-i)*factorial(L1-m1-i)*factorial(L2+m2-i)*factorial(j-L2+m1+i)*factorial(j-L1-m2+i))
		}
	}
	ret <- tmp1*tmp2*tmp3
	return(ret)
}


my.CGcoef.wiki <- function(j1,j2,m1,m2,j){
	return(my.CGcoef(j1,m1,j2,j,m1+m2))
}

my.CGcoef.wiki(1/2,1/2,1/2,1/2,1)
my.CGcoef.wiki(1/2,1/2,-1/2,-1/2,1)
my.CGcoef.wiki(1/2,1/2,1/2,-1/2,1)
my.CGcoef.wiki(1/2,1/2,1/2,-1/2,0)
my.CGcoef.wiki(1/2,1/2,-1/2,1/2,1)
my.CGcoef.wiki(1/2,1/2,-1/2,1/2,0)
my.CGcoef.wiki(1,1,1,1,2)
my.CGcoef.wiki(1,1,1,0,2)
my.CGcoef.wiki(1,1,1,0,1)
my.CGcoef.wiki(1,1,0,1,2)
my.CGcoef.wiki(1,1,0,1,1)
my.CGcoef.wiki(1,1,1,-1,2)
my.CGcoef.wiki(1,1,1,-1,1)
my.CGcoef.wiki(1,1,1,-1,0)
my.CGcoef.wiki(1,1,0,0,2)
my.CGcoef.wiki(1,1,0,0,1)
my.CGcoef.wiki(1,1,0,0,0)
my.CGcoef.wiki(1,1,-1,1,2)
my.CGcoef.wiki(1,1,-1,1,1)
my.CGcoef.wiki(1,1,-1,1,0)


my.CGcoef.wiki(3/2,1,-1/2,1,5/2)

my.CGcoef.wiki(1,1,-1,1,1)
my.CGcoef.wiki(1,1,-1,1,0)

my.CGcoef(0,0,0,0,0)
my.CGcoef(1/2,1/2,1/2,1,1)
my.CGcoef(1/2,-1/2,1/2,1,-1)


my.comp.complex.mom.form <- function(L1,L2,k,j,C){
	ret <- 0
	for(i in 1:(2*L1+1)){
		m <- i-L1-1
		tmp <- my.CGcoef(L1,m,L2,j,k)
		ret <- ret + tmp * C[[L1+1]][i] * C[[L2+1]][k-m+L2+1]
	}
	return(ret)
}

C <- list()
C[[1]] <- 1
C[[2]] <- rep(1,3)
C[[3]] <- rep(1,5)
C[[4]] <- rep(1,7)
my.comp.complex.mom.form(2,2,0,0,C)

my.form2.product <- function(L1,L2,L3,L4,j){
	ret <- 0
	for(k in (-j):j){
		ret <- ret + (-1)^(j-k)*my.comp.complex.mom.form(L1,L2,k,j) * my.comp.complex.mom.form(L3,L4,-k,j)
	}
	return(ret/sqrt(2*j+1))
}

my.form.mom.product <- function(L1,L2,j,C){
	ret <- 0
	for(k in (-j):j){
		ret <- ret + (-1)^(j-k) * my.comp.complex.mom.form(L1,L2,k,j) * C[[j+1]][-k+j+1]
	}
	return(ret/sqrt(2*j+1))
}

vec2list <- function(vec){
  L <- ceiling(sqrt(length(vec))) - 1
  vec0 <- append(vec, rep(0, (L+1)^2 - length(vec)))
  vec1 <- tapply(vec0, unlist(mapply(rep, 0:L, each=2*(0:L)+1)), c)
  return(vec1)
}

data <- read.table("rotatedCoef.txt")
data <- as.matrix(data)
coef.lists <- list()
for(i in 1:100){
	coef.lists[[i]] <- vec2list(data[i,])
}

sapply(coef.lists,my.comp.complex.mom.form,L1=0,L2=0,k=0,j=0)
sapply(coef.lists,my.comp.complex.mom.form,L1=0,L2=1,k=0,j=0)
