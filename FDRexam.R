m <- 100
u <- 50
v <- m-u

n1 <- 100
n2 <- 100

my.means <- function(m,u,M,S){
	m1 <- rep(0,m)
	m2 <- c(rep(0,m-u),rnorm(u,M,S))
	return(cbind(m1,m2))
}

M <- 0.1
S <- 1

ms <- my.means(m,u,M,S)

my.X <- function(m,u,M,S,n1,n2){
	ms <- my.means(m,u,M,S)
	ret <- matrix(0,n1+n2,m)
	for(i in 1:m){
		tmp1 <- rnorm(n1,ms[i,1],1)
		tmp2 <- rnorm(n2,ms[i,2],1)
		ret[,i] <- c(tmp1,tmp2)
	}
	return(ret)
}

my.t.test.multi <- function(X,P){
	p <- rep(0,length(X[1,]))
	for(i in 1:length(p)){
		p[i] <- t.test(X[,i]~P)$p.value
	}
	q <- p.adjust(p,method="BH")
	return(list(p=p,q=q))
}
X <- my.X(m,u,M,S,n1,n2)
P <- c(rep(0,n1),rep(1,n2))

pq <- my.t.test.multi(X,P)

