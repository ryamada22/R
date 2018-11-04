# tab 3 genotype count

my.hwe <- function(tab){
	p.exp <- (2*tab[1]+tab[2])/(2*sum(tab))
	g.exp <- c(p.exp^2,2*p.exp*(1-p.exp),(1-p.exp)^2)
	N <- sum(tab)
	chisq <- sum((tab-N*g.exp)^2/(N*g.exp))
	return(chisq)
}

n.iter <- 10^4

chisq.stock <- rep(0,n.iter)

for(i in 1:n.iter){
	p <-runif(1) * 0.6 + 0.2
	N <- sample(500:1000,1)
	obs <- sample(0:2,N,replace=TRUE,prob=g)
	tab <- table(obs)
	chisq <- my.hwe(tab)
	chisq.stock[i] <- chisq
}

hist(chisq.stock)

k <- 1
chisq.theory <- rchisq(n.iter,df=k)

plot(sort(chisq.stock),sort(chisq.theory))
