n1 <- 1000
n2 <- 50
m1 <- 0
m2 <- 2
s1 <- 2
s2 <- 4

r1 <- rnorm(n1,m1,s1)
r2 <- rnorm(n2,m2,s2)

hist(c(r1,r2))
hist(r1,add=TRUE,col=2,density=17)
hist(r2,add=TRUE,col=3,density=19)
v <- c(rep(1,n1),rep(2,n2))

ord <- order(c(r1,r2))

plot(c(r1,r2)[ord]-min(c(r1,r2)),col=v[ord],type="h")

max.benefit <- sum(rank(c(r1,r2))*c(r1,r2))


n.iter <- 1000
bn <- rep(0,n.iter) # simple shuffle
bn2 <- bn # intra group shuffle
bn3 <- bn # scoring according to distributions
for(i in 1:n.iter){
	bn[i] <- sum(sample(1:(n1+n2))*c(r1,r2))
	bn2[i] <- sum(c(sample(1:n1),n1+sample(1:n2)) * c(r1,r2))
	tmp.val <- c(rnorm(n1,m1,s1),rnorm(n2,m2,s2))
	bn3[i] <- sum(rank(tmp.val)*c(r1,r2))
}

hist(bn)
hist(bn2,add=2,col=2)

hist(c(bn,bn2,bn3))
hist(bn,add=TRUE,col=2,density=17)
hist(bn2,add=TRUE,col=3,density=22)
hist(bn3,add=TRUE,col=4,density=25)





plot(c(r1,r2),v)

kari <- c(sample(r1),sample(r2))

plot(c(sample(1:n1),n1+sample(1:n2)),c(r1,r2))

cor(c(r1,r2),kari)

n.iter <- 1000
cors <- rep(0,n.iter)
for(i in 1:n.iter){
	kari2 <- c(rnorm(n1,m1,s1),rnorm(n2,m2,s2))

	#plot(c(r1,r2),kari2)
	cors[i] <- cor(c(r1,r2),kari2)
}

hist(cors)
abline(v=cor(c(r1,r2),kari),col=2)



