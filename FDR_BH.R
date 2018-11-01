p <- 10000
k <- 200

z <- rnorm(p)
for(i in 1:k){
	z[i] <- rnorm(1,rnorm(5,2),1)*sample(c(-1,1),1)
}

t <- seq(from=0,to=5,length=1000)

q <- 0.05

x <- rep(0,length(t))
y <- rep(0,length(t)) # no. selected features
w <- rep(0,length(t)) # no. true positive
for(i in 1:length(x)){
	x[i] <- p * pnorm(t[i],lower.tail=FALSE)*2/max(sum(abs(z)>=t[i]),1)
	
	y[i] <- sum(abs(z)>=t[i])
	tmp <- which(abs(z)>=t[i])
	w[i] <- y[i]-sum(tmp %in% 1:k)
	
}

plot(t,x)
abline(h=q,col=2)
y. <- y
y.[which(y==0)] <- 1
points(t,w/y.,col=3)

