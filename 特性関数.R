my.charFx <- function(x,t){
	y <- rep(0,length(t))
	for(i in 1:length(y)){
		y[i] <- sum(exp(1i * t[i] * x))/length(x)
	}
	return(y)
}

mu <- 3
s <- 2
x <- rnorm(1000000,mu,s)

library(rgl)

t <- seq(from=-3,to=3,length=101)
y <- my.charFx(x,t)

plot3d(t,Re(y),Im(y))

y2 <- exp(1i * t * mu - s^2*t^2/2)

plot(Re(y),Re(y2))

y3 <- my.charFx(x^2,t)
