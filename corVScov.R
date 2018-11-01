library(mvtnorm)

sigma <- matrix(c(4,2,2,3), ncol=2)
n <- 10^4
x <- rmvnorm(n=n, mean=c(1,2), sigma=sigma)
colMeans(x)
var(x)
cor(x)
x. <- x
x.[,1] <- x.[,1] + rnorm(n,0.2)

var(x.)
cor(x.)

x.. <- x.
x.. <- x.. + rnorm(length(x..),0,10)

var(x..)
cor(x..)

plot(x..,asp=TRUE)
points(x.,col=2)
points(x,col=3)

