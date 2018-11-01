# Dirichlet process
alpha <- 5

n <- 1000000
v <- rbeta(n,1,alpha)
hist(v)

ps <- rep(0,n)


logv <- log(v)
log1v <- log(1-v)
cumsumlog1ps <- cumsum(log1v)
ps[1] <- logv[1]
logv. <- logv[-1]
tmp <- logv. + cumsumlog1ps[-n]
ps <- c(logv[1], tmp)
#for(i in 2:n){
#	ps[i] <- logv[i] + sum(log1v[1:(i-1)])
#}
ps <- exp(ps)


hist(ps)
sum(ps)
plot(cumsum(ps),type="l")

plot(ps,type="l")

thetas <- rgamma(n,2,4)

par(mfcol=c(1,2))
plot(thetas,ps,type="h")
hist(thetas,col=2,freq=TRUE,density=20)
par(mfcol=c(1,1))



