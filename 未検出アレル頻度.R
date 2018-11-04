# KŽí—Þ
# ŠÏŽ@‚Ík < KŽí—Þ
# ŽŸ‚ªVŽí‚Å‚ ‚éŠm—¦‚ÍH

# K=2
# n=c(n1,...,nk)

library(MCMCpack)
my.prob.new.type <- function(n,N=100,N2=10000){
	p <- seq(from=0,to=1,length=N)
	ret <- rep(0,N)
	for(i in 1:N){
		q <- (1-p[i]) * rdirichlet(N2,rep(1,length(n)))
		ret[i] <- mean(n* apply(log(q),1,sum))
	}
	ret
}


out <- my.prob.new.type(n=c(1,1))
out2 <- my.prob.new.type(n=c(2,0))
out3 <- my.prob.new.type(n=c(4,3,2,1))
out4 <- my.prob.new.type(n=c(4,3,2,1,0))

plot(exp(out3))
points(exp(out2),col=2)


R1 <- rdirichlet(10^5,c(5,3,2,1)+1)
R2 <- rdirichlet(10^5,c(5,3,1,1,1)+1)
R3 <- rdirichlet(10^5,c(5,3,2,1,0)+1)
R4 <- rdirichlet(10^5,c(5,1,1,1,1,1,1,0,0,0,0,0,0,0)+1)


R3. <- R3[,5]
R4. <- apply(R4[,8:14],1,sum)

plot(sort(R3.),sort(R4.))



mean(R1[,1])
mean(R2[,1])
mean(R3[,1])
mean(R4[,1])

plot(sort(R1[,1]),sort(R4[,1]))