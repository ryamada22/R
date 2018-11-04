t <- seq(from=0,to=10,length=1000)

N0 <- 1
Ninf <- 10
gamma <- 2

Nt <- Ninf/(1+((Ninf/N0)-1) * exp(-gamma*t))

plot(t,Nt,type="l",ylim = c(0,Ninf+1))
abline(h=Ninf,col=2)

