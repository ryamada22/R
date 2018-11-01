library(MCMCpack)
R <- rdirichlet(10^3,rep(1,3))

n <- 4
R1 <- rdirichlet(10^5,c(rep(10,n),1,0)+1)
R2 <- rdirichlet(10^5,c(10*n-(n-1),rep(1,n-1),1,rep(0,100))+1)

plot(sort(R1[,n+1]),sort(R2[,n+1]))

abline(0,1,col=2)

R3 <- rdirichlet(10^5,c(3,3,3,3,1)+1)
R4 <- rdirichlet(10^5,c(3,3,3,3,1,rep(0,100))+1)

m3 <- apply(R3,2,mean)
m4 <- apply(R4,2,mean)


R5 <- rdirichlet(10^5,c(100,100,0)+1)
R6 <- rdirichlet(10^5,c(100,rep(1,100),rep(0,100))+1)

plot(sort(R5[,3]),sort(R6[,102]))

