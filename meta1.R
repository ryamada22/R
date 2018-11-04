library(rgl)
# sample ”
n <- 10^4
# locus ”
n.locus <- 10^2
trues <- matrix(0,n.locus,8)
# regression results
lmoutXY <- lmoutXZ <- lmoutYZ <-list()
# n.locus‚É‹¤’Ê‚ÌX‚ğg‚¤
mx <- rnorm(1,0,1)
sx <- abs(rnorm(1,0,1))
X <- rnorm(n,mx,sx)
for(i in 1:n.locus){

ay <- rnorm(1,0,1)
by <- rnorm(1,0,1)
sy <- abs(rnorm(1,0,1))
Y <- ay * X + by + rnorm(n,0,sy)

az <- rnorm(1,0,abs(ay))
bz <- rnorm(1,0,1)
sz <- abs(rnorm(1,0,1))
Z <- az * X + by + rnorm(n,0,sz)

trues[i,] <- c(mx,sx,ay,by,sy,az,bz,sz)

lmoutXY[[i]] <- lm(Y~X)
lmoutXZ[[i]] <- lm(Z~X)
lmoutYZ[[i]] <- lm(Z~Y)
#plot3d(X,Y,Z)
}
