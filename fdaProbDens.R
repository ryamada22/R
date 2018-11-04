rangeval <- c(-3,3)
#  set up some standard normal data
x <- rnorm(50)
#  make sure values within the range
x[x < -3] <- -2.99
x[x >  3] <-  2.99
#  set up basis for W(x)
basisobj <- create.bspline.basis(rangeval, 11)
#  set up initial value for Wfdobj
Wfd0 <- fd(matrix(0,11,1), basisobj)
WfdParobj <- fdPar(Wfd0)
#  estimate density
denslisty <- density.fd(x, WfdParobj)

y <- rnorm(50)
#  make sure values within the range
y[y < -3] <- -2.99
y[y >  3] <-  2.99
#  set up basis for W(x)
basisobj <- create.bspline.basis(rangeval, 11)
#  set up initial value for Wfdobj
Wfd0 <- fd(matrix(0,11,1), basisobj)
WfdParobj <- fdPar(Wfd0)
#  estimate density
denslisty <- density.fd(y, WfdParobj)

inprod(denslistx[[1]],denslisty[[1]])

# https://rstudio-pubs-static.s3.amazonaws.com/170685_57aa963230cb4adb8d2684fdfdecca0c.html#/4 

library(fda.usc)

n.sample <- 10
n.pts <- 200

m <- matrix(runif(n.sample*n.pts)^2,ncol=n.pts)
m <- t(apply(m,1,cumsum))
m <- m/m[,n.pts]

N <- 200
m <- cbind(matrix(0,n.sample,N),m,matrix(1,n.sample,N))
matplot(t(m))

fm <- fdata(m)
par(mfcol=c(2,2))

plot(fm)

out1<-min.np(fm,type.CV=CV.S,type.S=S.NW)

sm.fm <- out1$fdata.est

plot(sm.fm)

m.d1 <- fdata.deriv(fm, nderiv = 1)

sm.m.d1 <- fdata.deriv(sm.fm, nderiv = 1)
plot(m.d1)

plot(sm.m.d1)

