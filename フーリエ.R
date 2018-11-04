t <- seq(from=0,to=1,length=10000)
x1 <- x2 <- rep(0,length(t))
x1[1:(length(t)/2+150)] <- 2
x1 <- x1/sum(x1)
x2[1:(length(t)/2+50)] <- 2
x2 <- x2/sum(x2)
x3 <- x1 + t*0.1
x3 <- x3/sum(x3)
fft1 <- fft(x1)
fft2 <- fft(x2)
fft3 <- fft(x3)

x11 <- x1 * rnorm(length(t),1,0.01)
x11 <- x11/sum(x11)
x22 <- x2 * rnorm(length(t),1,0.01)
x22 <- x22/sum(x22)
x33 <- x3 * rnorm(length(t),1,0.01)
x33 <- x33/sum(x33)
fft11 <- fft(x11)
fft22 <- fft(x22)
fft33 <- fft(x33)

y1 <- round((t^2),1)
y2 <- round((t^2.5),1)
y1 <- round(sin(t^2),10)
y2 <- round(sin(t^2.1),10)
ffty1 <- fft(y1)
ffty2 <- fft(y2)

par(mfcol=c(3,2))
matplot(cbind(y1,y2),type="l")
plot(y1,y2)
plot(Re(ffty1),Re(ffty2))
plot(Im(ffty1),Im(ffty2))


plot(ffty1,pch=20)
points(ffty2,pch=20,col=2)
sum((y1-y2)^2)*(t[2]-t[1])

sum(Mod(ffty1-ffty2)^2)

nt <- 100
t <- seq(from=0,to=3,length=nt)

y1 <- round(sin(t^2)*0.5,1)
y2 <- round(sin(t^2.1)*0.5,1)
ffty1 <- fft(y1)
ffty2 <- fft(y2)

par(mfcol=c(3,2))
matplot(cbind(y1,y2),type="l")
plot(y1,y2)
matplot(Re(cbind(ffty1,ffty2)),type="l")
matplot(Im(cbind(ffty1,ffty2)),type="l")

plot(ffty1,pch=20)
points(ffty2,pch=20,col=2)
sum((y1-y2)^2)*(t[2]-t[1])

sum(Mod(ffty1-ffty2)^2)


