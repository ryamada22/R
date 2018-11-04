# P
p <- c(0, -50, 50)

x <- 2^(seq(from=0,to=15,length=100))
x <- c(-x,0,x)
y <- x^2/4
y2 <- x^2
z <- rep(0,length(x))

# canvas‚Í y=0
photo <- matrix(0,length(x),2)
photo2 <- photo
for(i in 1:length(x)){
	vx <- x[i]-p[1]
	vy <- y[i]-p[2]
	vy2 <- y2[i]-p[2]
	vz <- z[i]-p[3]
	t <- (-p[2])/vy
	t2 <- (-p[2])/vy2
	photo[i,1] <- p[1] + vx*t
	photo[i,2] <- p[3] + vz*t
	photo2[i,1] <- p[1] + vx*t2
	photo2[i,2] <- p[3] + vz*t2
}

plot(photo,type="l")
points(photo2,col=2,type="l")


