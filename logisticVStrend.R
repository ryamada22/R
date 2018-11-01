lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}



n.iter <- 100
N <- 100
lg.p <- trend.p <- rep(0,n.iter)

for(i in 1:n.iter){
	ph <- sample(0:1,N,replace=TRUE,prob=c(0.5,0.5))
	g <- sample(0:1,N,replace=TRUE)
	
	ph <- sort(ph)
	g <- sort(g)
	g[1:(7/10*N)] <- sample(g[1:(7/10*N)])
	
	logistic.out <- glm(ph~g,family=binomial())
	lg.p[i] <- summary(logistic.out)$coefficients[2,4]
	trend.p[i] <- lmp(lm(ph~g))
}

plot(lg.p,trend.p)
abline(0,1,col=2)

plot(log(lg.p),log(trend.p))
abline(0,1,col=2)

