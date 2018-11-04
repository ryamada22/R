n <- 100
x <- 1:n
v <- sample(0:1,n,replace=TRUE)

v <- rep(0,n)
k <- 10
L <- 4
ss <- sample(x,k)
for(i in 1:k){
	v[ss[i]:(ss[i]+L)] <- 1
}


dx <- as.matrix(dist(x))
dv <- as.matrix(dist(v))

image(dx)
image(dv)

dxval <- sort(unique(c(dx)))

prdxval <- (dxval/max(dxval))^2

plot(dxval,prdxval)

n.iter <- 10

hxg <- list()
hxp <- list()
hxg[[1]] <- dv
hxp[[1]] <- prdxval
for(i in 1:n.iter){
	hxg[[i+1]] <- matrix(0,n,n)
	for(j in 1:n){
		for(k in 1:n){
			tmpdx <- dx[j,k]
			r <- runif(1)
			if(hxg[[i]][j,k]==1){
				if(r > hxp[[i]][which(dxval==tmpdx)]){
					hxg[[i+1]][j,k] <- 1
				}
			}else{
				if(r <= hxp[[i]][which(dxval==tmpdx)]){
					hxg[[i+1]][j,k] <- 1
				}
			}

		}
	}
	hxp[[i+1]] <- prdxval
	for(j in 1:length(prdxval)){
		tmpadd <- hxg[[i]][which(dx==dxval[j])]
		hxp[[i+1]][j] <- sum(tmpadd)/length(tmpadd)
	}
}

matplot(matrix(unlist(hxp),ncol=n.iter+1),type="l")

