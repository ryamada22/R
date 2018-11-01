n <- 10

k1 <- 3
k2 <- 5
k12 <- 2

V1 <- matrix(rnorm(n*k1),k1,n)
V2 <- rbind(V1[1:k12,],matrix(rnorm(n*(k2-k12)),k2-k12,n))

qr(V1)$rank
qr(V2)$rank
qr(rbind(V1,V2))$rank

n.iter <- 50

for(i in 1:n.iter){
	k1 <- sample(2:10,1)
	k2 <- sample(k2:10,1)
	k12 <- sample(0:min(k1,k2),1)
	n <- max(k1,k2)*2
	V1 <- matrix(rnorm(n*k1),k1,n)
	if(k12==0){
		V2 <- matrix(rnorm(n*k2),k2,n)
	}else{
		V2 <- rbind(V1[1:k12,],matrix(rnorm(n*(k2-k12)),k2-k12,n))
	}
	
	
	out <- my.unionrank(V1,V2)
	
	print(out$ru - k12)
	#print(out$ru)
	#print(k12)
	#print("---")
	#r1 <- qr(V1)$rank
	#r2 <- qr(V2)$rank
	#r12 <- qr(rbind(V1,V2))$rank
	#print("----")
	#print(n)
	#print(k1)
	#print(k2)
	#print(k12)
	#print(k1-r1)
	#print(k2-r2)
	#print(k1+k2-r12-k12)
}

my.unionrank <- function(V1,V2){
	r1 <- qr(V1)$rank
	r2 <- qr(V2)$rank
	r12 <- qr(rbind(V1,V2))$rank
	unionr <- r1 + r2 - r12
	
	return(list(r1=r1,r2=r2,r12=r12,ru=unionr))
}


