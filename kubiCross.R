n1 <- 4
n2 <- 3
n3 <- 5
M1 <- matrix((1:(n1*n2)),n1,n2)
M2 <- matrix((1:(n2*n3)),n3,n2)
my.sh.cross.kubi <- function(M1,M2){
	n1 <- length(M1[,1])
	n2 <- length(M2[,1])
	ret <- matrix(0,n1,n2)
	for(i in 1:n1){
		for(j in 1:n2){
			ret[i,j] <- min(M1[i,]+M2[j,])
		}
	}
	return(ret)
}

my.sh.cross.kubi(M1,M2)

