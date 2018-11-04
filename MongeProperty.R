AB.0 <- matrix(sample(1:1000,6),ncol=2)
ord <- rbind(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
n <- length(ord[,1])

for(k in 1:n){
	for(kk in 1:n){
		AB <- cbind(AB.0[ord[k,],1],AB.0[ord[kk,],2])
		abcd <- matrix(0,0,4)
		val <- c()
		for(i in 1:2){
			for(ii in (i+1):3){
				for(j in 1:2){
					for(jj in (j+1):3){
						tmp <- d.[c(AB[i,1],AB[ii,1]),c(AB[j,2],AB[jj,2])]
						abcd <- rbind(abcd,c(i,ii,j,jj))
						#abcd <- rbind(abcd,c(AB[i,1],AB[ii,1],AB[j,2],AB[jj,2]))
						val <- c(val,tmp[1,1] + tmp[2,2] - (tmp[1,2]+tmp[2,1]))
					}
				}
			}
		}
		print(prod(sign(val)==1))
		if(prod(sign(val)==1)){
			print(val)
		}
	}
}


#cbind(abcd,val)

library(gtools)
n <- 4
vs <- sample(1:1000,n)
pm <- permutations(n,n)

pm2 <- permutations(n,2)
pm2 <- t(apply(pm2,1,sort))


for(i in 1:length(pm[,1])){
	tmp <- d.[vs[pm[i,]],vs[pm[i,]]]
	ret <- c()
	for(j in 1:length(pm2[,1])){
		for(jj in 1:length(pm2[,1])){
			tmp2 <- tmp[pm2[j,],pm2[jj,]]
			tmp3 <- tmp2[1,1] + tmp2[2,2] - (tmp2[1,2]+tmp2[2,1])
			ret <- c(ret,tmp3)
		}
	}
	print(prod(sign(ret)==1))
		if(prod(sign(ret)==1)){
			print(ret)
		}
}

