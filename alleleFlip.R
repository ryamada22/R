base <- c("A","T","G","C")

n.snp <- 10

X1 <- X2 <- matrix("",n.snp,2)
for(i in 1:n.snp){
	X1[i,] <- sample(base,2)
	X2[i,] <- sample(X1[i,])
}
X1
X2
# refrence base を基準にフリップの有無を確認する
same.snp <- which(X1[,1]==X2[,1])
flip.snp <- which(X1[,1]!=X2[,1])

X1. <- X1
X2. <- X2
X2.[flip.snp,] <- cbind(X2[flip.snp,2],X2[flip.snp,1])

X1.
X2.

n.sample = 5
G1 <- G2 <- matrix("",n.snp,n.sample)
genotype <- c("0/0","0/1","1/1")
G1 <- matrix(sample(genotype,n.snp*n.sample,replace=TRUE),ncol=n.sample)
G2 <- matrix(sample(genotype,n.snp*n.sample,replace=TRUE),ncol=n.sample)

# G2のジェノタイプを変換する
G2.flip <- G2
for(i in 1:length(flip.snp)){
	tmp <- G2[flip.snp[i],]
	tmp.flip <- rep("",n.sample)
	g00 <- which(tmp==genotype[1])
	g01 <- which(tmp==genotype[2])
	g11 <- which(tmp==genotype[3])
	
	tmp.flip[g00] <- genotype[3]
	tmp.flip[g01] <- genotype[2] # 変えず
	tmp.flip[g11] <- genotype[1]
	G2.flip[flip.snp[i],] <- tmp.flip
}

G1
G2
G2.flip

