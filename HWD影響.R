# フェノタイプが大小二峰性
# ジェノタイプがHWE,HWDの２通り
# 帰無仮説
# fig1
n <- c(100,30)
P <- c(rnorm(n[1],0,1),rnorm(n[2],5,1))

hist(P)

hwd <- c(0.45,0.1,0.45)
hwe <- c(0.25,0.5,0.25)

n.snp <- 10^3
p.hwd <- p.hwe <- rep(0,n.snp)
for(i in 1:n.snp){
	G.hwd <- sample(0:2,sum(n),prob=hwd,replace=TRUE)
	G.hwe <- sample(0:2,sum(n),prob=hwe,replace=TRUE)
	
	tmp1 <- lm(P ~ G.hwd)
	tmp2 <- lm(P ~ G.hwe)
	p.hwd[i] <- anova(tmp1)$Pr[1]
	p.hwe[i] <- anova(tmp2)$Pr[1]
}

matplot(-log(cbind(sort(p.hwd),sort(p.hwe))),type="l")

# 対立仮説
# fig2
n <- c(100,30)
P <- c(rnorm(n[1],0,1),rnorm(n[2],5,1))

hist(P)

hwd <- c(0.45,0.1,0.45)
hwe <- c(0.25,0.5,0.25)

n.snp <- 10^3
p.hwd <- p.hwe <- rep(0,n.snp)
for(i in 1:n.snp){
	G.hwd <- sample(0:2,sum(n),prob=hwd,replace=TRUE)
	G.hwe <- sample(0:2,sum(n),prob=hwe,replace=TRUE)
	# 対立仮説。ジェノタイプが「線形に」フェノタイプ値に影響
	P.hwd <- rnorm(sum(n),G.hwd,5)
	P.hwe <- rnorm(sum(n),G.hwe,5)
	tmp1 <- lm(P.hwd ~ G.hwd)
	tmp2 <- lm(P.hwe ~ G.hwe)
	p.hwd[i] <- anova(tmp1)$Pr[1]
	p.hwe[i] <- anova(tmp2)$Pr[1]
}

matplot(-log(cbind(sort(p.hwd),sort(p.hwe))),type="l")

# 対立仮説
# fig3
n <- c(100,30)
P <- c(rnorm(n[1],0,1),rnorm(n[2],5,1))

hist(P)

hwd <- c(0.2,0.0,0.8)
hwe <- c(0.04,0.32,0.64)

n.snp <- 10^3
p.hwd <- p.hwe <- rep(0,n.snp)
for(i in 1:n.snp){
	G.hwd <- sample(0:2,sum(n),prob=hwd,replace=TRUE)
	G.hwe <- sample(0:2,sum(n),prob=hwe,replace=TRUE)
	# 対立仮説。ジェノタイプが「線形に」フェノタイプ値に影響
	P.hwd <- rnorm(sum(n),G.hwd,5)
	P.hwe <- rnorm(sum(n),G.hwe,5)
	tmp1 <- lm(P.hwd ~ G.hwd)
	tmp2 <- lm(P.hwe ~ G.hwe)
	p.hwd[i] <- anova(tmp1)$Pr[1]
	p.hwe[i] <- anova(tmp2)$Pr[1]
}

matplot(-log(cbind(sort(p.hwd),sort(p.hwe))),type="l")
