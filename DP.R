library(DPpackage)
library(MCMCpack)

################ 背景分布
n.cl <- 5 # クラスタ数
rs <- rdirichlet(1,rep(5,n.cl)) # クラスタごとの割合
ps <- sort(runif(n.cl)*0.4+0.2) # 各クラスタの表確率

plot(ps,rs,type="h",xlim=c(0,1))
################
################推定条件
prior4<-list(a0=1,
						b0=1,
						a1=1,
						b1=1)
prior3<-list(alpha=100,
						a1=1,
						b1=1)
# Initial state
state <- NULL
# MCMC parameters
ngrid <- 100
nsave <- 1000
mcmc2 <- list(nburn=500,
						 nsave=nsave,
						 nskip=3,
						 ndisplay=100,
						 ngrid=ngrid)
################
################観察条件

N.samples <- c(10,100,1000)
N.trials <- c(10,100,10000)

NN.st <- expand.grid(N.samples,N.trials)
###############

###############推定結果記録
densp.ms <- matrix(0,length(NN.st[,1]),ngrid)
cl.states <- matrix(0,length(NN.st[,1]),nsave)
alphas <- matrix(0,length(NN.st[,1]),nsave)
###############
for(ii in 1:length(NN.st[,1])){
	N.sample <- NN.st[ii,1] # サンプル数
	N.trial <- NN.st[ii,2] # 各サンプルでの試行回数

	# クラスタを割り当てる
	cl <- sample(1:n.cl,N.sample,replace=TRUE,prob=rs)
	# 観察データを作る
	y <- matrix(0,N.sample,2)
	for(i in 1:N.sample){
		tmp <- sample(1:0,N.trial,replace=TRUE,prob=c(ps[cl[i]],1-ps[cl[i]]))
		y[i,1] <- sum(tmp)
		y[i,2] <- N.trial
	}


	# Fitting the model

	fit4 <- DPbetabinom(y=y,ngrid=100,
	                       prior=prior3,
	                       mcmc=mcmc2,
	                       state=state,
	                       status=TRUE)

	densp.ms[ii,] <- fit4$densp.m
	cl.states[ii,] <- fit4$save.state$thetasave[,1]
	alphas[ii,] <- fit4$save.state$thetasave[,2]

}

par(mfrow=c(3,3))
for(ii in 1:length(NN.st[,1])){
	ttl <- paste("Nsample",NN.st[ii,1]," Niter ",NN.st[ii,2])
	plot(fit4$grid,densp.ms[ii,],type="l",main=ttl)
}
