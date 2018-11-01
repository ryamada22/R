# Data
library(DPpackage)
data(rolling)
y <- cbind(rolling$y1,rolling$y2)
# Prior information
prior<-list(alpha=1,
            a1=1,
            b1=1)
prior3<-list(alpha=10^(-4),
            a1=1,
            b1=1)
prior4<-list(a0=1,
            b0=1,
            a1=1,
            b1=1)
# Initial state
state <- NULL
# MCMC parameters
mcmc2 <- list(nburn=500,
             nsave=1000,
             nskip=3,
             ndisplay=100)
# Fitting the model
fit3 <- DPbetabinom(y=y,ngrid=100, 
                       prior=prior3, 
                       mcmc=mcmc2, 
                       state=state, 
                       status=TRUE)
fit4 <- DPbetabinom(y=y,ngrid=100, 
                       prior=prior4, 
                       mcmc=mcmc2, 
                       state=state, 
                       status=TRUE)

fit
summary(fit)
# density estimate
plot(fit,output="density")
# parameters
plot(fit,output="param")

# pの分布
plot(fit4$grid,fit4$densp.m,type="l")

# No. clusters の推移と分布
plot(fit4$save.state$thetasave[,1],type="l")
hist(fit4$save.state$thetasave[,1])

# alpha の推移と分布

plot(fit4$save.state$thetasave[,2],type="l")
plot(density(fit4$save.state$thetasave[,2]))

# 尤度の推移
plot(fit4$save.state$thetasave[,3],type="l")


# 条件を厳しくして回してみる
prior5<-list(alpha=10^(-4),
            a1=1,
            b1=1)


# Fitting the model
fit5 <- DPbetabinom(y=y,ngrid=100, 
                       prior=prior5, 
                       mcmc=mcmc2, 
                       state=state, 
                       status=TRUE)

# pの分布
plot(fit5$grid,fit5$densp.m,type="l")

# No. clusters の推移と分布
plot(fit5$save.state$thetasave[,1],type="l")
hist(fit5$save.state$thetasave[,1])


# 条件を厳しくして回してみる
prior6<-list(alpha=10^(-4),
            a1=100000,
            b1=10000)


# Fitting the model
fit6 <- DPbetabinom(y=y,ngrid=100, 
                       prior=prior6, 
                       mcmc=mcmc2, 
                       state=state, 
                       status=TRUE)

# pの分布
plot(fit6$grid,fit6$densp.m,type="l")

# No. clusters の推移と分布
plot(fit6$save.state$thetasave[,1],type="l")
hist(fit6$save.state$thetasave[,1])

# クラスタ数を巨大に
prior9<-list(a0=1,
            b0=1,
            a1=1,
            b1=1)


# Fitting the model
fit9 <- DPbetabinom(y=y,ngrid=100, 
                       prior=prior9, 
                       mcmc=mcmc2, 
                       state=state, 
                       status=TRUE)

# pの分布
plot(fit9$grid,fit9$densp.m,type="l")

# No. clusters の推移と分布
plot(fit9$save.state$thetasave[,1],type="l")
hist(fit9$save.state$thetasave[,1])
hist(fit9$save.state$thetasave[,2])

n.cl <- 5
library(MCMCpack)
rs <- rdirichlet(1,rep(5,n.cl))
ps <- runif(n.cl)*0.2+0.4

N.sample <- 100
N.trial <- 1000
cl <- sample(1:n.cl,N.sample,replace=TRUE,prob=rs)

y <- matrix(0,N.sample,2)
for(i in 1:N.sample){
	tmp <- sample(1:0,N.trial,replace=TRUE,prob=c(ps[cl[i]],1-ps[cl[i]]))
	y[i,1] <- sum(tmp)
	y[i,2] <- N.trial
}

hist(y[,1])
# Prior information
prior<-list(alpha=1,
            a1=1,
            b1=1)
prior3<-list(alpha=10^(-4),
            a1=1,
            b1=1)
prior4<-list(a0=1,
            b0=1,
            a1=1,
            b1=1)
# Initial state
state <- NULL
# MCMC parameters
mcmc2 <- list(nburn=500,
             nsave=1000,
             nskip=3,
             ndisplay=100)
# Fitting the model
fit3 <- DPbetabinom(y=y,ngrid=100, 
                       prior=prior3, 
                       mcmc=mcmc2, 
                       state=state, 
                       status=TRUE)
fit4 <- DPbetabinom(y=y,ngrid=100, 
                       prior=prior4, 
                       mcmc=mcmc2, 
                       state=state, 
                       status=TRUE)

fit
summary(fit)
# density estimate
plot(fit,output="density")
# parameters
plot(fit,output="param")

# pの分布
plot(fit4$grid,fit4$densp.m,type="l")

# No. clusters の推移と分布
plot(fit4$save.state$thetasave[,1],type="l")
hist(fit4$save.state$thetasave[,1])

# alpha の推移と分布

plot(fit4$save.state$thetasave[,2],type="l")
plot(density(fit4$save.state$thetasave[,2]))

# 尤度の推移
plot(fit4$save.state$thetasave[,3],type="l")


# 条件を厳しくして回してみる
prior5<-list(alpha=10^(-4),
            a1=1,
            b1=1)


# Fitting the model
fit5 <- DPbetabinom(y=y,ngrid=100, 
                       prior=prior5, 
                       mcmc=mcmc2, 
                       state=state, 
                       status=TRUE)

# pの分布
plot(fit5$grid,fit5$densp.m,type="l")

# No. clusters の推移と分布
plot(fit5$save.state$thetasave[,1],type="l")
hist(fit5$save.state$thetasave[,1])


# 条件を厳しくして回してみる
prior6<-list(alpha=10^(-4),
            a1=100000,
            b1=10000)


# Fitting the model
fit6 <- DPbetabinom(y=y,ngrid=100, 
                       prior=prior6, 
                       mcmc=mcmc2, 
                       state=state, 
                       status=TRUE)

# pの分布
plot(fit6$grid,fit6$densp.m,type="l")

# No. clusters の推移と分布
plot(fit6$save.state$thetasave[,1],type="l")
hist(fit6$save.state$thetasave[,1])

# クラスタ数を巨大に
prior9<-list(a0=1,
            b0=1,
            a1=1,
            b1=1)


# Fitting the model
fit9 <- DPbetabinom(y=y,ngrid=100, 
                       prior=prior9, 
                       mcmc=mcmc2, 
                       state=state, 
                       status=TRUE)

# pの分布
plot(fit9$grid,fit9$densp.m,type="l")

# No. clusters の推移と分布
plot(fit9$save.state$thetasave[,1],type="l")
hist(fit9$save.state$thetasave[,1])
hist(fit9$save.state$thetasave[,2])
