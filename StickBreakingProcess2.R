# Dirichlet process
alpha <- 10

n <- 10000
v <- rbeta(n,1,alpha)
hist(v)

ps <- rep(0,n)


logv <- log(v)
log1v <- log(1-v)
cumsumlog1ps <- cumsum(log1v)
ps[1] <- logv[1]
logv. <- logv[-1]
tmp <- logv. + cumsumlog1ps[-n]
ps <- c(logv[1], tmp)
#for(i in 2:n){
#	ps[i] <- logv[i] + sum(log1v[1:(i-1)])
#}
ps <- exp(ps)


hist(ps)
sum(ps)
plot(cumsum(ps[1:100]),type="l")

plot(ps,type="l")

thetas <- rgamma(n,2,4)

par(mfcol=c(1,2))
plot(thetas,ps,type="h")
hist(thetas,col=2,freq=TRUE,density=20)
par(mfcol=c(1,1))


# ƒ‰ƒ“ƒ_ƒ€‚É•ªŠ„‚·‚é

N <- 10000
R <- runif(N-1)
R. <- c(0,sort(R),1)
R.. <- diff(R.)
par(mfcol=c(1,2))
hist(R..)
plot(sort(R..))

x.norm<-rnorm(n=200,m=10,sd=2)
hist(x.norm,main="Histogram of observed data")
plot(density(x.norm),main="Density estimate of data")
plot(ecdf(x.norm),main="Empirical cumulative distribution function")

z.norm<-(x.norm-mean(x.norm))/sd(x.norm) ## standardized data
qqnorm(z.norm) ## drawing the QQplot
abline(0,1) ## drawing a 45-degree reference line
x.wei<-rweibull(n=200,shape=2.1,scale=1.1) ## sampling from a Weibull distribution with parameters shape=2.1 and scale=1.1
x.teo<-rweibull(n=200,shape=2, scale=1) ## theorical quantiles from a Weibull population with known paramters shape=2 e scale=1
qqplot(x.teo,x.wei,main="QQ-plot distr. Weibull") ## QQ-plot
abline(0,1) ## a 45-degree reference line is plotted

x.poi<-rpois(n=200,lambda=2.5)
hist(x.poi,main="Poisson distribution")
curve(dnorm(x,m=10,sd=2),from=0,to=20,main="Normal distribution")
curve(dgamma(x, scale=1.5, shape=2),from=0, to=15, main="Gamma
distribution")
curve(dweibull(x, scale=2.5, shape=1.5),from=0, to=15, main="Weibull
distribution")

library(fBasics) ## package loading
skewness(x.norm) ## skewness of a normal distribution

kurtosis(x.norm) ## kurtosis of a normal distribution

skewness(x.wei) ## skewness of a Weibull distribution

kurtosis(x.wei) ## kurtosis of a Weibull distribution

mean.hat<-mean(x.norm)
mean.hat


x.gam<-rgamma(200,rate=0.5,shape=3.5) ## sampling from a gamma distribution with ƒÉ=0.5 (scale parameter 12) and ƒ¿=3.5 (shape parameter)
med.gam<-mean(x.gam) ## sample mean
var.gam<-var(x.gam) ## sample variance
l.est<-med.gam/var.gam ## lambda estimate (corresponds to rate)
a.est<-((med.gam)^2)/var.gam ## alfa estimate

l.est
a.est

library(stats4) ## loading package stats4
ll<-function(lambda,alfa) {n<-200
x<-x.gam
-n*alfa*log(lambda)+n*log(gamma(alfa))-(alfa-1)*sum(log(x))+lambda*sum(x)}
## -log-likelihood function
est<-mle(minuslog=ll, start=list(lambda=2,alfa=1))
summary(est)

library(MASS) ## loading package MASS
fitdistr(x.gam,"gamma") ## fitting gamma pdf parameters

fitdistr(x.wei,densfun=dweibull,start=list(scale=1,shape=2))## fittingWeibull pdf parameters

fitdistr(x.norm,"normal") ## fitting gaussian pdf parameters
lambda.est<-mean(x.poi) ## estimate of parameter lambda
tab.os<-table(x.poi)## table with empirical frequencies
freq.os<-vector()
for(i in 1: length(tab.os)) freq.os[i]<-tab.os[[i]] ## vector of emprical frequencies
freq.ex<-(dpois(0:max(x.poi),lambda=lambda.est)*200) ## vector of fitted (expected) frequencies

freq.os
freq.ex

acc<-mean(abs(freq.os-trunc(freq.ex))) ## absolute goodness of fit index
acc
acc/mean(freq.os)*100 ## relative (percent) goodness of fit index


h<-hist(x.norm,breaks=15)
xhist<-c(min(h$breaks),h$breaks)
yhist<-c(0,h$density,0)
xfit<-seq(min(x.norm),max(x.norm),length=40)
yfit<-dnorm(xfit,mean=mean(x.norm),sd=sd(x.norm))
plot(xhist,yhist,type="s",ylim=c(0,max(yhist,yfit)), main=hNormal pdf and
histogramh)
lines(xfit,yfit, col=hredh)


library(vcd)## loading vcd package
gf<-goodfit(x.poi,type= "poisson",method= "MinChisq")
summary(gf)

plot(gf,main="Count data vs Poisson distribution")

x.gam.cut<-cut(x.gam,breaks=c(0,3,6,9,12,18)) ##binning data
table(x.gam.cut) ## binned data table

x.gam.cut

(pgamma(3,shape=a.est,rate=l.est)-pgamma(0,shape=a.est,rate=l.est))*200
(pgamma(6,shape=a.est,rate=l.est)-pgamma(3,shape=a.est,rate=l.est))*200
(pgamma(9,shape=a.est,rate=l.est)-pgamma(6,shape=a.est,rate=l.est))*200
(pgamma(12,shape=a.est,rate=l.est)-pgamma(9,shape=a.est,rate=l.est))*200
(pgamma(18,shape=a.est,rate=l.est)-pgamma(12,shape=a.est,rate=l.est))*200
f.ex<-c(20,71,61,31,17) ## expected frequencies vector


f.os<-vector()
for(i in 1:5) f.os[i]<- table(x.gam.cut)[[i]] ## empirical frequencies
X2<-sum(((f.os-f.ex)^2)/f.ex) ## chi-square statistic
gdl<-5-2-1 ## degrees of freedom
1-pchisq(X2,gdl) ## p-value


## computing relative expected frequencies
p<-c((pgamma(3,shape=3.5,rate=0.5)-pgamma(0,shape=3.5,rate=0.5)),
(pgamma(6,shape=3.5,rate=0.5)-pgamma(3,shape=3.5,rate=0.5)),
(pgamma(9,shape=3.5,rate=0.5)-pgamma(6,shape=3.5,rate=0.5)),
(pgamma(12,shape=3.5,rate=0.5)-pgamma(9,shape=3.5,rate=0.5)),
(pgamma(18,shape=3.5,rate=0.5)-pgamma(12,shape=3.5,rate=0.5)))
chisq.test(x=f.os,p=p) ## chi-square test

ks.test(x.wei,"pweibull", shape=2,scale=1)
x<-seq(0,2,0.1)
plot(x,pweibull(x,scale=1,shape=2),type="l",col="red", main="ECDF and
Weibull CDF")
plot(ecdf(x.wei),add=TRUE)

shapiro.test(x.norm)
library(tseries) ## package tseries loading
jarque.bera.test(x.norm)
zz<-rnorm(n=200,m=0,sd=1) ## sampling random numbers from N(0,1)
r<-zz[200]
q<-sd(zz[-200])
m<-mean(x.norm)
s<-sqrt(var(x.norm))
y<-q*((x.norm-m)/s)+(r/sqrt(200))
ks.test(y,hpnormh,m=0,sd=1)
library(nortest) ## package loading
sf.test(x.norm)
library(nortest) ## package loading
ad.test(x.norm)

library(nortest) ## package loading
cvm.test(x.norm)
library(nortest) ## package loading
lillie.test(x.norm)
library(nortest) ## package loading
pearson.test(x.norm)

