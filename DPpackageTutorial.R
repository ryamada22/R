data(galaxy)
galaxy<-data.frame(galaxy,speeds=galaxy$speed/1000)
attach(galaxy)
# Initial state
state <- NULL
# MCMC parameters
nburn<-100
nsave<-100
nskip<-10
ndisplay<-100
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)
# Prior
prior<-list(aa0=2.01,
ab0=0.01,
kmax=10,
a0=1,
b0=1)
# Fitting the model
fit <- BDPdensity(y=speeds,prior=prior,mcmc=mcmc,
state=state,status=TRUE)
plot(fit)
