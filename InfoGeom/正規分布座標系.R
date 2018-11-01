# 正規分布、(m,s)座標に(theta1,theta2)格子
# s^2 = m/theta1
# s^2 = -1/(2theta2)
theta1 <- seq(from=1/2,to=1000,by=1/2)
theta2 <- -seq(from=1/2,to=1000,by=1/2)
fr <- matrix(c(-1,0,1,1),byrow=TRUE,2,2)
plot(fr,asp=1,col=0,ylim=c(0,2))

for(i in 1:length(theta1)){
	#abline(0,theta1[i])
	#abline(0,-theta1[i])
	segments(0,0,2,2/theta1[i])
	segments(0,0,-2,2/theta1[i])
}

for(i in 1:length(theta2)){
	abline(h=sqrt(-1/(2*theta2[i])))
}
abline(h=0)

# 正規分布、(m,s^2)座標に(eta1,eta2)格子
# m = eta1
# m^2+s^2=eta2
eta1 <- seq(from=-1,to=1,by=1/(2^4))
eta2 <- seq(from=0,to=2,by=1/(2^4))

fr <- matrix(c(-1,0,1,1),byrow=TRUE,2,2)
plot(fr,asp=1,col=0,ylim=c(0,2))

for(i in 1:length(eta1)){
	abline(v=eta1[i])
}
t <- seq(from=0,to=1,length=100) * 2*pi
for(i in 1:length(eta2)){
	points(eta2[i] * cos(t),eta2[i] * sin(t),type="l")
}

##

fr <- matrix(c(-0.5,0,0.5,0.5),byrow=TRUE,2,2)
plot(fr,asp=1,col=0,ylim=c(0,1))
for(i in 1:length(theta1)){
	#abline(0,theta1[i])
	#abline(0,-theta1[i])
	segments(0,0,2,2/theta1[i],col=2)
	segments(0,0,-2,2/theta1[i],col=2)
}

for(i in 1:length(theta2)){
	abline(h=sqrt(-1/(2*theta2[i])),col=2)
}
for(i in 1:length(eta1)){
	abline(v=eta1[i],col=3)
}
t <- seq(from=0,to=1,length=100) * 2*pi
for(i in 1:length(eta2)){
	points(eta2[i] * cos(t),eta2[i] * sin(t),type="l",col=3)
}
