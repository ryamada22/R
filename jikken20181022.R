
spout2 <- my.spherization.d(g,w)
k <- 1

ddd <-spout2$dcos
diag(ddd) <- diag(ddd) * 10
outsq <- eigen(abs(ddd)^k*sign(ddd))

outsq <- eigen(abs(spout2$dcos)^k*sign(spout2$dcos))

cosun <- cos(spout2$dst.un)
outsq.un <- eigen(abs(cosun)^k*sign(cosun))



outsq.x <- outsq[[2]][,1:3]
outsq.x.st <- outsq.x/sqrt(apply(outsq.x^2,1,sum))

outsq.xun <- outsq.un[[2]][,1:3]
outsq.x.stun <- outsq.xun/sqrt(apply(outsq.xun^2,1,sum))

plot3d(outsq.x.st)
segments3d(outsq.x.st[c(t(xxx$edge)), ])

sss <- sample(1:length(outsq.x[,1]),10)
for(i in 1:length(sss)){
	spheres3d(outsq.x.st[sss[i],],color=3,radius=0.1)
	spheres3d(outsq.x.st[spout2$taiseki[sss[i]],],color=4,radius=0.1)

}	
vid <- 100
sh <- shortest_paths(g,vid,weights=w)

len <- sapply(sh[[1]],length)

plot3d(outsq.x.st)
spheres3d(outsq.x.st[sh[[1]][[which.max(len)]],],color = "blue",radius=0.1)
segments3d(outsq.x.st[c(t(xxx$edge)), ])


iptaiseki <- rep(0,length(outsq.x.st[,1]))

for(i in 1:length(outsq.x.st[,1])){
	a <- outsq.x.st[i,]
	b <- outsq.x.st[spout2$taiseki[i],]
	iptaiseki[i] <- sum(a*b)
}

