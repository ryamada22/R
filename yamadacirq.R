# Yamada-cirq

my.cirq <- function(x,N=1000){
	n <- length(x)
	theta <- (0:(n-1))/n * 2*pi
	X <- cbind(cos(theta),sin(theta))
	d <- as.matrix(dist(x))
	diag(d) <- 1
	pairs <- which(d==0,arr.ind=TRUE)
	pairs <- unique(t(apply(pairs,1,sort)))

	alphas <- theta[pairs[,1]]
	betas <- theta[pairs[,2]]
	
	rs <- 1/cos((betas-alphas)/2)
	
	ctrs <- rs * cbind(cos((alphas+betas)/2),sin((alphas+betas)/2))
	
	Rs <- tan((betas-alphas)/2)
	
	phis <- cbind(betas+pi/2, alphas+3*pi/2)
	Theta <- (0:(N-1))/N * 2 * pi
	plot(cos(Theta),sin(Theta),type="l")
	points(X,pch=20)
	for(i in 1:length(pairs[,1])){
		Phi <- seq(from = phis[i,1],to=phis[i,2],length=N)
		arc <- Rs[i] * rbind(cos(Phi),sin(Phi)) + ctrs[i,]
		points(t(arc),type="l")
	}
	
}

x <- c(0,4,1,2,10,3,10,1,4,5,6,7,8,4,9)

my.cirq(x)
