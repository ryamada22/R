n <- 10 # •ª•z‚Ì”
d <- 2 # €”
library(MCMCpack)
R <- rdirichlet(n,rep(0.6,d)) 

H <- matrix(0,n,n)

for(i in 1:n){
	for(j in 1:n){
		H[i,j] <- log(sum(R[i,]*R[j,]))/2
	}
}
H. <- -H
eigen.out <- eigen(H.)

V <- eigen.out[[2]] + 0*1i
S <- diag((eigen.out[[1]] + 0*1i)^0.5)

W <- V %*% S

range(Re(W%*% t(W)) -H.)
range(Im(W%*% t(W)) -H.)

Psi <- apply(W^2,1,sum)

P. <- log(t(R)) + Psi
Theta. <- t(cbind(W,rep(1,n)))

F. <- P.%*% t(Theta.) %*% solve(Theta.%*% t(Theta.))

log(t(R)) - (F. %*% Theta. - Psi)