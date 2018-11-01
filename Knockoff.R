p=200; n=100; k=15
p=40; n=100; k=15
mu = rep(0,p); Sigma = diag(p)
X = matrix(rnorm(n*p),n)

L2 <- sqrt(apply(X^2,2,sum))
#X <- t(t(X)/L2)
apply(X^2,2,sum)
t(X) %*% X

nonzero = sample(p, k)
beta = 3.5 * (1:p %in% nonzero)
y = X %*% beta + rnorm(n)

# Basic usage with default arguments
result = knockoff.filter(X, y,knockoffs = create.fixed)
print(result$selected)




######
n.iter <- 100
X.mean <- matrix(0,length(X[,1]),length(X[1,]))
for(i in 1:n.iter){
	result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
	X.mean <- result$Xk/n.iter
}



######
str(result)
X <- result$X
X. <- result$Xk

X. <- X.mean

XX <- cbind(result$X,result$X)
XX. <- cbind(result$X,result$Xk)

image(t(XX.) %*% XX.)
image(t(XX) %*% XX)

image(t(XX.) %*% XX. - t(XX) %*% XX)


image(t(X)%*% X - t(X.) %*% X.)
image(t(X)%*% X - t(X.) %*% X)
image(t(X.)%*% X. - t(X.) %*% X)
