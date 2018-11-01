N <- 10 # number of advisers

p <- runif(N) # probability of advisers to predict answer right

k <- 10000 # number of trials

truths <- sample(0:1,k,replace=TRUE)

w <- matrix(0,k,N)
w[1,] <- 1

for(i in 2:k){

}
