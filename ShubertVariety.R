n <- 8
k <- 3

# Vは正規直交基底でないとダメらしい
V <- diag(rep(1,n))
position <- sort(sample(1:n,k,replace=FALSE))
position

# シューベルトセルを作りつつ、その条件よりもゆるい
# シューベルト多様体条件に合致する、シューベルトセルも作る
Wcell <- matrix(0,k,n)
Wvar <- matrix(0,k,n)
for(i in 1:k){
	Wcell[i,] <- rnorm(position[i]) %*% V[1:position[i],] 
	t <- sample(c(TRUE,FALSE),n,replace=TRUE,prob=c(0.1,0.8))
	t[1:position[i]] <- TRUE
	t2 <- which(t)
	Wvar[i,] <- rnorm(length(t2)) %*% V[t2,] 
}

Ucell <- my.shubertcell(Wcell)
Uvar <- my.shubertcell(Wvar)

Ucell$U
Uvar$U

Ucell$position
Uvar$position

# シューベルトセルとシューベルト多様体の関係は
# ポジションの不等号が揃っていること(らしい)
Ucell$position <= Uvar$position

for(i in 1:n){
	tmp <- my.unionrank(V[1:i,],Wcell)
	print(tmp$ru)
}
for(i in 1:n){
	tmp <- my.unionrank(V[1:i,],Wvar)
	print(tmp$ru)
}


