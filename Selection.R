# (n1,...,nK) から (x1,...,xK)が取られる確率

my.choose <- function(n){
	lfactorial(sum(n))-sum(lfactorial(n))
}

my.multi.selection <- function(n,x,log=TRUE){
	if(any(n-x<0)){
		if(log){
			return(-Inf)
		}else{
			return(0)
		}	
	}
	ret <- - my.choose(c(sum(n)-sum(x),sum(x)))
	for(i in 1:length(n)){
		ret <- ret + my.choose(c(n[i]-x[i],x[i]))
	}
	if(log){
		return(ret)
	}else{
		return(exp(ret))
	}
}
my.multi.selection.asympt <- function(n,x,log=TRUE){

	ret <- lfactorial(sum(x)) - sum(lfactorial(x))
	
	ret <- ret + sum(x*log(n))
	if(log){
		return(ret)
	}else{
		return(exp(ret))
	}
}

exp(my.choose(c(3,2)))
my.multi.selection(c(1,2),c(1,0),log=FALSE)

library(gtools)
my.multi.selection.unlabel <- function(n,x,log=FALSE){
	N <- length(n)
	prms <- permutations(N,N)
	x.prms <- matrix(x[prms],ncol=length(prms[1,]))
	x.prms <- unique(x.prms)
	ret <- rep(0,length(x.prms[,1]))
	for(i in 1:length(ret)){
		ret[i] <- my.multi.selection(n,x.prms[i,],log=TRUE)
	}
	ret.all <- log(sum(exp(ret)))
	if(log){
		return(list(p=ret.all,ps=ret))
	}else{
		return(list(p=exp(ret.all),ps=exp(ret)))
	}
}

my.multi.selection.unlabel(c(1,2),c(1,1))

# replace = TRUE バージョン
# 指定したタイプ数 Nで一様ディリクレ乱数を発生し
# それの観察xの下でのunlabel尤度を返す
my.multi.selection.unlabel.asympt <- function(N,x,log=FALSE){
	if(length(x)<N){
		x <- c(x,rep(0,N-length(x)))
	}
	prms <- permutations(N,N)
	x.prms <- matrix(x[prms],ncol=length(prms[1,]))
	x.prms <- unique(x.prms)
	ret <- rep(0,length(x.prms[,1]))
	n <- rdirichlet(1,rep(1,N))
	for(i in 1:length(ret)){
		ret[i] <- my.multi.selection.asympt(n,x.prms[i,],log=TRUE)
	}
	ret.all <- log(sum(exp(ret)))
	if(log){
		return(list(p=ret.all,ps=ret,r=n))
	}else{
		return(list(p=exp(ret.all),ps=exp(ret),r=n))
	}
}

my.multi.selection.unlabel.asympt(1,c(1))

# 一様ディリクレ乱数の、観察xの下でのunlabel尤度を多数のディリクレ乱数について
# 一括して返す
# 計算のオーバヘッド(観察xのunlabel統合)をまとめる(もしくは引数で与える)
my.multi.selection.unlabel.asympt2 <- function(N,x,x.prms=NULL,log=FALSE,n.iter=1000){
	if(length(x)<N){
		x <- c(x,rep(0,N-length(x)))
	}
	
	if(is.null(x.prms)){
		prms <- permutations(N,N)
		x.prms <- matrix(x[prms],ncol=length(prms[1,]))
		x.prms <- unique(x.prms)
	}
	#print(x.prms)
	ns <- rdirichlet(n.iter,rep(1,N))
	tmp <- x.prms %*% t(log(ns))
	tmp <- log(apply(exp(tmp),2,sum))
	tmp2 <- lfactorial(sum(x)) - sum(lfactorial(x))
	tmp <- tmp2 + tmp
	if(log){
		return(list(pr=tmp,r=ns))
	}else{
		return(list(pr=exp(tmp),r=ns))
	}
	
}
out.r <- my.multi.selection.unlabel.asympt2(N=3,x=c(2))

plot(sort(out.r[[1]])/mean(out.r[[1]]))

# さまざまなNについて、観察xの下での、一様ディリクレ乱数の尤度を返す

my.Ns.multi.selection.unlabel.asympt2 <- function(Ns = c(1,2,3),x,x.prms=NULL,log=FALSE,n.iter=1000){
	ret <- list()
	for(i in 1:length(Ns)){
		ret[[i]] <- my.multi.selection.unlabel.asympt2(N=Ns[i],x=x,log=log,n.iter=n.iter)
	}
	return(list(out=ret,n=Ns))
}
out.rs <- my.Ns.multi.selection.unlabel.asympt2(Ns=1:4,x=c(2))

# もしあるNであると決まっていたら、尤度の高いディリクレ乱数だっただろう尤度は高く
# その逆も真
# 結局、尤度の和が1になるように尤度を標準化することで、個々のディリクレ乱数の「割合」が決まる
# 複数のNにわたって、尤度を比較するときには
# Nについて事前確率で重みづけをした上で、上記で行った同一N内の尤度重みづけを考慮した
# 2段階尤度を比較してそれが高い方がより、それらしいことになる

my.like.Ns.multi.selection.unlabel.asympt2 <- function(Ns = c(1,2,3),x,x.prms=NULL,log=FALSE,n.iter=1000,prior=rep(1/length(Ns),length(Ns))){
	out <- my.Ns.multi.selection.unlabel.asympt2(Ns = Ns,x=x,x.prms=x.prms,log=log,n.iter=n.iter)
	for(i in 1:length(out[[1]])){
		out[[1]][[i]][[1]] <- out[[1]][[i]][[1]]/sum(out[[1]][[i]][[1]])*prior[i]
	}
	return(out)
}

outrrr <- my.like.Ns.multi.selection.unlabel.asympt2(Ns=1:5,x=c(2))

sapply(outrrr[[1]],function(x){sum(x[[1]])})

# 検算
# 総数 sum(n)に対して、nの内訳になっているときに
# N個を取り出すことを考える
# N個の内訳がどうなるかは整数分割でパターン分けできる
# どの内訳がどのくらいの確率でおきるかを計算する
# そべの内訳パターンの生起確率を足し合わせると1になることを検算する

library(partitions)
n <- sample(0:4,6,replace=TRUE)
N <- sample(1:(sum(n)-1),1)
prts <- parts(N)
if(length(n)>N){
	prts <- rbind(prts,matrix(0,length(n)-N,length(prts[1,])))
}

prr <- rep(0,length(prts[1,]))
for(i in 1:length(prts[1,])){
	tmp <- my.multi.selection.unlabel(n,prts[,i])
	print(tmp)
	prr[i] <- tmp[[1]]
}
prr
sum(prr)

# ある総数分布 nのときに、N個(x1+...+xk=N)を取り出した時に、どんなNの内訳がおきるかの
# 計算はmy.multi.selection.unlabel(n,x)でできる
# 逆に、xの観察のもとでの、nの尤度も同じ計算
# nらしさはこの計算で行える
############################

library(MCMCpack)
my.random.multi.selection.unlabel <- function(N,x,k){
	R <- rdirichlet(N,rep(1,k))
	ret <- 
}
