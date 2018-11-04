n.paper <- 20 # 論文数
n.citation <- sample(1:100,20) # 引用数
n.cite.dec <- sort(n.citation,decreasing=TRUE)
plot(n.cite.dec)

tmp <- n.cite.dec - 1:n.paper
plot(1:n.paper,tmp)
abline(h=0)
tmp.nonneg <- which(tmp >= 0)
h.index <- tmp.nonneg[length(tmp.nonneg)]
abline(v=h.index)

library(scholar)
myid <- 'VG3xe_IAAAAJ'
pb <- get_publications(myid)
str(pb)
n.citation <- pb$cites
n.paper <- length(n.citation)
n.cite.dec <- sort(n.citation,decreasing=TRUE)
tmp <- pb$cites - 1:n.paper
plot(1:n.paper,tmp)
abline(h=0)
tmp.nonneg <- which(tmp >= 0) # 差が0以上の論文のランク列
this.paper <- length(tmp.nonneg) # この論文のランク・引用数がh-indexを決める
# 当該論文のランクか、その引用数の小さい方
h.index <- min(this.paper,n.cite.dec[this.paper])
abline(v=h.index)

my.h.index <- function(v){
	sort.v <- sort(v,decreasing=TRUE)
	tmp <- sort.v - 1:length(sort.v)
	tmp.nonneg <- which(tmp>=0)
	this.paper <- length(tmp.nonneg) # この論文のランク・引用数がh-indexを決める
	# 当該論文のランクか、その引用数の小さい方
	h.index <- min(this.paper,n.cite.dec[this.paper])
	return(h.index)
}

year <- 2005:2018
h.index.year <- rep(0,length(year))
for(i in 1:length(year)){
	n.citation.tmp <- pb$cites[which(pb$year<=year[i])]
	h.index.year[i] <- unlist(my.h.index(n.citation.tmp))

}

plot(h.index.year,ylim=c(0,max(h.index.year)))

future <- predict_h_index(myid)

plot(future,ylim=c(0,max(future)))

# 過去と未来
pastFuture <- c(h.index.year,future$h_index[-1])
plot(pastFuture,ylim = c(0,max(pastFuture)))
abline(v=length(h.index.year))
