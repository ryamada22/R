# 集合Sがあったとき
# Sの部分集合の集合Tがある
# 今、Sの別のタイプの部分集合Kは必ずTのただ一つ要素の部分集合である
# Kの全要素をTに帰属させたい
# そのときKのあつ要素を調べることで、それが帰属するTの要素が確認できると言う
# またTが決まると、それに帰属するKの要素は列挙可能である

# まずペアをパスに帰属させる…

# Kのリストがあったときに、K-T帰属関係をどのように探索するのがよいだろうか？

my.STK <- function(g,w,k=1){
	nv <- length(V(g)) # 頂点数
	d <- distances(g,weights=w) # エッジ重み付き頂点間距離
	ad <- get.adjacency(g) # 隣接行列
	d.noWeight <- distances(g) # エッジ重みなし頂点間距離
	#ret <- matrix(0,nv,nv) # i番ノードにとってj番ノードが極遠点なら1を立てる行列
	sh <- my.WarshallFloyd(g,w) # Warshall-Floyd 行列を作っておく
	trios <- expand.grid(1:nv,1:nv,1:nv)

	loop <- TRUE

	while(loop){
		this.trio <- trios[1,]
		shout12 <- extractPath(sh,this.trio[1],this.trio[2])
		shout23 <- extractPath(sh,this.trio[2],this.trio[3])
		shout31 <- extractPath(sh,this.trio[3],this.trio[1])
		full <- c(shout12,shout23,shout31)
		un <- unique(full)
		if(length(full) == length(un)+3){
			tmp.cycle <- 
		}else{
			trios <- trios[-1,]
		}
		if(length(trios[,1])==0){
			loop <- FALSE
		}
	}
	
}