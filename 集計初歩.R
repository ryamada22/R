# cell number
n <- 1000
# SH decomposition coefficient matrix (l+1)^2 * 3
l <- 7

sh.list <- list()
for(i in 1:n){
	sh.list[[i]] <- matrix(rnorm((l+1)^2*3),ncol=3)
}

sh.arr <- array(0,c(n,(l+1)^2,3))
for(i in 1:n){
	sh.arr[i,,] <- sh.list[[i]]
}

# cellごとの属性(たとえば体積、表面積)を列に格納したもの
# カテゴリーも数値で表すことにして
# 数値行列扱いする
# １行目：体積
# ２行目：表面積
# ３行目：エラーメッセージタイプ {0,1,2,...}
# ４行目：実験ID
# ５行目：時刻
# ６行目：頂点数

n.atr <- 6
cell.atr <- matrix(0,n,n.atr)

cell.atr[,1] <- runif(n)
cell.atr[,2] <- runif(n)
cell.atr[,3] <- sample(0:2,n,replace=TRUE,prob=c(0.9,0.05,0.05))
cell.atr[,4] <- sample(1:10,n,replace=TRUE)
cell.atr[,5] <- sample(1:20,n,replace=TRUE)
cell.atr[,6] <- sample(1000:5000,n)

# 単純に集計する

# カテゴリは集計する

table(cell.atr[,3])

# 量的変数はヒストグラム
hist(cell.atr[,1])

# たくさんあってヒストグラムがやっかいなときは
# 値をソートしてプロットすれば
# その単調増加曲線が分布を表すので
# その単調増加曲線を合わせて描く

sh.mat <- matrix(0,n,(l+1)^2*3)
for(i in 1:n){
	sh.mat[i,] <- sh.list[[i]]
}
matplot(apply(sh.mat,2,sort),type="l")

# 群別ヒストグラム
# https://stats.biopapyrus.jp/r/ggplot/geom_histogram.html 

# 群別カテゴリ集計
table(cell.atr[,3],cell.atr[,4])

# 量的変数のペアワイズ関係
quant.col <- c(1,2,6)
plot(as.data.frame(cell.atr[,quant.col]))

# 全部の変数でとにかく関係を計算しておく

sh.atr.mat <- cbind(sh.mat,cell.atr)
# cor.mat

cor.mat <- cor(sh.atr.mat)
image(cor.mat)


# エラーメッセージで選別して作っておく

non.errors <- which(cell.atr[,3]==0)

sh.arr.ne <- sh.arr[non.errors,,]
cell.atr.ne <- cell.atr[non.errors,]

#####
# 回転同一視
# いい加減に回転同一視での、細胞ペア間内積が取れたとする。
ip.cells <- matrix(1,n,n) # 対角成分の内積は1
for(i in 1:(n-1)){
	for(j in i:n){
		ip.cells[i,j] <- ip.cells[j,i] <- runif(1) * sample(c(-1,1),1)
	}
}
image(ip.cells)

hclust.out <- hclust(as.dist(acos(ip.cells)))
plot(hclust.out)

# 対称行列は固有値分解
eigen.out <- eigen(ip.cells)
plot(eigen.out[[1]])

