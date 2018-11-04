# 直方体ボクセルを作る(ボクセルの１頂点座標でボクセルを表させる)

my.rectVoxel <- function(n,st = rep(0,3)){
	ns <- rep(n,3)
	ret <- as.matrix(expand.grid(0:ns[1],0:ns[2],0:ns[3]))
	ret <- t(t(ret) + st)
	ret
}

# ２つの連結していない直方体ボクセル集合を併せて、１が立っているボクセルを列挙

Vox.list <- unique(rbind(my.rectVoxel(3),my.rectVoxel(3,c(5,5,5))))

# グラフ理論パッケージ
library(igraph)
# ボクセル座標間の位置を計算
d <- as.matrix(dist(Vox.list))
# 距離が１のボクセル同士にエッジを置く
el <- which(d==1,arr.ind=TRUE)
# エッジリストからグラフを作る
g <- graph.edgelist(el)
# グラフをプロット
plot(g)
# 連結成分取り出し

clusters(g)