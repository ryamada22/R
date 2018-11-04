# 球化したobjファイルには、球面上の座標 x,y,zがある
# それを極座標にする
# 一応書きますが、どこかしらから関数を取ってくる方が安全で賢い
# 処理を具体的にイメージするためのコード
my.kyokuzahyou <- function(v){
	v <- v/sqrt(sum(v^2)) # 念のため長さを１に
	xy <- sqrt(v[1]^2+v[2]^2) # xy 平面への写像の長さ
	x <- v[1]
	y <- v[2]
	z <- v[3]
	theta <- atan(z/xy) # xy==0の場合の場合わけとかを本当はしないといけない
	phi <- atan(y/x)
	return(c(theta,phi))
}
# 少しは検算する
v <- c(1,1,0)
my.kyokuzahyou(v)/pi
v <- c(1/2,-1/2,1/sqrt(2))
my.kyokuzahyou(v)/pi


# theta,phiの三角関数でできた関数を作る
my.f <- function(ang){
	theta <- ang[1]
	phi <- ang[2]
	cos(theta*6)^2 + 3* cos(theta) * sin(phi*4)
}

# 球上の一様乱点は、３次元標準正規乱数から作れる(なぜでしょう？)
my.sprand <- function(n){
	X <- matrix(rnorm(n*3),ncol=3) # 3D正規乱数
	R <- sqrt(apply(X^2,1,sum)) # 各点の原点からの距離
	X <- X/R
	return(X)
}

library(rgl)
X <- my.sprand(10000)
plot3d(X)

ThPhi <- t(apply(X,1,my.kyokuzahyou))
new.r <- apply(ThPhi,1,my.f)

k <- 0.1

new.X <- X * (1+k*new.r)
plot3d(new.X)



