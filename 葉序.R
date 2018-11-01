# 黄金比
psi <- (1+sqrt(5))/2
# 黄金角
psi.angle <- (2*pi)/(1+psi)
# シリアルナンバー
t <- 0:5
# シリアル角
angles <- psi.angle * t
# 角は y = ax 的に増える
plot(t,angles)
# 角に対応する単位円周上の座標
x <- cos(angles)
y <- sin(angles)
# プロットする
plot(x,y)
# 順番を見せる
plot(x,y,type="b")
# シリアルナンバーを増やす
t <- 0:100

angles <- psi.angle * t

plot(t,angles)

x <- cos(angles)
y <- sin(angles)

plot(x,y)

plot(x,y,type="b")

# 回りながら半径を増やす
t <- 0:5
r <- t
angles <- psi.angle * t

plot(t,angles)

x <- cos(angles)
y <- sin(angles)

x2 <- x * r
y2 <- y * r
plot(x2,y2,type="b")

# 回りながら半径を増やす
t <- 0:100
r <- t
angles <- psi.angle * t

plot(t,angles)

x <- cos(angles)
y <- sin(angles)

x2 <- x * r
y2 <- y * r
plot(x2,y2,type="b")

library(rgl)
k <- 1
z <- t * k

plot3d(x2,y2,z)
