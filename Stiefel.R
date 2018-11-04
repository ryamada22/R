n <- 5
p <- 3
A <- matrix(rnorm(n*p),ncol=p)
A
# qr 分解
qrout <- qr(A)
Q <- qr.Q(qrout)
R <- qr.R(qrout)
Q
R
# A = QR
Q %*% R - A

# t(R) %*% R はAのペアワイズ距離二乗行列
t(R) %*% R
t(A) %*% A

# 結局
# A=QRと分解できて、Rは、Aのペアワイズ距離を保ったp次元座標
# Qは nxp行列であって、これは、p本のn次元ベクトルで、それらが正規直交：Stiefel行列
t(Q) %*% Q