##  First zero on the critical line s = 0.5 + i t
# Rではアドオンパッケージをダウンロードすることがあります
install.packages("pracma")
# うまく行かなければ、https://www.kunihikokaneko.com/free/r/rpackage.html　などを参考にしてください
# ダウンロードしたパッケージを使えるようにするには、以下のようにします
library(pracma)

# [0,20]の細かい数値列を作る
x <- seq(0, 20, len=1001)
head(x)
plot(x)
# 複素数を作る
# 色々なxの値に対して、zも色々作る
z <- 0.5 + x*1i # 1i が虚数単位
z
# 実部と虚部
Re(z)
Im(z)

# zeta()という函数でリーマンのゼータ函数の値を計算する
# その実部と虚部を格納
fr <- Re(zeta(z))
fi <- Im(zeta(z))
# 絶対値
fa <- abs(zeta(z))
# プロットする
plot(x, fa, type="n", xlim = c(0, 20), ylim = c(-1.5, 2.5),
     xlab = "Imaginary part (on critical line)", ylab = "Function value",
     main = "Riemann's Zeta Function along the critical line")
# 線を加える
lines(x, fr, col="blue")
lines(x, fi, col="darkgreen")
lines(x, fa, col = "red", lwd = 2)
# 点を加える
points(14.1347, 0, col = "darkred")
# 説明書きを加える
legend(0, 2.4, c("real part", "imaginary part", "absolute value"),
       lty = 1, lwd = c(1, 1, 2), col = c("blue", "darkgreen", "red"))
# 格子を加える
grid()