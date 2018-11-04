# 2x2分割表のカイ二乗値はサンプル総数に比例する

x <- runif(4)
t1 <- matrix(x,2,2)

# chisq.test()のデフォルトは、イェーツ補正をするようになっており
# それだと、カイ二乗値とサンプル数との関連が現れないので、correct = FALSEとしています
# また、遺伝統計では、イェーツの補正のように、「１個の検定」のための(賢しらな)補正法などを使うと、たくさんの出力の分布や性質がゆがむので、素のカイ二乗統計量を使います。

out1 <- chisq.test(t1,correct=FALSE)

k <- runif(1)
t2 <- k * t1
out2 <- chisq.test(t2,correct=FALSE)

out1$stat
out2$stat

out2$stat/out1$stat
k
