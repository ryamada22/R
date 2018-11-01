n <- 5
k <- 2
W <- matrix(rnorm(n*k)+rnorm(n*k)*1i,k,n)
W # n次元複素ベクトルのk本の組
Wk <- W[1:k,1:k] # 先頭のk列が作る正方行列
Wkinv <- solve(Wk) # その逆行列
Wst <- Wkinv %*% W # 標準化する
Wst # 先頭のk列部分が単位行列化している。残りのn-k列部分が自由度のある部分
library(complexplus)
library(gtools)
cmb <- combinations(n,k)
pl <- pl.st <- rep(0,length(cmb[,1]))
for(i in 1:length(pl)){
	pl[i] <- Det(W[,cmb[i,]])
	pl.st[i] <- Det(Wst[,cmb[i,]])
}
pl
pl.st # 1,2列で標準化したWstのプリュッカー座標なので、第１成分が1

pl / pl.st # WとWstは同一視できるものなので、そのプリュッカー座標も射影座標としては同一
pl/pl[1] - pl.st # プリュッカー座標にしてから標準化してもよい

# n列から、k-1列をとるすべての組み合わせ
cmb_ <- combinations(n,k-1)
# n列から、k+1列をとるすべての組み合わせ
cmb. <- combinations(n,k+1)

# k-1列の組み合わせと、k+1列の組み合わせの総当り
# どの場合にも、小行列に関する式の値が0になることを以下で確認する
for(i in 1:length(cmb_[,1])){
	for(j in 1:length(cmb.[,1])){
		print("#####")
		# k-1列の組み合わせID
		print(paste("i=",i))
		# k+1列の組み合わせID
		print(paste("j=",j))
		a <- cmb_[i,]
		b <- cmb.[j,]
		# k-1列、k+1列の列番号を提示
		print(paste("k-1 elems=",toString(sort(a))))
		print(paste("k-2 elems=",toString(sort(b))))
		# 加算値の初期値
		tmp <- 0
		# 加算式の初期オブジェクト
		tmpformula <- ""
		print("======")
		# k+1列の要素を1番からk+1番まで順々に取り出してk-1側に移行しつつ積算する
		for(p in 1:(k+1)){
			jL <- b[p] # k+1要素のp番目要素
			a. <- c(a,jL) # k-1側の末尾に加える(行列式では行番・列番は大事なので、末尾に)
			b. <- b[-p] # k+1側から抜き出す
			W1 <- Det(W[,a.]) # k-1側からできたk列の小行列式
			W2 <- Det(W[,b.]) # k+1側からできたk列の小行列式
			print("----")
			# k-1側、k+1側からできたk要素を並べる。同一要素があるかどうかも確認できる
			print(paste("k elems=",toString(sort(a.))))
			#print(W1)
			print(paste("k elems=",toString(sort(b.))))
			#print(W2)
			# 加算にあたっては、重複列があるためにW1=0となっているかもしれないが
			# 気にせず加算する。0を加えても変化なし
			tmp <- tmp + (-1)^p * W1 * W2
			print("----")
			# 加算式の文字列オブジェクトを作る場合には
			# 0となる小行列式が入らないように書く
			# こうすることでプリュッカー座標に登場しない列の組は加算式に現れなくなる
			if(W1!=0 & W2!=0){
				# 加算の際の符号のための文字列
				ss <- "+"
				if((-1)^p == (-1)){
					ss <- "-"
				}
				# k列の表示は、a.,b.を作ったとおりにする
				# 列順が変わると小行列式の符号が変わってしまうから
				#tmpformula <- paste(tmpformula,ss,"W", toString(sort(a.)), "W", toString(sort(b.)))
				tmpformula <- paste(tmpformula,ss,"W", toString(a.), "W", toString(b.))
				#print(tmpformula)
			}
		}
		#print(tmp)
		print(paste(tmpformula, "=", tmp))
		print("#####")
	}
}

