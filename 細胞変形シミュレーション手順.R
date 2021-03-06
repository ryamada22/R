# パッケージの準備
library(devtools)
# Ronlyryamadaというパッケージをインストールするのに、devtoolsというパッケージのインストールが必要です
install_github("ryamada22/Ronlyryamada") 
library(Ronlyryamada)
library(rgl)
library(RFOC)

# 使うのはRonlyryamada パッケージ
help(package="Ronlyryamada")
# 使うのは、そのmy.cell.sim()関数
# まずは使ってみる
help(my.cell.sim)
# たくさんのファイルができるので、作業場所を選んでおく
# setwd()関数を使って作業場所を指定する
#examples をコピーペーストして実行

# 50枚の３次元プロットができる。そのプロットをさせるファイルを50ファイル作成する

# 作業の内容
# my.cell.sim()には色々な変数があります。それらの取り方によっては、３次元プロットが、よくない形を作ります
# よい形というのは、細胞の形のように、でこぼこはするものの、全体が閉じていて、内側にめり込んでいたりしないことをいいます。
# パラメタを色々と変えて実行し、「よい形になった場合」がどれで、駄目な形がどれかを分けること
# setwd()を使うと、出力先を変更できるので、パラメタを変えるたびに、出力先のフォルダ名を変えて、「よい形」なら、「よい形フォルダ」の下にそれを移し、「悪い形」なら、「悪い形フォルダ」に移します
# 各回の実行パラメタもhoge_param.txtというファイル名で出ているはずなので、それも併せて保管しておきます
# もし、うまく行ったパラメタセットとうまく行かないパラメタセットをエクセルなどに管理して保存してもらえれば、さらにうれしいです
