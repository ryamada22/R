# n 箇所で連結問題が発生しているとする
# 各箇所 i では、2つの入路pi1,pi2と2つの出路qi1,qi2とがあるものとし
# 組み合わせ問題 (pi1 ->qi1, pi2 ->qi2) vs. (pi1 ->qi2, pi2 ->qi1)の選択　が
# 発生しているものとする
# 組み合わせを、Si = {si1, si2}={(1,2),(2,1)} のように、順列表現として特定することとする

# 今、箇所iにてsijiなる組み合わせのコストを
# f(siji)という関数で計算できるとする
# ただし、関数 f(siji)は、各箇所で採用される組み合わせによって変化するものであるという
# そのことを、f(siji|s1j1,s2j2,...,snjn) と書くことにする

# 計算コストを考えないのであれば、次のようにして最適解が得られる
# すべての箇所のすべての組み合わせ 2^n通りを総当たりで試し
# Sigma f(siji | s1ji,s2j2,...,snjn) が最小になるような組み合わせを解とする

# 今、反復法が有効である、とは
# ある箇所iの組み合わせを
# f(si1 | s1j1,s2j2,...,s(i-1)j(i-1), s(i+1)j(i+1),...,snjn) と
# f(si2 | s1j1,s2j2,...,s(i-1)j(i-1), s(i+1)j(i+1),...,snjn) とを比べ
# 小さい方の組み合わせ sijiを採用する、という手順を繰り返したときに
# 前述の総当たり最適解に収束する
# ということである

# さて。
# これを満足するようなコスト関数 fにはどのような条件があるのだろうか？

n <- 50

p <- matrix(rnorm(n*2),ncol=2)
q <- matrix(rnorm(n*2),ncol=2)

# この関数が考えどころ!
make.f <- function(x,y){
  m <- mean(x-y)
  f <- function(x,y){
    return(abs(abs(x-y) - m))
  }
  return(f)
}

pre.s <- cbind(rep(1,n),rep(2,n))

n.iter <- 100
for(i in 1:n.iter){
  tmp.pre.s <- pre.s
  for(j in 1:n){
    x1 <- p[,1]
    x2 <- p[,2]
    y1 <- q[,1]
    y2 <- q[,2]
    tmp <- which(pre.s[,1]==2)
    y1[tmp] <- q[tmp,2]
    y2[tmp] <- q[tmp,1]
    x. <- c(x1[-j],x2[-j])
    y. <- c(y1[-j],y2[-j])
    tmp.f <- make.f(x.,y.)
    tmp.cost1 <- tmp.f(x1[j],y1[j]) + tmp.f(x2[j],y2[j])
    tmp.cost2 <- tmp.f(x1[j],y2[j]) + tmp.f(x2[j],y1[j])
    if(tmp.cost1 <= tmp.cost2){
      tmp.s <- c(pre.s[j,1],pre.s[j,2])
    }else{
      tmp.s <- c(pre.s[j,2],pre.s[j,1])
    }
    pre.s[j,] <- tmp.s
  }
  print(sum(abs(pre.s-tmp.pre.s)))

}
