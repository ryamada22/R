# Yamada's observation of punctuality
# When you plan to be in time for something,
# you can not realize the probability to be late zero.
# Yamada knows it as a statistician.
# Then how does he observe others' (and his own) sicerity not to be late?

my.sincere.punc <- function(x){
	m <- mean(x)
	v <- var(x)
	p <- pnorm(9,m,sqrt(v),lower.tail=FALSE)
	return(p)
}

x.yamada <- c(8,8.25,7.75,8,8.15,7.85)

yamada.sincere <- my.sincere.punc(x.yamada)
yamada.sincere

x.hoge <- c(8.92,8.93,8.95,8.9,8.93,8.99)

hoge.sincere <- my.sincere.punc(x.hoge)

hoge.sincere

x.hoge2 <- c(8.92,8.93,8.95,8.9,8.93,9.01)

hoge2.sincere <- my.sincere.punc(x.hoge2)

hoge2.sincere
