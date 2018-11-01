n.cell <- 20

k <- 10
d <- 3

m <- list()

for(i in 1:n.cell){
	m[[i]] <- matrix(rnorm(k*d),ncol=d)
}


D <- matrix(0,n.cell,n.cell)

for(i in 1:n.cell){
	for(j in 1:n.cell){
		D[i,j] <- sqrt(sum(m[[i]]-m[[j]])^2)
	}
}

image(D)

heatmap(D)

install.packages("tsne")

library(tsne)

tsne.out <- tsne(D)

plot(tsne.out)
