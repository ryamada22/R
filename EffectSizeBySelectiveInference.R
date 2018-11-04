install.packages("selectiveInference")
library(selectiveInference)

# generate a LD block and case-control genotype(G) and phenotype(P) data set
# where n.sample = No. individuals,
#       n.snp = No. snps,
# k SNPs are strongly associated

n.sample <- 1000
n.snp <- 10

G <- matrix(sample(0:2,n.sample*n.snp,replace=TRUE),ncol=n.snp)

as <- rnorm(n.snp,0,0.1)
k <- 2
as[1:k] <- 1:k + 10

R <- G %*% as

P <- (R > mean(R)) + 0

# In GWAS/Genetic epidemiology setting
# type 1 threshold is not determined by forward stepwise nor Lasso
# but it is determined by family wise error rate.
# Actually the spherization method is the method to determine
# the threshold "nominal p value" to reject null hypothesis
# using truncated normal distribution calculator by Botev.
# In the following,
# we just bollow the functions in selectiveInference
# that returns confidence interval of effect size of linear coefficients.

# run forward stepwise

fsfit = fs(x,y)

# compute sequential p-values and confidence intervals

# The type 1 cutoff values in GWAS is conventionally somewhere around 10^(-8)
# As mentioned above the cutoff will be given by Botev method for the set of SNPs.
# Anyway, in the following, confidence intervals and their relation to the cutoff values are evaluated.

alphas <- 10^((-1):(-10))

outs <- list()
for(i in 1:length(alphas)){
	outs[[i]] = fsInf(fsfit,alpha=alphas[i])

}

cis <- matrix(0,length(alphas),2)
for(i in 1:length(alphas)){
	cis[i,] <- outs[[i]]$ci[1,]
}

matplot(log10(alphas),cis,type="l")

# As seen in the plot
# the cofidence interval becomes wider where cutoff p values become lower.
# When nominal p values are 10^(-8) and cutoff p value is 10^(-8), 
# it means its effect is barely accepted as "non-zero",
# i.e., the effect size 1 should be barely out of the confidence interval.
# This corresponds to the relation between cutoff p values and width of confidence interval.

