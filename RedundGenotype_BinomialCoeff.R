## redundancy.R
## Katie E Lotterhos
## Feb 12, 2019
## Illustrate concept of genetic redundancy for different phenotypic values

library(gtools)
library(dplyr)

### Function to calculate Redundancy for `n=10` genotypes - NUMERICALLY ####

redundancy.numerical <- function(nloci, alpha) {
		
	np <- round(nloci/2,0) #number of mutations with positive effect sizes
	nn <- np #number of mutations with negative effect sizes
	
	p_names <- paste0("p", 1:np)
	n_names <- paste0("n", 1:nn)
	
	all_loci <- data.frame(names=(c(p_names, n_names)), 
	                       alpha = c(rep(alpha, length(p_names)),
	                                 rep(-alpha, length(n_names))))
	all_loci$names <- as.character(all_loci$names)
	
	# could have 1 to 10 mutations evolve
	# need to loop through all possible mutational states "r" from 1 to 10
	# each one is a different size matrix
	# idea: collapse genotype into single character with it's phenotype value
	
	genotype <- c()
	phenotype <- c()
	
	for (r in 1:nrow(all_loci)){
	  print(r)
	  combo.df <- data.frame(combinations(n=length(all_loci$names), r=r, 
	                           all_loci$names, repeats.allowed = FALSE),
	                         stringsAsFactors = FALSE)
	  
	  head(combo.df)
	  colnames(combo.df) <- paste0("A", 1:ncol(combo.df))
	  
	  combo.df2 <- combo.df
	  for (i in 1:ncol(combo.df)){
	    combo.df[,i][grep("n", combo.df[,i], perl=TRUE)] <- -alpha
	      # get indexes with n and replace with minus alpha
	    combo.df[,i][grep("p", combo.df[,i], perl=TRUE)] <- alpha
	    combo.df[,i] <- as.numeric(combo.df[,i])
	  }
	  
	  
	  genotype <- c(genotype, apply( combo.df2 , 1 , paste , collapse = "_" ))
	  phenotype <- c(phenotype, rowSums(combo.df))
	}
	final.df <- data.frame(genotype, phenotype)
	return(final.df)
}

### Function to calculate Redundancy for `n=10` genotypes - BINOMIAL COEFFICIENT ####

redundancy <- function(nloci, alpha) {
  # assume equal number of + loci and - loci
  Num <- choose(nloci, 0:nloci) # probability surface for redundancy
  phen <- round(seq(-alpha*nloci/2, alpha*nloci/2, alpha),3)
  return(data.frame(phen, Num, log10Num=log10(Num)))
}


alpha = 0.2
n=10
final.df.num <- redundancy.numerical(n, alpha)
hist(final.df.num$phenotype, main = "Numerical number of genotypes for each phenotype", xlab="Phenotype", breaks=seq(-(n/2)*alpha-alpha/2, (n/2)*alpha+alpha/2, by=alpha))


final.df.binom <- redundancy(n, alpha)
points(final.df.binom$phen, final.df.binom$Num, col="blue")

alpha = 0.3
n=20
final.df.num <- redundancy.numerical(n, alpha)
hist(final.df.num$phenotype, main = "Numerical number of genotypes for each phenotype", xlab="Phenotype", breaks=seq(-(n/2)*alpha-alpha/2, (n/2)*alpha+alpha/2, by=alpha))


final.df.binom <- redundancy(n, alpha)
points(final.df.binom$phen, final.df.binom$Num, col="blue")

