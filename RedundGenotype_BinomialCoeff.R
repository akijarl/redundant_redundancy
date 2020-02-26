## redundancy.R
## Katie E Lotterhos
## Feb 12, 2019
## Illustrate concept of genetic redundancy for different phenotypic values

library(gtools)
library(dplyr)

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
  
  genotype <- c(genotype, apply(combo.df2 , 1 , paste , collapse = "_" ))
  phenotype <- c(phenotype, rowSums(combo.df))
}

final.df <- data.frame(genotype, phenotype)
head(final.df, 100)
tail(final.df, 100)

# Plot example of fitness as a function of phenotype ####
# stabilizing selection

#pdf("RedundancyForOptimum_10genotype.pdf", width=5, height=4)  
par(mar=c(4,3,1,1))

x2<- seq(0,2,0.001)
y2 <- dnorm(x2, mean = 1, sd=0.5/10)
y2_scaled <- y2/max(y2)

plot(c(x2, x2[1]), c(y2_scaled, y2_scaled[1]), 
     ylim=c(-1,1.3), bty="n", yaxt="n",type="l",
     xlim=c(-2,2), xaxs="i", ylab="", xlab="Phenotype",
     col=rgb(0,0.4, 0.8))

polygon(c(x2, x2[1]), c(y2_scaled, y2_scaled[1]), 
        col=rgb(0,0.4, 0.8))


mtext("Density", side=2, line=1)

text(0, 0.5, "Fitness\nfunction", col=rgb(0,0.4, 0.8))
abline(h=0, lwd=1)

# Plot redundancy as a function of phenotype
# number of ways to get phenotype over evolutionary timescales
y1 <- hist(phenotype, breaks=seq(-alpha*nn-alpha/2,alpha*nn+alpha/2, alpha), plot=FALSE)
(y1df <- cbind(y1$mids, y1$counts))
y1_scaled =  y1$density/max(y1$density)
points(y1$mids, -y1_scaled, type="h", col=rgb(1,0.5,0.5), lwd=1)

text(0.7, -0.5, "Genetic\nredundancy", col=rgb(1,0.5,0.5))

#dev.off()
# end plot ####

y1df

redundancy <- function(nloci, alpha) {
  # assume equal number of + loci and - loci
  Num <- choose(nloci, 0:nloci) # probability surface for redundancy
  phen <- round(seq(-alpha*nloci/2, alpha*nloci/2, alpha),3)
  return(data.frame(phen, Num, log10Num=log10(Num)))
}

redundancy(10, 0.2)
redundancy(40, 0.2)
