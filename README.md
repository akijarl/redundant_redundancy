# redundant_redundancy

Code used to generate simulations of multigenic quantitative trait redundancy, and quantify levels of said redundancy.

Files:
Redund_clean.R - R script containing simulation code written by Katie Lotterhos and code calculating C_chisquared scores written by Sam Yeaman with a wrapper written by Áki Jarl Láruson

Redund_clean.RData - R data file containing output results from four simulations (20 loci representing 2% of sampled genome, 30 loci representing 3% of genome, 40 loci as 4% of genome, and 50 loci as 5% of genome)

RedundGenotype_BinomialCoeff.R - A stand alone script validating the use of the binomial coefficient to quantify genotypes underlying phenotypes in the model from Redund_clean.R