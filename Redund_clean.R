## Originally: redundancy.R
## Katie E Lotterhos
## Feb 12, 2019
## Illustrate concept of genetic redundancy for different phenotypic values
## Mod. by Áki Jarl 
## Added C_chisq calculation with function from Sam Yeaman and visualization
## Aug. 28, 2019
################################################################
# Calculate C_chisq and show how it changes with number of loci 
################################################################
library(gtools)
library(dplyr)

# Don't go below 20 loci
redundancy <- function(nloci, alpha) {
  # assume equal number of + loci and - loci
  Num <- choose(nloci, 0:nloci) # probability surface for redundancy
  phen <- round(seq(-alpha*nloci/2, alpha*nloci/2, alpha),3)
  return(data.frame(phen, Num, log10Num=log10(Num)))
}


w.zik <- function(zik, theta, omega.sq){
  # Used in simulation
  # zik is phenotype of individual i in patch k
  # theta is the optimum in each patch k
  # omega.sq is the strength of stabilizing selection
  1 - (zik - theta)^2/omega.sq
}

# simulatePop modified by ÁJL
simulatePop <- function(N, mu, alpha, nloci, omega.sq, m, seed=NULL){
  opt <- c(-1, 1) 
  burnin <- 4*N
  ### loop through generations ####
  for (t in 1:(8*N)){
    
    if(t==1){ # initialize
      pop <- matrix(0, N, nloci)
      tot <- N*nloci
      patch <- c(rep(-1, N/2), rep(1, N/2))
      if (length(seed)>0){
        i <- seed
      }else{
        i <- round(runif(1)*10e06, 0)
        print(paste("seed = ", i))
      }
    }
    
    if (t %% 1000 == 0){print(paste("Sim. Generation", t))}
    # Mutation
    set.seed(i)
    a<- which(runif(tot)<mu)
    i <- i+1
    #(a<- which(runif(tot)<mu))
    if (length(a)>0){
      pop[a] <- pop[a]+1
    }
    
    # Back mutation
    b <- which(pop>1)
    if (length(b)>0){
      pop[b] <- 0
    }
    
    # migration
    set.seed(i)
    c <- which(runif(N) < m)
    patch[c] <- patch[c]*-1
    i <- i + 1
    # if the individual migrates, simply switch patches by changing 
    # the patch ID
    
    # phenotype and selection on phenotype
    # first 5 loci have negative effects
    # next 5 loci have positive effects
    phenotype <- rowSums(pop[,1:(nloci/2)]*-alpha + pop[,(nloci/2+1):nloci]*alpha)

    # if (t %% 1000 == 0){
    #   par(mar=c(4,4,4,1))
    #   hist(phenotype, breaks=seq((-2-alpha/2),(2+alpha/2),  by=alpha),
    #        main=paste("Gen",t))}
    # 
    if(t<burnin){
      fitness <- 1-(phenotype - 0)^2/200#rep(1, N)#
    }
    if(t >= burnin & t < (burnin+100)){
      scale <- (t-burnin)/100
      fitness <- 1-(phenotype - (scale*patch))^2/omega.sq
    }
    if(t >(burnin+100)){
      fitness <- 1-(phenotype - patch)^2/omega.sq
    }
    fitness[which(fitness<0)] <- 0
    
    # reproduction
    # sample parents according to their fitnesses
    # patch 1
    whichp1 <- which(patch==1)
    whichp2 <- which(patch==-1)
    
    set.seed(i)
    newp1 <- sample(whichp1, size = N/2, replace=TRUE, prob = fitness[whichp1])
    i <- i+1
    set.seed(i)
    newp2 <- sample(whichp2, size = N/2, replace=TRUE, prob = fitness[whichp2])
    i <- i+1
    
    newpop <- rbind(pop[newp1,], pop[newp2,])
    patch <- c(rep(1, N/2), rep(-1, N/2))
    pop <- newpop
  } # end loop through generations
  
  # calculations for final dataframe
  genotype <- apply(pop, 1, paste, collapse="")
  phenotype <- rowSums(pop[,1:(nloci/2)]*-alpha + pop[,(nloci/2+1):nloci]*alpha)
  (redund_tots <- table(genotype, phenotype))
  (segregating_redund <- colSums(redund_tots>0))
  seg_redund_df <- data.frame(phen=round(as.numeric(names(segregating_redund)),3), 
                              num_seg=segregating_redund)
  
  (genotypic_redund <- redundancy(nloci, alpha))
  
  
  final_df <- left_join(genotypic_redund, seg_redund_df, by="phen")
  final_df$num_seg[which(is.na(final_df$num_seg))] <- 0
  
  final_df$fitness_p1 <- w.zik(final_df$phen,theta = -1, omega.sq)
  final_df$fitness_p2 <- w.zik(final_df$phen,theta = 1, omega.sq)
  final_df$fitness_p1[final_df$fitness_p1<0] <- 0
  final_df$fitness_p2[final_df$fitness_p2<0] <- 0
  
  out <- hist(phenotype, breaks=c((final_df$phen-alpha/2), final_df$phen[length(final_df$phen)]+alpha/2), plot=FALSE)
  outdf <- data.frame(phen=out$mids, num_phen=out$counts)
  
  nrow(outdf)
  nrow(final_df)
  
  final_df <- cbind(final_df, outdf)
  ind_df <- data.frame(genotype, phenotype, patch)
  
  home_fitness <- w.zik(phenotype, patch, omega.sq)
  away_fitness <- w.zik(phenotype, -1*patch, omega.sq)
  
  LA <- mean(home_fitness) - mean(away_fitness)
  
  return(list(pheno_dist=final_df, geno_counts=redund_tots, ind_df=ind_df, LA=LA))
} # end function

# Sam Yeaman's pairwise C_chisquared function
pairwise_c_chisq <- function (input,num_permute = 100000,na.rm = F){
  
  numcol <- ncol (input)
  c_score <- array(NA, c(numcol,numcol))
  
  for (loop1 in 1:(numcol-1)){
    for (loop2 in (loop1+1):numcol){
      results_chisq <- array (NA,num_permute)
      input_sub <- as.matrix (input[,c(loop1,loop2)])
      obs1 <- rowSums (input_sub,na.rm = na.rm)
      exp1 <- array (mean(obs1), length(obs1))
      chisq1 <- (obs1 - exp1)^2 / exp1
      
      for (loop3 in 1:num_permute){
        input_sub[,1] <- sample (input_sub[,1],nrow (input_sub),replace = F)
        input_sub[,2] <- sample (input_sub[,2],nrow (input_sub),replace = F)
        
        obs2 <- rowSums(input_sub,na.rm = na.rm)
        exp2 <- array (mean(obs2),length(obs2))
        chisq2 <- (obs2 - exp2)^2/exp2
        
        results_chisq[loop3] <- sum (chisq2)
      }
      c_score [loop1,loop2]  <- (sum (chisq1) - mean (results_chisq)) / sd (results_chisq)
    }
  }
  mean(c_score,na.rm = T)
}

#############################################################################################################################
#Code to redund. simulate quant traits and generate C_chisquared values with full pairwise comparisons of all generated sims
#############################################################################################################################

# Sim / C score function 
# argument 're' indicates how many simulations to accumulate prior to calculating C score, default is 10
Sim_C_score<-function(N,m,Nmu,nloci,per_g,alpha,os,re=10){
  comp=1
  for(i in rep(nloci,re)){
    # Run simulation until there are 're' many replicates with a local adaptation of at least 0.1
    count=1
    nyet=FALSE
    while(nyet==FALSE){
      nsim1 <- simulatePop(N, Nmu/N, alpha, nloci=i, os, m)
      if(nsim1$LA>=0.1){
        nyet=TRUE
        nsim<<-append(nsim,list(nsim1))
        print(paste("L.A. >= 0.1, simulation ",comp,"/",re," complete",sep=""))
        comp=comp+1
      }
      else{
        count=count+1
        print(paste("L.A. < 0.1, re-running simulation ",comp,": attempt number ",count,sep="")) 
      }
    }
  }
  # create objects that will be placeholders for the output below
  sim<-NULL
  sima<-NULL
  simb<-NULL
  allelesa<-NULL
  allelesb<-NULL
  freqsa <- NULL
  freqsb <- NULL
  the_d <- NULL
  
  for(j in 1:length(nsim)){
    sim<-append(sim,list(nsim[[j]]$ind_df))
    #sim[[j]]$genotype_neut<-NA
    # loop below adds neutral loci to the genotype list from the simulation above, the number of neutral genotypes is based on the 'per_g' parameter
    #for(l in 1:nrow(sim[[j]])){
      #sim[[j]]$genotype_neut[l]<-paste(c(paste(sim[[j]]$genotype[l],collapse=""),paste(rbinom(nloci/per_g - nloci, 1, 0),collapse="")),collapse="")
    #}
    #sim[[j]]$genotype_neut<-as.factor(sim[[j]]$genotype_neut)
    sima <- append(sima,list(sim[[j]][sim[[j]][,3] > 0,]))
    simb <- append(simb,list(sim[[j]][sim[[j]][,3] < 0,]))
    gtypesa <<- append(gtypesa,list(strsplit (as.character(sima[[j]][,1]),split = "")))
    gtypesb <<- append(gtypesb,list(strsplit (as.character(simb[[j]][,1]),split = "")))
    allelesa <- append(allelesa,list(array (NA, c (nrow(sima[[j]]),length (gtypesa[[j]][[1]])))))
    allelesb <- append(allelesb,list(array (NA, c (nrow(simb[[j]]),length (gtypesb[[j]][[1]])))))
    for (k in 1:length (gtypesa[[j]][[1]])){
      allelesa[[j]][,k] <- as.numeric (sapply (gtypesa[[j]],"[[",k))
      allelesb[[j]][,k] <- as.numeric (sapply (gtypesb[[j]],"[[",k))
    }
    freqsa <- append(freqsa,list(colSums (allelesa[[j]]) / nrow (allelesa[[j]])))
    freqsb <- append(freqsb,list(colSums (allelesb[[j]]) / nrow (allelesb[[j]])))
    the_d <- append(the_d,list(abs(freqsa[[j]] - freqsb[[j]])) )
  }
  
  # Add some number of neutral "no difference between patches" loci (bunch of zeros) - number of loci is determined by 'per_g' variable 
  for(l in 1:length(the_d)){
    the_d[[l]]<-append(the_d[[l]], rep(0,nloci/per_g - nloci))
  }
  
  comb<- combinations(length(the_d),2,1:length(the_d))

  count=1
  for (i in 1:nrow(comb)){
    print(paste("Calculating C score for combination ",count,sep=""))
    the_mat <- cbind (the_d[[comb[i,1]]],the_d[[comb[i,2]]])
    C_score <<- c(C_score,pairwise_c_chisq(the_mat))
    count=count+1
  }
}

## Set parameters
N <- 1000 # Population size
Nmu <- 0.1 # Population scaled mutation rate
(mu <- Nmu/N) # Mutation rate
alpha <- 0.1 # Effect on trait
os <- 5 # omega.sq
m <- 0.20 # Mutation rate
nloci <- 50 # Number of loci under selection
per_g <- 0.5 # Percentage of total genome that loci under selection represent, if set to 1 then 100% of loci are under selection

# Receiving object creation
nsim<-list()
C_score<-NULL
gtypesa<-NULL
gtypesb<-NULL

# Run function
Sim_C_score(N,m,Nmu,nloci,per_g,alpha,os)

# Visualize C score
hist(C_score, xlab="C score", main=paste(nloci," locus case, ",per_g*100,"% of the genome",sep =""))

#nsim_20_20<-nsim
#C_score_20_20<-C_score
#nsim_50_50<-nsim
#C_score_50_50<-C_score

