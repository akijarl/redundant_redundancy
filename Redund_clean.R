## Originally: redundancy.R
## Katie E Lotterhos
## Feb 12, 2019
## Illustrate concept of genetic redundancy for different phenotypic values
## Mod. by Áki Jarl 
## Added C_chisq calculation with function from Sam Yeaman and visualization
## Aug. 28, 2019
##################################################################################
# Calculate C_chisq and show how it changes with number of loci affecting a trait
##################################################################################
packages<-c("gtools","dplyr","ggplot2","cowplot")

for(j in packages){
  if(!(j %in% installed.packages())){install.packages(j)}
}
lapply(packages,require,character.only = TRUE)

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
Sim_C_score<-function(N,m,Nmu,nloci,per_g,alpha,os,LA=0.1,re=10){
  print(paste("Running ",nloci," loci (",per_g*100,"%) simulation",sep=""))
  comp=1
  for(i in rep(nloci,re)){
    # run simulation until there are 're' many replicates with a local adaptation of at least 0.1
    cnt=1
    nyet=FALSE
    while(nyet==FALSE){
      nsim1 <- simulatePop(N, Nmu/N, alpha, nloci=i, os, m)
      if(nsim1$LA>=LA){
        nyet=TRUE
        nsim<<-append(nsim,list(nsim1))
        print(paste("L.A. >= ",LA,", simulation ",comp,"/",re," complete",sep=""))
        comp=comp+1
      }
      else{
        cnt=cnt+1
        print(paste("L.A. < ", LA,", re-running simulation ",comp,": attempt number ",cnt,sep="")) 
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
  
  # add some number of neutral "no difference between patches" loci (bunch of zeros) - number of loci is determined by 'per_g' variable 
  for(l in 1:length(the_d)){
    the_d[[l]]<-append(the_d[[l]], rep(0,nloci/per_g - nloci))
  }
  
  comb<- combinations(length(the_d),2,1:length(the_d))

  cnt=1
  for (i in 1:nrow(comb)){
    print(paste("Calculating C score for combination ",cnt,sep=""))
    the_mat <- cbind (the_d[[comb[i,1]]],the_d[[comb[i,2]]])
    C_score <<- c(C_score,pairwise_c_chisq(the_mat))
    cnt=cnt+1
  }
}

########################################
## set parameters and run simulation
########################################

N <- 1000 # population size
Nmu <- 0.1 # population scaled mutation rate
(mu <- Nmu/N) # mutation rate
alpha <- 0.1 # effect on trait
os <- 5 # omega.sq
m <- 0.20 # migration rate
genome <- 1000 # genome size
nloci <- c(20,30,40,50) # set of four different number of loci under selection
per_g <- nloci/genome # Percentage of total genome that loci under selection represent, if set to 1 then 100% of loci are under selection

# run function
for(r in 1:length(nloci)){
  nsim<-list()
  C_score<-NULL
  gtypesa<-NULL
  gtypesb<-NULL
  Sim_C_score(N,m,Nmu,nloci[r],per_g[r],alpha,os,LA=0.1,re=10)
  assign(paste("nsim_",nloci[r],sep=""),nsim)
  assign(paste("C_score_",nloci[r],sep=""),C_score)
}

# collect output for each loci scenario into one file
nsim_all<-mget(ls(pattern = "nsim_[0-9]"))
C_score_all<-mget(ls(pattern = "C_score_[0-9]"))

#nsim_all<-mget(ls(pattern = "nsim_[0-9][0-9][0-9]"))
#C_score_all<-mget(ls(pattern = "C_score_[0-9][0-9][0-9]"))

# Visualize C score
for(i in 1:length(C_score_all)){
  hist(C_score_all[[i]], xlab="C score", main=paste(nloci[i]," locus case, ",per_g[i]*100,"% of the genome",sep =""))
}

############################################
# Visualize effect of redundancy on C score
############################################

comb<- combinations(length(nsim),2,1:length(nsim))
gen_redun_final<-NULL

for(r in 1:length(nloci)){
  (genotypic_redund <- redundancy(nloci[r], alpha))
  phen<-genotypic_redund$phen
  med_phen<-NULL
  num_seg_tot<-NULL
  for(s in 1:nrow(comb)){
    med_phen<<-c(med_phen,median(abs(c(nsim_all[[r]][[comb[s,1]]]$ind_df$phenotype,nsim_all[[r]][[comb[s,2]]]$ind_df$phenotype))))
    num_seg_tot<<-c(num_seg_tot,sum(c(nsim_all[[r]][[comb[s,1]]]$pheno_dist$num_seg,nsim_all[[r]][[comb[s,2]]]$pheno_dist$num_seg)))
  }
  gen_redun<-data.frame(genotypic_redund[match(round(med_phen,10),genotypic_redund$phen),])
  row.names(gen_redun)<-1:nrow(comb)
  gen_redun$num_seg<-num_seg_tot
  gen_redun$C_score<-C_score_all[[r]]
  gen_redun_final<-rbind(gen_redun_final,gen_redun)
}

gen_redun_final$Nloci<-as.factor(c(rep(paste(per_g[1]*100,"%"),nrow(gen_redun)),rep(paste(per_g[2]*100,"%"),nrow(gen_redun)),rep(paste(per_g[3]*100,"%"),nrow(gen_redun)),rep(paste(per_g[4]*100,"%"),nrow(gen_redun))))

#can save table of simulation summary below
#write.csv(gen_redun_final, "gen_redun_final.csv",row.names = F, quote=F)

#save a plot of genotypic redundancy vs C score with its own legend
lc1<-ggplot(data=gen_redun_final,aes(x=log10Num,y=C_score))+
  geom_point(aes(col=Nloci,shape=Nloci))+
  xlab(expression(paste("Log"[10], " genotypic redundancy (at median evolved phenotypes)")))+
  ylab(expression(paste("C"[chi^2], " score")))+
  geom_smooth(method='lm',se=F)+
  labs(color='Percent loci\n affecting trait',shape='Percent loci\n affecting trait')+
  theme_classic()

#save just the legend from above
leg<-get_legend(lc1+theme(legend.margin=margin(t=5, r=0, b=0, l=0, unit="cm")))

#save plot of the number of segregating sites vs C score with no legend
sc<-ggplot(data=gen_redun_final,aes(x=num_seg,y=C_score))+
  geom_point(aes(col=Nloci,shape=Nloci))+
  xlab(expression(paste("Number of segregating sites")))+
  ylab(expression(paste("C"[chi^2], " score")))+
  geom_smooth(method='lm',se=F)+
  labs(color='Percent of genome\n under selection')+
  theme_classic()+
  guides(col=FALSE, shape=FALSE)

#save a plot of genotypic redundancy vs C score with no legend
lc1_nl<-ggplot(data=gen_redun_final,aes(x=log10Num,y=C_score))+
  geom_point(aes(col=Nloci,shape=Nloci))+
  xlab(expression(paste("Log"[10], " genotypic redundancy (at median evolved phenotypes)")))+
  ylab(expression(paste("C"[chi^2], " score")))+
  geom_smooth(method='lm',se=F)+
  labs(color='Percent loci\n affecting trait')+
  theme_classic()+
  guides(col=FALSE,shape=FALSE)

#plot combined figure of genotypic redundancy vs C score and number of segregating sites vs C score with a shared legend
plot_grid(leg,lc1_nl, NULL, sc,ncol = 2, nrow = 2,rel_widths = c(1,3,3))

#linear regression stats
summary(lm(gen_redun_final$C_score~gen_redun_final$log10Num))
cor(gen_redun_final$C_score, gen_redun_final$log10Num, method = "pearson")

summary(lm(gen_redun_final$C_score~gen_redun_final$num_seg))
cor(gen_redun_final$C_score, gen_redun_final$num_seg, method = "pearson")
#############################################################
# Visualize genotypic redund at optima
#############################################################
N <- 1000 # Population size
Nmu <- 0.1 # Population scaled mutation rate
(mu <- Nmu/N) # Mutation rate
alpha <- 0.1 # Effect on trait
os <- 5 # omega.sq
m <- 0.20 # Migration rate
nloci = 20
nsim20 <- simulatePop(N, Nmu/N, alpha, nloci, os, m, seed=2098493)
nloci = 50
nsim50 <- simulatePop(N, Nmu/N, alpha, nloci, os, m, seed=6667198)
os <- 1.45 # omega.sq
m <- 0.50 # Migration rate
nloci = 20
nsim20Hm <- simulatePop(N, Nmu/N, alpha, nloci, os, m, seed=6969281)
nloci = 50
nsim50Hm <- simulatePop(N, Nmu/N, alpha, nloci, os, m,seed=2141839)

pdf(paste0("Redund_twocases.pdf"), width=8, height= 11)

par(mar=c(2,4,2,0.5), oma=c(3,0,2,0), mfcol=c(4,2), cex=1.2)

# Left column, lower migration scenario

plot(nsim50$pheno_dist$phen, nsim50$pheno_dist$log10Num, bty="l", type="l", ylim=c(0, max(nsim50$pheno_dist$log10Num)), las=2, xlim=c(-2.5,2.5), col="blue", lty=2,
     xlab="Phenotype", ylab="Log (# genotypes)", main="A) Genotypic redundancy")
points(nsim20$pheno_dist$phen, nsim20$pheno_dist$log10Num, bty="l", type="l", ylim=c(0, max(nsim20$pheno_dist$log10Num)), las=2, xlim=c(-2.5,2.5), col="magenta")

text(0,6.5,"20 loci", col="magenta")
text(1.75,12.5,"50 loci", col="blue")

plot(nsim50$pheno_dist$phen, nsim50$pheno_dist$fitness_p1, bty="l", xlim=c(-2.5,2.5), lty=4,
     xlab="Phenotype", ylab="Individual Fitness", main="B) Fitness landscape", type="l", las=2)
points(nsim50$pheno_dist$phen, nsim50$pheno_dist$fitness_p2, type="l", lty=3)
text(-1, 0.8, "Patch\n1")
text(1, 0.8, "Patch\n2")

plot(nsim20$pheno_dist$phen, nsim20$pheno_dist$num_phen, bty="l", xlim=c(-2.5,2.5), col="magenta",
     xlab="Phenotype", ylab="Count", main="C) Evolved phenotypes", type="l", las=2)
points(nsim50$pheno_dist$phen, nsim50$pheno_dist$num_phen, bty="l", xlim=c(-2.5,2.5), col="blue", lty=2, type="l", las=2)
text(1.5,350,"20 loci", col="magenta")
text(-1.5,200,"50 loci", col="blue")

plot(nsim20$pheno_dist$phen, nsim20$pheno_dist$num_seg, bty="l", xlim=c(-2.5,2.5), ylim=c(0, max(nsim50$pheno_dist$num_seg)),
     xlab="Phenotype", ylab="Count", main="D) Segregating redundancy", type="l", las=2, col="magenta")
points(nsim50$pheno_dist$phen, nsim50$pheno_dist$num_seg, bty="l", xlim=c(-2.5,2.5), lty=2, col="blue", type="l", las=2)
text(1.5,2,"20 loci", col="magenta")
text(-1.5,4,"50 loci", col="blue")

mtext("Phenotype",side = 1, outer=TRUE, line=1, adj=0.54, cex=2.5)
mtext("Low migration",side = 3, outer=TRUE, line=0, adj=0.2, cex=2.5)
mtext("High migration",side = 3, outer=TRUE, line=0, adj=0.95, cex=2.5)

# Right column, lower migration scenario

plot(nsim20Hm$pheno_dist$phen, nsim20Hm$pheno_dist$log10Num, bty="l", type="l", ylim=c(0, max(nsim50Hm$pheno_dist$log10Num)), las=2, xlim=c(-2.5,2.5), col="magenta",
     xlab="Phenotype", ylab="", main="E) Genotypic redundancy")
points(nsim50Hm$pheno_dist$phen, nsim50Hm$pheno_dist$log10Num, bty="l", type="l", ylim=c(0, max(nsim50Hm$pheno_dist$log10Num)), las=2, xlim=c(-2.5,2.5), col="blue", lty=2, main="E) Genotypic redundancy")
text(0,6.5,"20 loci", col="magenta")
text(1.75,12.5,"50 loci", col="blue")

plot(nsim50Hm$pheno_dist$phen, nsim50Hm$pheno_dist$fitness_p1, bty="l", xlim=c(-2.5,2.5), lty=4,
     xlab="Phenotype", ylab="", main="F) Fitness landscape", type="l", las=2)
points(nsim50Hm$pheno_dist$phen, nsim50Hm$pheno_dist$fitness_p2, type="l", lty=3)
text(-1, 0.7, "Patch\n1")
text(1, 0.7, "Patch\n2")

plot(nsim20Hm$pheno_dist$phen, nsim20Hm$pheno_dist$num_phen, bty="l", xlim=c(-2.5,2.5), col="magenta",
     xlab="Phenotype", ylab="", main="G) Evolved phenotypes", type="l", las=2)
points(nsim50Hm$pheno_dist$phen, nsim50Hm$pheno_dist$num_phen, bty="l", xlim=c(-2.5,2.5), col="blue", lty=2, type="l", las=2)
text(0.6,800,"20 loci", col="magenta")
text(1.2,400,"50 loci", col="blue")

plot(nsim20Hm$pheno_dist$phen, nsim20Hm$pheno_dist$num_seg, bty="l", xlim=c(-2.5,2.5), 
     #ylim=c(0, max(nsim_50_highm$pheno_dist$num_seg)),
     xlab="Phenotype", ylab="", main="H) Segregating redundancy", type="l", las=2, col="magenta", ylim=c(0,8))
points(nsim50Hm$pheno_dist$phen, nsim50Hm$pheno_dist$num_seg, bty="l", xlim=c(-2.5,2.5), lty=2, col="blue", type="l", las=2)
text(0,6,"20 loci", col="magenta")
text(1.5,4.5,"50 loci", col="blue")

dev.off()

###########################################################################
# Visualize relationship of Gen. redundancy to Num. Seg. sites
###########################################################################
ggplot(data=gen_redun_final,aes(x=log10Num,y=num_seg))+
  geom_point(aes(col=Nloci,shape=Nloci))+
  xlab(expression(paste("Log"[10], " genotypic redundancy (at median evolved phenotypes)")))+
  ylab(expression(paste("Number of segregating sites")))+
  geom_smooth(method='lm',se=F)+
  labs(color='Percent loci\n affecting trait',shape='Percent loci\n affecting trait')+
  theme_classic()

summary(lm(gen_redun_final$num_seg~gen_redun_final$log10Num))
cor(gen_redun_final$num_seg, gen_redun_final$log10Num, method = "pearson")
