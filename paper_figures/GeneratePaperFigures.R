# figure generation


# For FDR/Sens LFC/Sens Classification Barplots
# Create Random Dataset with Ground Truth 
# analyse like normal

########################################
# G E N E R A T E   R A N D O M    D F #
########################################

library(dplyr)
library(MASS)
library(cellsigsyn)
library(parallel)

compute_ratiosR <- function(data, nreps, pseudo_count=1, lowess_norm=1, col_ctrl=NULL,
                            col_condA=NULL, col_condB=NULL, col_condAB=NULL){
  
  
  col_ctrl <- 1:nreps
  col_condA <- 1:nreps + 1*nreps 
  col_condB <- 1:nreps + 2*nreps 
  col_condAB <- 1:nreps + 3*nreps 
  
  DataPseudo <- data + pseudo_count
  reps <- length(col_ctrl)
  
  LogRatios <- matrix(nrow=length(DataPseudo[,1]), ncol=reps*5) #empty matrix
  LogIntensity <- matrix(nrow=length(DataPseudo[,1]), ncol=reps*5) #empty matrix
  
  for (i in 1:reps){
    LogRatios[,1+ (5*(i-1))] <- log2(DataPseudo[, col_condA[i]] / DataPseudo[, col_ctrl[i]])
    LogRatios[,2+ (5*(i-1))] <- log2(DataPseudo[, col_condAB[i]] / DataPseudo[, col_condB[i]])
    LogRatios[,3+ (5*(i-1))] <- log2(DataPseudo[, col_condB[i]] / DataPseudo[, col_ctrl[i]])
    LogRatios[,4+ (5*(i-1))] <- log2(DataPseudo[, col_condAB[i]] / DataPseudo[, col_condA[i]])
    LogRatios[,5+ (5*(i-1))] <- log2(DataPseudo[, col_condAB[i]] / DataPseudo[, col_ctrl[i]])
    
    LogIntensity[,1+ (5*(i-1))] <- log2(DataPseudo[, col_condA[i]] * DataPseudo[, col_ctrl[i]])
    LogIntensity[,2+ (5*(i-1))] <- log2(DataPseudo[, col_condAB[i]] * DataPseudo[, col_condB[i]])
    LogIntensity[,3+ (5*(i-1))] <- log2(DataPseudo[, col_condB[i]] * DataPseudo[, col_ctrl[i]])
    LogIntensity[,4+ (5*(i-1))] <- log2(DataPseudo[, col_condAB[i]] * DataPseudo[, col_condA[i]])
    LogIntensity[,5+ (5*(i-1))] <- log2(DataPseudo[, col_condAB[i]] * DataPseudo[, col_ctrl[i]])
  }
  n_ratios <- length(LogRatios[1,])
  
  if (lowess_norm==0) {
    message("no lowess normalization")
    rownames(LogRatios) <- rownames(data)
    return(LogRatios)}
  
  else{
    message("performing lowess normalization")
    Ratios <- matrix(ncol=n_ratios, nrow=length(LogRatios[,1]))
    pb <- txtProgressBar(min = 0, max = n_ratios, initial = 0, char = "=", style = 3)
    
    for (i in 1:n_ratios){
      Ratios[,i] <- RNAseqLowess(LogIntensity[,i], LogRatios[,i])
      setTxtProgressBar(pb,i)
    }
    rownames(Ratios) <- rownames(data)
    return(Ratios)
  }
}

parfunct <- function(ratios, nclus){
  clust <- makeCluster(nclus)
  strucvec <- rep(c(0,1,0,0,1,1,1,0,0,1,0,1,1,1,1), (length(ratios[1,])/5))
  X <- matrix(strucvec, ncol=3, byrow=T)
  clusterExport(clust, "X", envir=environment())
  mstats <- parApply(clust, ratios, 1, function(x){
    r<- summary(lm(x ~ X))$r.squared
    cov<- diag(vcov(lm(x ~ X)))
    return(c(cov,r))
  })
  stopCluster(clust)
  return(t(mstats))
}

ratios.fitR <- function(ratios, CompThreshold=1.5, n_rep=2, nclus=NULL){
  strucvec <- rep(c(0,1,0,0,1,1,1,0,0,1,0,1,1,1,1), (length(ratios[1,])/5))
  X <- matrix(strucvec, ncol=3, byrow=T)
  
  multiplefit <- lm(t(ratios)~X)
  msums <- summary(multiplefit)
  B <- t(multiplefit$coefficients)
  
  rsquared <- vector(length = length(msums))
  ngen <- length(msums)
  covB <- data.frame(cov_int=double(ngen), cov_x1=double(ngen),
                     cov_x2=double(ngen), cov_x3=double(ngen))
  
  if (is.null(nclus)){
    message("computing covariance of coefficients for random dataset...")
    pb <- txtProgressBar(min = 0, max = length(msums), initial = 0, char = "=", style = 3)
    for (i in 1:length(msums)){
      rsquared[i] <- msums[[i]]$r.squared
      covB[i,] <- diag(vcov(lm(ratios[i,]~X)))
      setTxtProgressBar(pb,i)
    }
  }
  
  if (!is.null(nclus)){
    message("computing covariance of coefficients... [parallel]")
    
    mstats <- data.frame(parfunct(ratios, nclus=nclus))
    colnames(mstats)[ncol(mstats)] <- "R2"
    rsquared <- mstats$R2
    covB <- mstats[,-ncol(mstats)]
  }
  
  n_samples <- 4*n_rep
  n_var <- 3
  DoF <- n_samples - n_var - 1
  ttest_stat <- (abs(B) - log2(CompThreshold)) / sqrt(covB) #CompThreshold is H0 hypothesis
  ttest_stat <- data.frame(ttest_stat)
  
  Pvalue = 1 - apply(ttest_stat, 2, pt, df=DoF)
  
  Datafit = B %*% t(cbind(rep(1,10), X)) #matrix multiplication to fit model
  Residual = ratios - Datafit
  Expression_variation = rowMeans(ratios**2)
  Expression_residual = rowMeans(Residual**2)
  Varexplain = 100*(Expression_variation - Expression_residual) / Expression_variation
  
  Fit <- data.frame(Expression_variation, Varexplain)
  outputdf <- data.frame(cbind(B,Pvalue,rsquared))
  colnames(outputdf)[9] <- "R2"
  return(outputdf)
}

gen_randomstats <- function(sourcedf, size=20000, nrep=2, lowessn=0, onlyDF=0,nclus=NULL){
  normed <- MedianNorm(sourcedf[rowMeans(sourcedf)>10,]) # filter out low counts
  
  message("You will be asked to enter column names for each condition.")
  message("For replicates, separate the column numbers with a comma. Example: 1,2")
  
  ctrl  <- readline("Enter column number(s) for condition ctrl:")
  condA <- readline("Enter column number(s) for condition A:")
  condB <- readline("Enter column number(s) for condition B:")
  condAB<- readline("Enter column number(s) for condition AB:")
  
  col_ctrl <- as.numeric(unlist(strsplit(ctrl, ",")))
  col_condA <- as.numeric(unlist(strsplit(condA, ",")))
  col_condB <- as.numeric(unlist(strsplit(condB, ",")))
  col_condAB <- as.numeric(unlist(strsplit(condAB, ",")))
  
  mean_per_condition <- matrix(ncol=4,nrow=nrow(normed))
  mean_per_condition[,1] <- rowMeans(normed[,col_ctrl])
  mean_per_condition[,2] <- rowMeans(normed[,col_condA])
  mean_per_condition[,3] <- rowMeans(normed[,col_condB])
  mean_per_condition[,4] <- rowMeans(normed[,col_condAB])
  
  sd_per_condition <- matrix(ncol=4,nrow=nrow(normed))
  sd_per_condition[,1] <- sqrt(rowSums((rowMeans(normed[,col_ctrl]) - normed[,col_ctrl])**2))
  sd_per_condition[,2] <- sqrt(rowSums((rowMeans(normed[,col_condA]) - normed[,col_condA])**2))
  sd_per_condition[,3] <- sqrt(rowSums((rowMeans(normed[,col_condB]) - normed[,col_condB])**2))
  sd_per_condition[,4] <- sqrt(rowSums((rowMeans(normed[,col_condAB]) - normed[,col_condAB])**2))
  
  mucov<-data.frame(mu=array(mean_per_condition), cov=array(sd_per_condition/mean_per_condition))
  lmucov <- log(mucov)
  
  #perform data binning on points variable
  olmucov <- lmucov[order(lmucov$mu),]
  binmucov <- olmucov %>% mutate(mu_bin = cut(mu, breaks=20))
  ncovs <- c()
  
  #sample cov between each bin
  for (i in unique(binmucov$mu_bin)){
    nsamp <- size
    bin_size <- nrow(binmucov[binmucov$mu_bin==i,])
    newsamplesize <- round((bin_size/nrow(binmucov))*nsamp)
    #print(newsamplesize)
    covs <- sample(binmucov[binmucov$mu_bin==i,]$cov, newsamplesize, replace=T)
    ncovs <- c(ncovs,covs)
  }
  
  mustat <-MASS::fitdistr(rowMeans(normed), "lognormal")
  new_row_means <- rlnorm(size, mustat[[1]][1], mustat[[1]][2])
  onrm <- new_row_means[order(new_row_means)]
  sampledMUCOV <- data.frame(cbind(log(onrm), ncovs))
  colnames(sampledMUCOV) <- c("mu", "cov")
  
  sampledMUCOV$mu
  sampledMUCOV$cov
  
  
  
  drawreps <- function(mucovr, rep=4*nrep){
    rnorm(n=rep, mean=exp(mucovr[1]), sd=exp(mucovr[1]+mucovr[2]))
  }
  
  gtt <- t(round(apply(sampledMUCOV,1, FUN=drawreps), 1))
  gtt[gtt<0] <- 1
  
  if(onlyDF==1){
    return(gtt)
  }
  
  #mediannorm
  normedR <- MedianNorm(gtt)
  #compute ratios
  rats <- compute_ratiosR(normedR, nreps=nrep, lowess_norm=lowessn)
  #fit model
  stats <- ratios.fitR(rats, CompThreshold=1,n_rep=nrep, nclus=nclus)
  
  # FP TN rates
  stats$SIGA <- stats[,6]<0.05 & stats$R2>0.8
  stats$SIGB <- stats[,7]<0.05 & stats$R2>0.8
  stats$SIGAB <- stats[,8]<0.05 & stats$R2>0.8
  
  FDRA <- stats[stats$SIGA,]
  FDRB <- stats[stats$SIGB,]
  FDRAB<- stats[stats$SIGAB,]
  
  FDRab <- FDRb <- FDRa <- c()
  for (i in 1:length(seq(0,4,0.001))){
    FDRab[i] <- mean(abs(FDRAB$X3)>(seq(0,4,0.001)[i]))
  }
  plot(seq(0,4,0.001), FDRab, log='x', type='l', col='green',lwd=1.5, xlim=c(0.001,2), ylim=c(0,0.99))
  for (i in 1:length(seq(0,4,0.001))){
    FDRb[i] <- mean(abs(FDRB$X2)>(seq(0,4,0.001)[i]))
  }
  lines(seq(0,4,0.001), FDRb, log='x', type='l', col='orange',lwd=1.5)
  for (i in 1:length(seq(0,4,0.001))){
    FDRa[i] <- mean(abs(FDRA$X1)>(seq(0,4,0.001)[i]))
  }
  lines(seq(0,4,0.001), FDRa, log='x', type='l', col='blue',lwd=1.5)
  abline(h=0.05, col="red", lwd=2)
  
  
  #output list of FDRs and Betas?
  
  Xs <- seq(0,4,0.001)
  FDRDF <- data.frame("Xs"=Xs, "FDRa"=FDRa,"FDRb"=FDRb, "FDRab"=FDRab)
  print("\n")
  print(nrow(FDRA))
  print(nrow(FDRB))
  print(nrow(FDRAB))
  print(paste("Threshold for A:",seq(0,4,0.001)[which.min(abs(FDRa-0.05))]))
  print(paste("Threshold for B:",seq(0,4,0.001)[which.min(abs(FDRb-0.05))]))
  print(paste("Threshold for AB:",seq(0,4,0.001)[which.min(abs(FDRab-0.05))]))
  print(seq(0,4,0.001)[which.min(abs(FDRa-0.1))])
  print(seq(0,4,0.001)[which.min(abs(FDRb-0.1))])
  print(seq(0,4,0.001)[which.min(abs(FDRab-0.1))])
  return("A B AB")
}

cts <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/IL6IL10combDF.csv", row.names=1)[,1:8]

rdf1 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf2 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf3 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf4 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf5 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf6 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf7 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf8 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf9 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf0 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)

# add ground truth counts

#scaling factor for different Betas
logFoldChangesFactor <- log2(c(1.5,1.625,1.75,1.875,2,2.125,2.25,2.375,2.5,2.625,2.75,2.875,3,3.125,3.25))
#logFoldChangesFactor2_3.25 <- log2(c(2,2,2.125,2.25,2,2.125,2.25,2.375,2.5,2.625,2.75,2.875,3,3.125,3.25))

logFoldChangesFactor2 <- c(logFoldChangesFactor, logFoldChangesFactor)

# simulate logic: A + A:B. In total 6000 ground truth genes
groundTruthScalingFactor <- matrix(rbind(matrix(rep(c(0,0,1,1,0,0,0,0),3000), ncol=8,byrow=T),
                      matrix(rep(c(0,0,-1,-1,0,0,0,0),3000), ncol=8, byrow=T)),ncol=8) * rep(logFoldChangesFactor2, 200)


emptyMatrix <- matrix(0, nrow=nrow(rdf1), ncol=8)

sampledRows <- sample(1:40000, 6000) # sample 6000 genes 
emptyMatrix[sampledRows,] <- groundTruthScalingFactor

# ground truth data set 2 signals, 2 replicates, variable ground truth LFCs
sig2rep2v0<-(2^(emptyMatrix)*rdf0)
sig2rep2v1<-(2^(emptyMatrix)*rdf1)
sig2rep2v2<-(2^(emptyMatrix)*rdf2)
sig2rep2v3<-(2^(emptyMatrix)*rdf3)
sig2rep2v4<-(2^(emptyMatrix)*rdf4)
sig2rep2v5<-(2^(emptyMatrix)*rdf5)
sig2rep2v6<-(2^(emptyMatrix)*rdf6)
sig2rep2v7<-(2^(emptyMatrix)*rdf7)
sig2rep2v8<-(2^(emptyMatrix)*rdf8)
sig2rep2v9<-(2^(emptyMatrix)*rdf9)

rownames(sig2rep2v0) <- as.character(1:nrow(sig2rep2v0))
rownames(sig2rep2v1) <- as.character(1:nrow(sig2rep2v1))
rownames(sig2rep2v2) <- as.character(1:nrow(sig2rep2v2))
rownames(sig2rep2v3) <- as.character(1:nrow(sig2rep2v3))
rownames(sig2rep2v4) <- as.character(1:nrow(sig2rep2v4))
rownames(sig2rep2v5) <- as.character(1:nrow(sig2rep2v5))
rownames(sig2rep2v6) <- as.character(1:nrow(sig2rep2v6))
rownames(sig2rep2v7) <- as.character(1:nrow(sig2rep2v7))
rownames(sig2rep2v8) <- as.character(1:nrow(sig2rep2v8))
rownames(sig2rep2v9) <- as.character(1:nrow(sig2rep2v9))

write.csv(sig2rep2v0, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf0_2sig_2rep.csv")
write.csv(sig2rep2v1, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf1_2sig_2rep.csv")
write.csv(sig2rep2v2, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf2_2sig_2rep.csv")
write.csv(sig2rep2v3, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf3_2sig_2rep.csv")
write.csv(sig2rep2v4, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf4_2sig_2rep.csv")
write.csv(sig2rep2v5, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf5_2sig_2rep.csv")
write.csv(sig2rep2v6, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf6_2sig_2rep.csv")
write.csv(sig2rep2v7, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf7_2sig_2rep.csv")
write.csv(sig2rep2v8, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf8_2sig_2rep.csv")
write.csv(sig2rep2v9, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf9_2sig_2rep.csv")

sampledRows <- as.character(sampledRows)
write(sampledRows, "sampledRows.txt") # vector with row names of ground truth genes

#all ground truth genes:
head(sig2rep2v0[sampledRows,])

all10DFs <- list(sig2rep2v0,sig2rep2v1,sig2rep2v2,sig2rep2v3,sig2rep2v4,
              sig2rep2v5,sig2rep2v6,sig2rep2v7,sig2rep2v8,sig2rep2v9)


###########################################
########## S T A T S ######################
###########################################


########################################
##### A N A L Y S E    L I N S I G #####
########################################
LinSigStats <- list()

# compute recommended thresholds -> later replace with FDR<0.05 calculation. Takes a long time though for each df...
gen_randomstats(sig2rep2v,nrep=2, size=100000, nclus=7, lowessn=0)

for (i in 1:10){
  randomDataFrame <- all10DFs[[i]]
  GTnorm <- MedianNorm(randomDataFrame)
  GTratios <- compute_ratios(GTnorm, lowess = 0, structureDataFrame = data.frame(col_ctrl=c(1,2),
                                                                                 col_condA=c(3,4),
                                                                                 col_condB=c(5,6),
                                                                                 col_condAB=c(7,8)))
  GTstat <- ratios.fit(GTratios, CompThreshold = 1)
  
  o4 <- (GTstat$cov_x3<0.05 & GTstat$R2>0.8 & abs(GTstat$X3)>0.773)#, na.rm=T)
  o4s<- (GTstat[sampledRows,]$cov_x3<0.05 & GTstat[sampledRows,]$R2>0.8 & abs(GTstat[sampledRows,]$X3)>0.773)#, na.rm=T)
  o3 <- (GTstat$cov_x2<0.05 & GTstat$R2>0.8 & abs(GTstat$X2)>0.552)#, na.rm=T)
  o3s<- (GTstat[sampledRows,]$cov_x2<0.05 & GTstat[sampledRows,]$R2>0.8 & abs(GTstat[sampledRows,]$X2)>0.552)#, na.rm=T)
  #o2 <- (GTstat$cov_x1<0.05 & GTstat$R2>0.8 & abs(GTstat$X1)>0.549)#, na.rm=T)
  #o2s<- (GTstat[sampledRows,]$cov_x1<0.05 & GTstat[sampledRows,]$R2>0.8 & abs(GTstat[sampledRows,]$X1)>0.549)#, na.rm=T)
  
  FDRv <- data.frame("FDR" = c((sum(o4, na.rm=T)-sum(o4s, na.rm=T))/sum(o4, na.rm=T),
                               (sum(o3, na.rm=T)-sum(o3s, na.rm=T))/sum(o3, na.rm=T)),
                     "term" = c("int", "main"),
                     "method" = c("own", "own"))
  SENv <- data.frame("Recall"= c(sum(o4s,na.rm=T)/6000,
                                 sum(o3s,na.rm=T)/6000),
                     "term"=c("int", "main"),
                     "method"=c("own", "own"))
  SR <- c()
  for (j in 1:30){
    SR <- append(SR, sum((GTstat[sampledRows[30*0:199+j],]$cov_x3<0.05 & GTstat[sampledRows[30*0:199+j],]$R2>0.8 & abs(GTstat[sampledRows[30*0:199+j],]$X3)>0.773), na.rm=T)/200)
  }
  for (j in 1:30){
    SR <- append(SR, sum((GTstat[sampledRows[30*0:199+j],]$cov_x2<0.05 & GTstat[sampledRows[30*0:199+j],]$R2>0.8 & abs(GTstat[sampledRows[30*0:199+j],]$X2)>0.552), na.rm=T)/200)
  }
  
  sensLFCdf <- data.frame("rate"=SR,
                          "method"=rep("own", 60),
                          "term" = rep(c("int", "main"), each=30),
                          "beta" = rep(logFoldChangesFactor2, 2),
                          "rep" = i)
  
  LinSigStats[[i]] <- list(FDRv, SENv, sensLFCdf)
}

########################################
##### A N A L Y S E     L I M M A  #####
########################################

library(limma)
library(edgeR)
cond <- c("ctrl", "LPS", "pH", "LPSpH")
stim <- factor(c("ctrl", "ctrl", "LPS", "LPS", "ctrl", "ctrl", "LPS", "LPS"), levels=c("ctrl", "LPS"))
pH <- factor(c("H","H","H", "H",  "L", "L", "L", "L"), levels=c("H", "L"))
design <- model.matrix(~stim*pH)
design
colnames(design)
cont.matrix <- cbind(HvsLinctrl=c(0,0,1,0),
                     HvsLinLPS=c(0,0,1,1),
                     Diff=c(0,0,0,1))

limmaStats <- list()
for (i in 1:10){
  randomDataFrame <- all10DFs[[i]]
  dge <- DGEList(counts=randomDataFrame)
  dge<- calcNormFactors(dge)
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  l4 <- (topTable(fit, number=Inf,coef=4)$adj.P.Val<0.05)   # over whole dataset
  l4s<- (topTable(fit, number=Inf,coef=4)[sampledRows,]$adj.P.Val<0.05) # over just ground truth
  l2 <- (topTable(fit, number=Inf,coef=2)$adj.P.Val<0.05)
  l2s<- (topTable(fit, number=Inf,coef=2)[sampledRows,]$adj.P.Val<0.05)
  
  FDRv <- data.frame("FDR" = c((sum(l4, na.rm=T)-sum(l4s, na.rm=T))/sum(l4, na.rm=T),
                               (sum(l2, na.rm=T)-sum(l2s, na.rm=T))/sum(l2, na.rm=T)),
                     "term" = c("int", "main"),
                     "method" = c("limma", "limma"))
  SENv <- data.frame("Recall"= c(sum(l4s,na.rm=T)/6000,
                                 sum(l2s,na.rm=T)/6000),
                     "term"=c("int", "main"),
                     "method"=c("limma", "limma"))
  
  SR <- c()
  for (j in 1:30){
    SR <- append(SR, (sum(topTable(fit, number=Inf,coef=4)[sampledRows[30*0:199+j],]$adj.P.Val<0.05)/200))
  }
  for (j in 1:30){
    SR <- append(SR, (sum(topTable(fit, number=Inf,coef=2)[sampledRows[30*0:199+j],]$adj.P.Val<0.05)/200))
  }
  
  sensLFCdf <- data.frame("rate"=SR,
                          "method"=rep("limma", 60),
                          "term" = rep(c("int", "main"), each=30),
                          "beta" = rep(logFoldChangesFactor2, 2),
                          "rep" = i)
  limmaStats[[i]] <- list(FDRv, SENv, sensLFCdf)
}



########################################
##### A N A L Y S E    D E S E Q 2 #####
########################################

#library(DESeq2)
metadata <- data.frame("condition"=c("Ctrl", "Ctrl", "Ctrl", "Ctrl", "Trt", "Trt", "Trt", "Trt"), 
                       "genotype"=c("WT", "WT", "MU", "MU", "WT","WT", "MU", "MU"))
rownames(metadata) <- colnames(sig2rep2v0)

deseq2Stats <- list()
for (i in 1:10){
  randomDataFrame <- all10DFs[[i]]
  
  ddsr <-  DESeqDataSetFromMatrix(round(randomDataFrame), colData=metadata, design =~ genotype + condition + genotype:condition)
  ddsr$genotype = relevel(ddsr$genotype, "WT")
  ddsr <- DESeq(ddsr)
  resrAB = results(ddsr, name="genotypeMU.conditionTrt", independentFiltering = F)
  resrAB$pvalue[is.na(resrAB$pvalue)] <- 1
  resrAB$padj[is.na(resrAB$padj)] <- 1
  sum(resrAB$pvalue<0.05, na.rm=T)
  resrB = results(ddsr, contrast=c("condition","Trt","Ctrl"), independentFiltering = F)
  resrB$pvalue[is.na(resrB$pvalue)] <- 1
  resrB$padj[is.na(resrB$padj)] <- 1
  print(paste("B", sum(resrB$pvalue<0.05, na.rm=T)))
  resrA = results(ddsr, contrast=c("genotype","MU","WT"), independentFiltering = F)
  resrA$pvalue[is.na(resrA$pvalue)] <- 1
  resrA$padj[is.na(resrA$padj)] <- 1
  print(sum(resrA$pvalue<0.05, na.rm=T))
  
  
  d4 <- (resrAB$padj<0.05)
  d4s<- (resrAB[sampledRows,]$padj<0.05)
  d2 <- (resrA$padj<0.05)
  d2s<- (resrA[sampledRows,]$padj<0.05)
  
  FDRv <- data.frame("FDR" = c((sum(d4, na.rm=T)-sum(d4s, na.rm=T))/sum(d4, na.rm=T),
                               (sum(d2, na.rm=T)-sum(d2s, na.rm=T))/sum(d2, na.rm=T)),
                     "term" = c("int", "main"),
                     "method" = c("deseq2", "deseq2"))
  SENv <- data.frame("Recall"= c(sum(d4s,na.rm=T)/6000,
                                 sum(d2s,na.rm=T)/6000),
                     "term"=c("int", "main"),
                     "method"=c("deseq2", "deseq2"))
  SR <-c()
  for (j in 1:30){
    SR <- append(SR, (sum(resrAB[sampledRows[30*0:199+j],]$padj<0.05)/200))
  }
  for (j in 1:30){
    SR <- append(SR, (sum(resrA[sampledRows[30*0:199+j],]$padj<0.05)/200))
  }
  
  sensLFCdf <- data.frame("rate"=SR,
                          "method"=rep("deseq2", 60),
                          "term" = rep(c("int", "main"), each=30),
                          "beta" = rep(logFoldChangesFactor2, 2),
                          "rep" = i)
  
  deseq2Stats[[i]] <- list(FDRv, SENv, sensLFCdf)
}


#LinSigStats
#limmaStats
#deseq2Stats

allFDRs <- rbind(LinSigStats[[1]][[1]], LinSigStats[[2]][[1]], LinSigStats[[3]][[1]],
                 LinSigStats[[4]][[1]], LinSigStats[[5]][[1]], LinSigStats[[6]][[1]],
                 LinSigStats[[7]][[1]], LinSigStats[[8]][[1]], LinSigStats[[9]][[1]],
                 LinSigStats[[10]][[1]],
                 limmaStats[[1]][[1]], limmaStats[[2]][[1]], limmaStats[[3]][[1]],
                 limmaStats[[4]][[1]], limmaStats[[5]][[1]], limmaStats[[6]][[1]],
                 limmaStats[[7]][[1]], limmaStats[[8]][[1]], limmaStats[[9]][[1]],
                 limmaStats[[10]][[1]],
                 deseq2Stats[[1]][[1]], deseq2Stats[[2]][[1]], deseq2Stats[[3]][[1]],
                 deseq2Stats[[4]][[1]], deseq2Stats[[5]][[1]], deseq2Stats[[6]][[1]],
                 deseq2Stats[[7]][[1]], deseq2Stats[[8]][[1]], deseq2Stats[[9]][[1]],
                 deseq2Stats[[10]][[1]])
allFDRs$rep <- rep(rep(1:10, each=2),3)

allSENs <- rbind(LinSigStats[[1]][[2]], LinSigStats[[2]][[2]], LinSigStats[[3]][[2]],
                 LinSigStats[[4]][[2]], LinSigStats[[5]][[2]], LinSigStats[[6]][[2]],
                 LinSigStats[[7]][[2]], LinSigStats[[8]][[2]], LinSigStats[[9]][[2]],
                 LinSigStats[[10]][[2]],
                 limmaStats[[1]][[2]], limmaStats[[2]][[2]], limmaStats[[3]][[2]],
                 limmaStats[[4]][[2]], limmaStats[[5]][[2]], limmaStats[[6]][[2]],
                 limmaStats[[7]][[2]], limmaStats[[8]][[2]], limmaStats[[9]][[2]],
                 limmaStats[[10]][[2]],
                 deseq2Stats[[1]][[2]], deseq2Stats[[2]][[2]], deseq2Stats[[3]][[2]],
                 deseq2Stats[[4]][[2]], deseq2Stats[[5]][[2]], deseq2Stats[[6]][[2]],
                 deseq2Stats[[7]][[2]], deseq2Stats[[8]][[2]], deseq2Stats[[9]][[2]],
                 deseq2Stats[[10]][[2]])
allSENs$rep <- rep(rep(1:10, each=2),3)

allFDRs$Recall <- allSENs$Recall
allfdrsens <- allFDRs

allfdrsens$col <- rep(c("red", "green", "blue"), each= 20)
allfdrsens$shape <- rep(c(18,19), 30)

write.csv(allfdrsens, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/FDR_SENS_plotdata_variableLFC1.5_4.csv")

#write.csv(allfdrsens, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/FDR_SENS_plotdata_variableLFC2_3.25.csv")


ggplot(data=allfdrsens, aes(x=Recall, y=1-FDR, fill=method)) +
  geom_point(aes(shape=term), size=2)+
  scale_shape_manual(values=c(21, 22, 24)) +
  coord_cartesian(ylim=c(0.92,0.99), xlim=c(0.7,0.99))+
  xlab("Sensitivity") +
  stat_ellipse(geom="polygon", level=0.95, aes(fill=method), alpha=0.25)


allSensVarLFC <- rbind(LinSigStats[[1]][[3]], LinSigStats[[2]][[3]], LinSigStats[[3]][[3]],
                 LinSigStats[[4]][[3]], LinSigStats[[5]][[3]], LinSigStats[[6]][[3]],
                 LinSigStats[[7]][[3]], LinSigStats[[8]][[3]], LinSigStats[[9]][[3]],
                 LinSigStats[[10]][[3]],
                 limmaStats[[1]][[3]], limmaStats[[2]][[3]], limmaStats[[3]][[3]],
                 limmaStats[[4]][[3]], limmaStats[[5]][[3]], limmaStats[[6]][[3]],
                 limmaStats[[7]][[3]], limmaStats[[8]][[3]], limmaStats[[9]][[3]],
                 limmaStats[[10]][[3]],
                 deseq2Stats[[1]][[3]], deseq2Stats[[2]][[3]], deseq2Stats[[3]][[3]],
                 deseq2Stats[[4]][[3]], deseq2Stats[[5]][[3]], deseq2Stats[[6]][[3]],
                 deseq2Stats[[7]][[3]], deseq2Stats[[8]][[3]], deseq2Stats[[9]][[3]],
                 deseq2Stats[[10]][[3]])


ggplot(data=allSensVarLFC[allSensVarLFC$beta>0,], aes(x=beta, y=rate, group=interaction(method, term), color=method, linetype=term))+
  geom_smooth(method="loess",se=T, span=0.7, aes(fill=method), alpha=0.3)+
  ylab("Sensitivity")+
  xlab("fold change") +
  coord_cartesian(ylim=c(0.4, 1))

#write.csv(allSensVarLFC, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/SENS_LFC_plotdata_variableLFC1.5_4.csv")
allSensVarLFC <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/SENS_LFC_plotdata_variableLFC1.5_4.csv")

################################################################################
################### R E G U L A T I O N    R E C O V E R Y #####################
################################################################################

#7 different types of regulations

#dataframe with all types of regulation.
#For each regulation compute the number of wrong logic/no logic/right logic?

f=log2(2)
A <-  c(0,0,f,f,0,0,f,f)
B <-  c(0,0,0,0,f,f,f,f)
AB<-  c(0,0,0,0,0,0,f,f)
BA<-  c(0,0,f,f,-f,-f,0,0)
AAB<- c(0,0,f,f,0,0,0,0)
BAB<- c(0,0,0,0,f,f,0,0)
ABAB<-c(0,0,f,f,f,f,f,f)

cts <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/IL6IL10combDF.csv", row.names=1)[,1:8]

rdf1 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf2 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf3 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf4 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf5 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf6 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf7 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf8 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf9 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)
rdf0 <- gen_randomstats(cts,nrep=2, size=40000, onlyDF=1)

GTmatrix <- matrix(rep(c(A,B,AB,BA,AAB,BAB,ABAB),750), ncol=8,byrow = T)
sampledRows <- sample(2000:40000, 5250)

emptyMatrix <- matrix(0, nrow=nrow(rdf1), ncol=8)
emptyMatrix[sampledRows,] <- GTmatrix * sample(c(-1,1), 5250, replace=TRUE)

# ground truth data set 2 signals, 2 replicates, variable ground truth LFCs
sig2rep2RR0<-(2^(emptyMatrix)*rdf0) # RR for Regulation Recovery
sig2rep2RR1<-(2^(emptyMatrix)*rdf1)
sig2rep2RR2<-(2^(emptyMatrix)*rdf2)
sig2rep2RR3<-(2^(emptyMatrix)*rdf3)
sig2rep2RR4<-(2^(emptyMatrix)*rdf4)
sig2rep2RR5<-(2^(emptyMatrix)*rdf5)
sig2rep2RR6<-(2^(emptyMatrix)*rdf6)
sig2rep2RR7<-(2^(emptyMatrix)*rdf7)
sig2rep2RR8<-(2^(emptyMatrix)*rdf8)
sig2rep2RR9<-(2^(emptyMatrix)*rdf9)

sig2rep2RR0[sig2rep2RR0<0] <-0
sig2rep2RR1[sig2rep2RR1<0] <-0
sig2rep2RR2[sig2rep2RR2<0] <-0
sig2rep2RR3[sig2rep2RR3<0] <-0
sig2rep2RR4[sig2rep2RR4<0] <-0
sig2rep2RR5[sig2rep2RR5<0] <-0
sig2rep2RR6[sig2rep2RR6<0] <-0
sig2rep2RR7[sig2rep2RR7<0] <-0
sig2rep2RR8[sig2rep2RR8<0] <-0
sig2rep2RR9[sig2rep2RR9<0] <-0

rownames(sig2rep2RR0) <- as.character(1:nrow(sig2rep2RR0))
rownames(sig2rep2RR1) <- as.character(1:nrow(sig2rep2RR1))
rownames(sig2rep2RR2) <- as.character(1:nrow(sig2rep2RR2))
rownames(sig2rep2RR3) <- as.character(1:nrow(sig2rep2RR3))
rownames(sig2rep2RR4) <- as.character(1:nrow(sig2rep2RR4))
rownames(sig2rep2RR5) <- as.character(1:nrow(sig2rep2RR5))
rownames(sig2rep2RR6) <- as.character(1:nrow(sig2rep2RR6))
rownames(sig2rep2RR7) <- as.character(1:nrow(sig2rep2RR7))
rownames(sig2rep2RR8) <- as.character(1:nrow(sig2rep2RR8))
rownames(sig2rep2RR9) <- as.character(1:nrow(sig2rep2RR9))

write.csv(sig2rep2RR0, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf0_2sig_2rep_RR.csv")
write.csv(sig2rep2RR1, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf1_2sig_2rep_RR.csv")
write.csv(sig2rep2RR2, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf2_2sig_2rep_RR.csv")
write.csv(sig2rep2RR3, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf3_2sig_2rep_RR.csv")
write.csv(sig2rep2RR4, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf4_2sig_2rep_RR.csv")
write.csv(sig2rep2RR5, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf5_2sig_2rep_RR.csv")
write.csv(sig2rep2RR6, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf6_2sig_2rep_RR.csv")
write.csv(sig2rep2RR7, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf7_2sig_2rep_RR.csv")
write.csv(sig2rep2RR8, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf8_2sig_2rep_RR.csv")
write.csv(sig2rep2RR9, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf9_2sig_2rep_RR.csv")

sampledRows <- as.character(sampledRows)
write(sampledRows, "sampledRows_RR.txt") # vector with row names of ground truth genes

sig2rep2RR0 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf0_2sig_2rep_RR.csv", row.names = 1)
sig2rep2RR1 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf1_2sig_2rep_RR.csv", row.names = 1)
sig2rep2RR2 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf2_2sig_2rep_RR.csv", row.names = 1)
sig2rep2RR3 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf3_2sig_2rep_RR.csv", row.names = 1)
sig2rep2RR4 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf4_2sig_2rep_RR.csv", row.names = 1)
sig2rep2RR5 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf5_2sig_2rep_RR.csv", row.names = 1)
sig2rep2RR6 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf6_2sig_2rep_RR.csv", row.names = 1)
sig2rep2RR7 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf7_2sig_2rep_RR.csv", row.names = 1)
sig2rep2RR8 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf8_2sig_2rep_RR.csv", row.names = 1)
sig2rep2RR9 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/rdf9_2sig_2rep_RR.csv", row.names = 1)

sampledRows <- readLines("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/sampledRows_RR.txt")

##### LIMMA #####
library(limma)
library(edgeR)
library(DESeq2)

cond <- c("ctrl", "LPS", "pH", "LPSpH")
stim <- factor(c("ctrl", "ctrl", "LPS", "LPS", "ctrl", "ctrl", "LPS", "LPS"), levels=c("ctrl", "LPS"))
pH <- factor(c("H","H","H", "H",  "L", "L", "L", "L"), levels=c("H", "L"))
design <- model.matrix(~stim*pH)
design
colnames(design)
cont.matrix <- cbind(HvsLinctrl=c(0,0,1,0),
                     HvsLinLPS=c(0,0,1,1),
                     Diff=c(0,0,0,1))

cTrueRegDF <- c()
cMisTrueDF <- c()
cregFDRDF <- c()
cTMFdf <- c()

dfs <- list(sig2rep2RR0, sig2rep2RR1, sig2rep2RR2, sig2rep2RR3, sig2rep2RR4, 
            sig2rep2RR5, sig2rep2RR6, sig2rep2RR7, sig2rep2RR8, sig2rep2RR9)


for (i in dfs){
  print(head(i))
  vregStats <- fit_intmodel(i, size=40000, nclus=7, structureDataFrame = data.frame(col_ctrl=c(1,2),
                                                                                    col_condA=c(3,4),
                                                                                    col_condB=c(5,6),
                                                                                    col_condAB=c(7,8)))
  
  onlyRegulatedGenes <- data.frame(vregStats[as.character(sampledRows),])
  
  otruereg <- c(sum(onlyRegulatedGenes$FDRA<0.05 & onlyRegulatedGenes$FDRB>0.05 & onlyRegulatedGenes$FDRAB>0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1<0.05 & onlyRegulatedGenes$cov_x2>0.00 & onlyRegulatedGenes$cov_x3>0.00 & rep(1:7, 750)==2, na.rm=T), #
                #just B rep(1:7, 750)==2
                sum(onlyRegulatedGenes$FDRA>0.05 & onlyRegulatedGenes$FDRB<0.05 & onlyRegulatedGenes$FDRAB>0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1>0.00 & onlyRegulatedGenes$cov_x2<0.05 & onlyRegulatedGenes$cov_x3>0.00 & rep(1:7, 750)==1, na.rm=T),
                #just A+B rep(1:7, 750)==3
                sum(onlyRegulatedGenes$FDRA<0.05 & onlyRegulatedGenes$FDRB<0.05 & onlyRegulatedGenes$FDRAB>0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1<0.05 & onlyRegulatedGenes$cov_x2<0.05 & onlyRegulatedGenes$cov_x3>0 & rep(1:7, 750)==4, na.rm=T),
                #just AB rep(1:7, 750)==4
                sum(onlyRegulatedGenes$FDRA>0.05 & onlyRegulatedGenes$FDRB>0.05 & onlyRegulatedGenes$FDRAB<0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1>0.00 & onlyRegulatedGenes$cov_x2>0.00 & onlyRegulatedGenes$cov_x3<0.05 & rep(1:7, 750)==3, na.rm=T),
                #just A+AB rep(1:7, 750)==5
                sum(onlyRegulatedGenes$FDRA<0.05 & onlyRegulatedGenes$FDRB>0.05 & onlyRegulatedGenes$FDRAB<0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1<0.05 & onlyRegulatedGenes$cov_x2>0.00 & onlyRegulatedGenes$cov_x3<0.05 & rep(1:7, 750)==6, na.rm=T),
                #just B+AB rep(1:7, 750)==6
                sum(onlyRegulatedGenes$FDRA>0.05 & onlyRegulatedGenes$FDRB<0.05 & onlyRegulatedGenes$FDRAB<0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1>0.00 & onlyRegulatedGenes$cov_x2<0.05 & onlyRegulatedGenes$cov_x3<0.05 & rep(1:7, 750)==5, na.rm=T),
                #just A+B+AB rep(1:7, 750)==7
                sum(onlyRegulatedGenes$FDRA<0.05 & onlyRegulatedGenes$FDRB<0.05 & onlyRegulatedGenes$FDRAB<0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1<0.05 & onlyRegulatedGenes$cov_x2<0.05 & onlyRegulatedGenes$cov_x3<0.05 & rep(1:7, 750)==7, na.rm=T))
  
  #regrec
  omistrue <- c(sum(onlyRegulatedGenes$FDRA<0.05 & onlyRegulatedGenes$FDRB>0.05 & onlyRegulatedGenes$FDRAB>0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1<0.05 & onlyRegulatedGenes$cov_x2>0.00 & onlyRegulatedGenes$cov_x3>0.00 ,na.rm=T), #
                #just B rep(1:7, 750)==2
                sum(onlyRegulatedGenes$FDRA>0.05 & onlyRegulatedGenes$FDRB<0.05 & onlyRegulatedGenes$FDRAB>0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1>0.00 & onlyRegulatedGenes$cov_x2<0.05 & onlyRegulatedGenes$cov_x3>0.00, na.rm=T),
                #just A+B rep(1:7, 750)==3
                sum(onlyRegulatedGenes$FDRA<0.05 & onlyRegulatedGenes$FDRB<0.05 & onlyRegulatedGenes$FDRAB>0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1<0.05 & onlyRegulatedGenes$cov_x2<0.05 & onlyRegulatedGenes$cov_x3>0.00, na.rm=T),
                #just AB rep(1:7, 750)==4
                sum(onlyRegulatedGenes$FDRA>0.05 & onlyRegulatedGenes$FDRB>0.05 & onlyRegulatedGenes$FDRAB<0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1>0.00 & onlyRegulatedGenes$cov_x2>0.00 & onlyRegulatedGenes$cov_x3<0.05, na.rm=T),
                #just A+AB rep(1:7, 750)==5
                sum(onlyRegulatedGenes$FDRA<0.05 & onlyRegulatedGenes$FDRB>0.05 & onlyRegulatedGenes$FDRAB<0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1<0.05 & onlyRegulatedGenes$cov_x2>0.00 & onlyRegulatedGenes$cov_x3<0.05, na.rm=T),
                #just B+AB rep(1:7, 750)==6
                sum(onlyRegulatedGenes$FDRA>0.05 & onlyRegulatedGenes$FDRB<0.05 & onlyRegulatedGenes$FDRAB<0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1>0.00 & onlyRegulatedGenes$cov_x2<0.05 & onlyRegulatedGenes$cov_x3<0.05, na.rm=T),
                #just A+B+AB rep(1:7, 750)==7
                sum(onlyRegulatedGenes$FDRA<0.05 & onlyRegulatedGenes$FDRB<0.05 & onlyRegulatedGenes$FDRAB<0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1<0.05 & onlyRegulatedGenes$cov_x2<0.05 & onlyRegulatedGenes$cov_x3<0.05, na.rm=T))
  #Tregrec
  oallreg <- c(sum(vregStats$FDRA<0.05 & vregStats$FDRB>0.05 & vregStats$FDRAB>0.05 & vregStats$R2>0.8 &
                     vregStats$cov_x1<0.05 & vregStats$cov_x2>0.00 & vregStats$cov_x3>0.00 ,na.rm=T), #
               #just B rep(1:7, 750)==2
               sum(vregStats$FDRA>0.05 & vregStats$FDRB<0.05 & vregStats$FDRAB>0.05 & vregStats$R2>0.8 &
                     vregStats$cov_x1>0.00 & vregStats$cov_x2<0.05 & vregStats$cov_x3>0.00, na.rm=T),
               #just A+B rep(1:7, 750)==3
               sum(vregStats$FDRA<0.05 & vregStats$FDRB<0.05 & vregStats$FDRAB>0.05 & vregStats$R2>0.8 &
                     vregStats$cov_x1<0.05 & vregStats$cov_x2<0.05 & vregStats$cov_x3>0.00, na.rm=T),
               #just AB rep(1:7, 750)==4
               sum(vregStats$FDRA>0.05 & vregStats$FDRB>0.05 & vregStats$FDRAB<0.05 & vregStats$R2>0.8 &
                     vregStats$cov_x1>0.00 & vregStats$cov_x2>0.00 & vregStats$cov_x3<0.05, na.rm=T),
               #just A+AB rep(1:7, 750)==5
               sum(vregStats$FDRA<0.05 & vregStats$FDRB>0.05 & vregStats$FDRAB<0.05 & vregStats$R2>0.8 &
                     vregStats$cov_x1<0.05 & vregStats$cov_x2>0.00 & vregStats$cov_x3<0.05, na.rm=T),
               #just B+AB rep(1:7, 750)==6
               sum(vregStats$FDRA>0.05 & vregStats$FDRB<0.05 & vregStats$FDRAB<0.05 & vregStats$R2>0.8 &
                     vregStats$cov_x1>0.00 & vregStats$cov_x2<0.05 & vregStats$cov_x3<0.05, na.rm=T),
               #just A+B+AB rep(1:7, 750)==7
               sum(vregStats$FDRA<0.05 & vregStats$FDRB<0.05 & vregStats$FDRAB<0.05 & vregStats$R2>0.8 &
                     vregStats$cov_x1<0.05 & vregStats$cov_x2<0.05 & vregStats$cov_x3<0.05, na.rm=T))
  
  regs<-c("A","B", "A+B", "A:B", "A+AB", "B+AB", "A+B+AB")
  
  # compare with other methods
  dge <- DGEList(counts=i)
  dge<- calcNormFactors(dge)
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  #sampledRows <- as.character(sampledRows)
  
  l4 <- (topTable(fit, number=Inf,coef=4, sort="none"))
  l4s<- (topTable(fit, number=Inf,coef=4)[as.character(sampledRows),])
  l3 <- (topTable(fit, number=Inf,coef=3, sort="none"))
  l3s<- (topTable(fit, number=Inf,coef=3)[sampledRows,])
  l2 <- (topTable(fit, number=Inf,coef=2, sort="none"))
  l2s<- (topTable(fit, number=Inf,coef=2)[sampledRows,])
  
  litruereg<- c(sum(l2s$adj.P.Val<0.05 & l3s$adj.P.Val>0.05 & l4s$adj.P.Val>0.05 & rep(1:7, 750)==1),
                sum(l2s$adj.P.Val>0.05 & l3s$adj.P.Val<0.05 & l4s$adj.P.Val>0.05 & rep(1:7, 750)==2),
                sum(l2s$adj.P.Val<0.05 & l3s$adj.P.Val<0.05 & l4s$adj.P.Val>0.05 & rep(1:7, 750)==4),
                sum(l2s$adj.P.Val>0.05 & l3s$adj.P.Val>0.05 & l4s$adj.P.Val<0.05 & rep(1:7, 750)==3),
                sum(l2s$adj.P.Val<0.05 & l3s$adj.P.Val>0.05 & l4s$adj.P.Val<0.05 & rep(1:7, 750)==5),
                sum(l2s$adj.P.Val>0.05 & l3s$adj.P.Val<0.05 & l4s$adj.P.Val<0.05 & rep(1:7, 750)==6),
                sum(l2s$adj.P.Val<0.05 & l3s$adj.P.Val<0.05 & l4s$adj.P.Val<0.05 & rep(1:7, 750)==7))
  
  #Aliregsens
  limistrue<- c(sum(l2s$adj.P.Val<0.05 & l3s$adj.P.Val>0.05 & l4s$adj.P.Val>0.05),
                sum(l2s$adj.P.Val>0.05 & l3s$adj.P.Val<0.05 & l4s$adj.P.Val>0.05),
                sum(l2s$adj.P.Val<0.05 & l3s$adj.P.Val<0.05 & l4s$adj.P.Val>0.05),
                sum(l2s$adj.P.Val>0.05 & l3s$adj.P.Val>0.05 & l4s$adj.P.Val<0.05),
                sum(l2s$adj.P.Val<0.05 & l3s$adj.P.Val>0.05 & l4s$adj.P.Val<0.05),
                sum(l2s$adj.P.Val>0.05 & l3s$adj.P.Val<0.05 & l4s$adj.P.Val<0.05),
                sum(l2s$adj.P.Val<0.05 & l3s$adj.P.Val<0.05 & l4s$adj.P.Val<0.05))
  
  #TAliregsens
  liallreg<- c(sum(l2$adj.P.Val<0.05 & l3$adj.P.Val>0.05 & l4$adj.P.Val>0.05),
               sum(l2$adj.P.Val>0.05 & l3$adj.P.Val<0.05 & l4$adj.P.Val>0.05),
               sum(l2$adj.P.Val<0.05 & l3$adj.P.Val<0.05 & l4$adj.P.Val>0.05),
               sum(l2$adj.P.Val>0.05 & l3$adj.P.Val>0.05 & l4$adj.P.Val<0.05),
               sum(l2$adj.P.Val<0.05 & l3$adj.P.Val>0.05 & l4$adj.P.Val<0.05),
               sum(l2$adj.P.Val>0.05 & l3$adj.P.Val<0.05 & l4$adj.P.Val<0.05),
               sum(l2$adj.P.Val<0.05 & l3$adj.P.Val<0.05 & l4$adj.P.Val<0.05))
  
  
  ##### DESeq2 #####
  
  metadata <- data.frame("condition"=c("Ctrl", "Ctrl", "Ctrl", "Ctrl", "Trt", "Trt", "Trt", "Trt"), "genotype"=c("WT", "WT", "MU", "MU", "WT","WT", "MU", "MU"))
  rownames(metadata) <- colnames(onlyRegulatedGenes[1:8])
  
  rownames(metadata) <- c("V1","V2","V3","V4","V5","V6","V7","V8")
  
  ddsr <-  DESeqDataSetFromMatrix(round(i), colData=metadata, design =~ genotype + condition + genotype:condition)
  
  ddsr$genotype = relevel(ddsr$genotype, "WT")
  ddsr <- DESeq(ddsr)
  resrAB = results(ddsr, name="genotypeMU.conditionTrt", independentFiltering = F)
  resrAB$pvalue[is.na(resrAB$pvalue)] <- 1
  resrAB$padj[is.na(resrAB$padj)] <- 1
  sum(resrAB$pvalue<0.05, na.rm=T)
  resrB = results(ddsr, contrast=c("condition","Trt","Ctrl"), independentFiltering = F)
  resrB$pvalue[is.na(resrB$pvalue)] <- 1
  resrB$padj[is.na(resrB$padj)] <- 1
  sum(resrB$pvalue<0.05, na.rm=T)
  resrA = results(ddsr, contrast=c("genotype","MU","WT"), independentFiltering = F)
  resrA$pvalue[is.na(resrA$pvalue)] <- 1
  resrA$padj[is.na(resrA$padj)] <- 1
  sum(resrA$pvalue<0.05, na.rm=T)
  
  d2truereg <-c(sum(resrA[sampledRows,]$padj<0.05 & resrB[sampledRows,]$padj>0.05 & resrAB[sampledRows,]$padj>0.05 & rep(1:7, 750)==1),
                sum(resrA[sampledRows,]$padj>0.05 & resrB[sampledRows,]$padj<0.05 & resrAB[sampledRows,]$padj>0.05 & rep(1:7, 750)==2),
                sum(resrA[sampledRows,]$padj<0.05 & resrB[sampledRows,]$padj<0.05 & resrAB[sampledRows,]$padj>0.05 & rep(1:7, 750)==4),
                sum(resrA[sampledRows,]$padj>0.05 & resrB[sampledRows,]$padj>0.05 & resrAB[sampledRows,]$padj<0.05 & rep(1:7, 750)==3),
                sum(resrA[sampledRows,]$padj<0.05 & resrB[sampledRows,]$padj>0.05 & resrAB[sampledRows,]$padj<0.05 & rep(1:7, 750)==5),
                sum(resrA[sampledRows,]$padj>0.05 & resrB[sampledRows,]$padj<0.05 & resrAB[sampledRows,]$padj<0.05 & rep(1:7, 750)==6),
                sum(resrA[sampledRows,]$padj<0.05 & resrB[sampledRows,]$padj<0.05 & resrAB[sampledRows,]$padj<0.05 & rep(1:7, 750)==7))
  
  #Ad2regsens
  d2mistrue <-c(sum(resrA[sampledRows,]$padj<0.05 & resrB[sampledRows,]$padj>0.05 & resrAB[sampledRows,]$padj>0.05),
                sum(resrA[sampledRows,]$padj>0.05 & resrB[sampledRows,]$padj<0.05 & resrAB[sampledRows,]$padj>0.05),
                sum(resrA[sampledRows,]$padj<0.05 & resrB[sampledRows,]$padj<0.05 & resrAB[sampledRows,]$padj>0.05),
                sum(resrA[sampledRows,]$padj>0.05 & resrB[sampledRows,]$padj>0.05 & resrAB[sampledRows,]$padj<0.05),
                sum(resrA[sampledRows,]$padj<0.05 & resrB[sampledRows,]$padj>0.05 & resrAB[sampledRows,]$padj<0.05),
                sum(resrA[sampledRows,]$padj>0.05 & resrB[sampledRows,]$padj<0.05 & resrAB[sampledRows,]$padj<0.05),
                sum(resrA[sampledRows,]$padj<0.05 & resrB[sampledRows,]$padj<0.05 & resrAB[sampledRows,]$padj<0.05))
  #TAd2regsens
  d2allreg <-c(sum(resrA$padj<0.05 & resrB$padj>0.05 & resrAB$padj>0.05),
               sum(resrA$padj>0.05 & resrB$padj<0.05 & resrAB$padj>0.05),
               sum(resrA$padj<0.05 & resrB$padj<0.05 & resrAB$padj>0.05),
               sum(resrA$padj>0.05 & resrB$padj>0.05 & resrAB$padj<0.05),
               sum(resrA$padj<0.05 & resrB$padj>0.05 & resrAB$padj<0.05),
               sum(resrA$padj>0.05 & resrB$padj<0.05 & resrAB$padj<0.05),
               sum(resrA$padj<0.05 & resrB$padj<0.05 & resrAB$padj<0.05))
  
  
  TrueRegDF <- data.frame("Recall"=c(otruereg, d2truereg, litruereg),
                          "RegLgc"=rep(factor(c("A", "B", "A+B", "A:B", "A+AB", "B+AB", "A+B+AB"), level=c("A", "B", "A+B", "A:B", "A+AB", "B+AB", "A+B+AB")),3),
                          "Method"= rep(c("own", "deseq2", "limma"),each=7))
  cTrueRegDF <- rbind(cTrueRegDF, TrueRegDF)
  
  MisTrueDF <- data.frame("Recall"=c(omistrue, d2mistrue, limistrue),
                          "RegLgc"=rep(factor(c("A", "B", "A+B", "A:B", "A+AB", "B+AB", "A+B+AB"), level=c("A", "B", "A+B", "A:B", "A+AB", "B+AB", "A+B+AB")),3),
                          "Method"= rep(c("own", "deseq2", "limma"),each=7))
  cMisTrueDF <- rbind(cMisTrueDF, MisTrueDF)
  
  regFDRDF <- data.frame("FDR"=c((oallreg-omistrue), (d2allreg-d2mistrue), (liallreg-limistrue)),
                         "RegLgc"=rep(factor(c("A", "B", "A+B", "A:B", "A+AB", "B+AB", "A+B+AB"), level=c("A", "B", "A+B", "A:B", "A+AB", "B+AB", "A+B+AB")),3),
                         "Method"= rep(c("own", "deseq2", "limma"),each=7))
  cregFDRDF <- rbind(cregFDRDF, regFDRDF)
  
  
  TMFdf <- data.frame("rate"=c(TrueRegDF$Recall, (MisTrueDF$Recall-TrueRegDF$Recall), regFDRDF$FDR),
                      "discovery"=c(rep(c("True Disc", "Mis Class", "False Disc"),each=21)),
                      "method"=rep(TrueRegDF$Method, 3),
                      "Regulation"=rep(TrueRegDF$RegLgc, 3),
                      "replicate"=rep(10, 63)
  )
  print(TMFdf)
  cTMFdf <- rbind(cTMFdf, TMFdf)
}


write.csv(cTMFdf, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/2sig_2rep_RegulationRecoveryData.csv")
cTMFdf <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/2sig_2rep_RegulationRecoveryData.csv")

RegulationGroups <- c("A/B", "A/B", "A+B","A:B", "A/B+A:B","A/B+A:B", "A/B+A:B")
cTMFdf$RegulationGroups <- rep(RegulationGroups, 90)



statTMFdf <- cTMFdf %>%
  group_by(discovery, method, RegulationGroups) %>%
  summarise(
    Weight = mean(rate),
    sd = sd(rate),
    n = n(),
    se = sd / sqrt(n)
  )

ggplot(data=statTMFdf, aes(x = method, y = Weight, fill = discovery)) +
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ RegulationGroups) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_hline(yintercept=750, linetype="dashed", color="red")

