# figure generation for 2 replicates


# For FDR/Sens LFC/Sens Classification Barplots
# Create Random Dataset with Ground Truth 
# analyse like normal

########################################
# G E N E R A T E   R A N D O M    D F #
########################################

library(dplyr)
library(MASS)
#library(cellsigsyn)
library(parallel)
library(ggplot2)

library(limma)
library(edgeR)
library(DESeq2)

source("~/Boston Internship/Github/Rsyn/paper_figures/basic_model_functions/data_standardization.R") # for MedianNorm function
source("~/Boston Internship/Github/Rsyn/paper_figures/gen_synthetic_data_helpers.R")
source("~/Boston Internship/Github/Rsyn/paper_figures/basic_model_functions/fit_model.R") # for RNAseqLowess function
source("~/Boston Internship/Github/Rsyn/paper_figures/compute_lfc_thresholds.R")
source("~/Boston Internship/Github/Rsyn/paper_figures/basic_model_functions/compute_ratios.R")

cts <- read.csv("~/Boston Internship/Github/Rsyn/paper_figures/IL6IL10combDF.csv", row.names=1)[,1:8]
cts <- read.csv("~/Boston Internship/Subdata_cts_q.csv", row.names=1)[,1:8]
dim(cts)
#### Functions ####

generateGroundTruthGenes <- function(
    n_ground_truth_genes, n_total_genes, 
    range, fold_change, regulation, nrep=2, randomize_pos_neg=T){
  
  sampled_rows <- sample(range, n_ground_truth_genes)
  log_fold_changes <- rep(fold_change, n_ground_truth_genes/length(fold_change)) # length(fold_change) should be a proper divisor of n_ground_truth_genes
  full_matrix <- matrix(0, nrow=n_total_genes, ncol=nrep*4)
  gt_matrix <- matrix(
    rep(regulation, 
        n_ground_truth_genes/(length(regulation)/(nrep*4)))
    ,ncol=nrep*4, byrow=T)
  
  if (randomize_pos_neg){ random_positive_negative_multiplier <- sample(c(-1,1), n_ground_truth_genes, replace=TRUE) }
  else { random_positive_negative_multiplier <- 1 }
  
  full_matrix[sampled_rows,] <- gt_matrix * log_fold_changes * random_positive_negative_multiplier
  
  ground_truth_data <- data.frame(
    sampled_rows = as.character(sampled_rows),
    lfc = random_positive_negative_multiplier * log_fold_changes,
    regulation = rep(1:(length(regulation)/(nrep*4)), 
                     each = n_ground_truth_genes/(length(regulation)/(nrep*4)))
  )
  
  return(list(emptyMatrix=full_matrix, ground_truth_data=ground_truth_data))
}

### PARAMETERS ###
params <- new.env(parent = emptyenv())
params$struc_vec <- c(0,1,0, # Only A
                      0,1,1, # A + A:B
                      1,0,0, # Only B
                      1,0,1, # B + A:B
                      1,1,1) # A + B + A:B
params$strucDF <- data.frame(
  col_ctrl = c(1, 2),      # Columns for CTRL (e.g., replicates 1 and 2)
  col_condA = c(3, 4),     # Columns for Condition A (e.g., replicates 3 and 4)
  col_condB = c(5, 6),     # Columns for Condition B (e.g., replicates 5 and 6)
  col_condAB = c(7, 8)     # Columns for Condition AB (e.g., replicates 7 and 8)
)
params$n_synth_dfs <- 10
params$n_genes_simulated <- 40000
params$lfc_reg_recovery_figure <- log2(2)
params$lfc_fdr_sensitivity_per_lfc_figure <- log2(seq(1.5, 3.25, by = 0.125))
params$x_axis_thresholds <- seq(0.001, 4.001, 0.001) # initialize an X-axis for the LFC threshold figure (LFC 0.001-4.001, with a small step size for high resolution curve)



###################################################################
datasets <- generate_multiple_datasets(source_df = cts, nrep = 2, size = params$n_genes_simulated)


gt_df <- generateGroundTruthGenes(n_ground_truth_genes = 6000,
                                 n_total_genes = 40000,
                                 range=1:40000,
                                 fold_change = params$lfc_fdr_sensitivity_per_lfc_figure,
                                 regulation = c(0, 0,  1,  1, 0, 0, 0, 0),
                                 randomize_pos_neg = T)


ground_truth_datasets <- lapply(datasets, function(df) {
  scaled_df <- 2^gt_df$emptyMatrix * df
  rownames(scaled_df) <- as.character(1:nrow(scaled_df))
  return(scaled_df)
})

sampledRows <- gt_df$ground_truth_data$sampled_rows


###########################################
########## S T A T S ######################
###########################################

##### A N A L Y S E    L I N S I G #####
LinSigStats <- list()

# compute recommended thresholds -> later replace with FDR<0.05 calculation.
compute_lfc_thresholds(ground_truth_datasets[[2]],nrep=2, size=100000, lowessn=0)

for (i in 1:params$n_synth_dfs){
  randomDataFrame <- ground_truth_datasets[[i]]
  GTnorm <- MedianNorm(randomDataFrame)
  GTratios <- compute_ratios(GTnorm, lowess = 0, structureDataFrame = params$strucDF)
  GTstat <- ratios.fit(GTratios, CompThreshold = 1)
  
  lfc_thresholds <- compute_lfc_thresholds(randomDataFrame,nrep=2, size=100000, lowessn=0)

  interaction_threshold <- lfc_thresholds["AB",]$Threshold_at_5pct
  main_threshold <- lfc_thresholds["B",]$Threshold_at_5pct
  
  o4 <- (GTstat$cov_x3<0.05 & GTstat$R2>0.8 & abs(GTstat$X3) > interaction_threshold)
  o4s<- (GTstat[sampledRows,]$cov_x3<0.05 & GTstat[sampledRows,]$R2>0.8 & abs(GTstat[sampledRows,]$X3)>interaction_threshold)
  o3 <- (GTstat$cov_x2<0.05 & GTstat$R2>0.8 & abs(GTstat$X2) >  main_threshold)
  o3s<- (GTstat[sampledRows,]$cov_x2<0.05 & GTstat[sampledRows,]$R2>0.8 & abs(GTstat[sampledRows,]$X2)>main_threshold)

  
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
                          "beta" = rep(params$lfc_fdr_sensitivity_per_lfc_figure, 2),
                          "rep" = i)
  
  LinSigStats[[i]] <- list(fdr_df=FDRv, sen_df=SENv, sen_lfc_df=sensLFCdf)
}

##### A N A L Y S E     L I M M A  #####
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
for (i in 1:params$n_synth_dfs){
  randomDataFrame <- ground_truth_datasets[[i]]
  dge <- DGEList(counts=randomDataFrame)
  dge<- calcNormFactors(dge)
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  int_effect <- topTable(fit, number=Inf,coef=4)
  main_effect <- topTable(fit, number=Inf,coef=2)
  
  l4 <- int_effect$adj.P.Val<0.05   # over whole dataset
  l4s<- int_effect[sampledRows,]$adj.P.Val<0.05 # over just ground truth
  l2 <- main_effect$adj.P.Val<0.05
  l2s<- main_effect[sampledRows,]$adj.P.Val<0.05
  
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
    SR <- append(SR, (sum(int_effect[sampledRows[30*0:199+j],]$adj.P.Val<0.05)/200))
  }
  for (j in 1:30){
    SR <- append(SR, (sum(main_effect[sampledRows[30*0:199+j],]$adj.P.Val<0.05)/200))
  }
  
  sensLFCdf <- data.frame("rate"=SR,
                          "method"=rep("limma", 60),
                          "term" = rep(c("int", "main"), each=30),
                          "beta" = rep(params$lfc_fdr_sensitivity_per_lfc_figure, 2),
                          "rep" = i)
  limmaStats[[i]] <- list(fdr_df=FDRv, sen_df=SENv, sen_lfc_df=sensLFCdf)
}



##### A N A L Y S E    D E S E Q 2 #####
metadata <- data.frame("condition"=c("Ctrl", "Ctrl", "Ctrl", "Ctrl", "Trt", "Trt", "Trt", "Trt"), 
                       "genotype"=c("WT", "WT", "MU", "MU", "WT","WT", "MU", "MU"))
#rownames(metadata) <- colnames(sig2rep2v0)

deseq2Stats <- list()
for (i in 1:params$n_synth_dfs){
  randomDataFrame <- ground_truth_datasets[[i]]
  
  ddsr <-  DESeqDataSetFromMatrix(round(randomDataFrame), colData=metadata, design =~ genotype + condition + genotype:condition)
  ddsr$genotype = relevel(ddsr$genotype, "WT")
  ddsr <- DESeq(ddsr)
  resrAB = results(ddsr, name="genotypeMU.conditionTrt", independentFiltering = F)
  resrAB$pvalue[is.na(resrAB$pvalue)] <- 1
  resrAB$padj[is.na(resrAB$padj)] <- 1

  resrB = results(ddsr, contrast=c("condition","Trt","Ctrl"), independentFiltering = F)
  resrB$pvalue[is.na(resrB$pvalue)] <- 1
  resrB$padj[is.na(resrB$padj)] <- 1
  
  resrA = results(ddsr, contrast=c("genotype","MU","WT"), independentFiltering = F)
  resrA$pvalue[is.na(resrA$pvalue)] <- 1
  resrA$padj[is.na(resrA$padj)] <- 1

  
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
                          "beta" = rep(params$lfc_fdr_sensitivity_per_lfc_figure, 2),
                          "rep" = i)
  
  deseq2Stats[[i]] <- list(fdr_df=FDRv, sen_df=SENv, sen_lfc_df=sensLFCdf) # TODO: named list to improve readability and repeatability
}


# Define method names and metrics
methods <- c("LinSigStats", "limmaStats", "deseq2Stats")

# Create combined data frames using list comprehension
combine_data <- function(metric) {
  do.call(rbind, lapply(methods, function(method) {
    do.call(rbind, lapply(1:10, function(i) {
      df <- get(method)[[i]][[metric]]
      # Add replicate number column
      df$rep <- i
      df
    }))
  }))
}

# Create combined FDR and Sensitivity data frames
allFDRs <- combine_data("fdr_df")
allSENs <- combine_data("sen_df")

# Combine metrics
allfdrsens <- allFDRs
allfdrsens$Recall <- allSENs$Recall

# Add visualization parameters
allfdrsens$col <- rep(c("red", "green", "blue"), each = 20)
allfdrsens$shape <- rep(c(18, 19), 30)

allSensVarLFC <- combine_data("sen_lfc_df")


## FDR ~ Sensitivity plot
ggplot(data=allfdrsens, aes(x=Recall, y=1-FDR, fill=method)) +
  geom_point(aes(shape=term), size=2)+
  scale_shape_manual(values=c(21, 22, 24)) +
  coord_cartesian(ylim=c(0.92,0.99), xlim=c(0.6,0.99))+
  xlab("Sensitivity") +
  stat_ellipse(geom="polygon", level=0.95, aes(fill=method), alpha=0.25)


## Sensitivity ~ Fold Change plot
ggplot(data=allSensVarLFC[allSensVarLFC$beta>0,], aes(x=beta, y=rate, group=interaction(method, term), color=method, linetype=term))+
  geom_smooth(method="loess",se=T, span=0.7, aes(fill=method), alpha=0.3)+
  ylab("Sensitivity")+
  xlab("fold change") +
  coord_cartesian(ylim=c(0.0, 1))





################################################################################
################### R E G U L A T I O N    R E C O V E R Y #####################
################################################################################

#7 different types of regulations

#dataframe with all types of regulation.
#For each regulation compute the number of wrong logic/no logic/right logic?

  # f=params$lfc_reg_recovery_figure
  # 
  # A <-  c(0,0,f,f,0,0,f,f)
  # B <-  c(0,0,0,0,f,f,f,f)
  # AB<-  c(0,0,0,0,0,0,f,f)
  # BA<-  c(0,0,f,f,-f,-f,0,0)
  # AAB<- c(0,0,f,f,0,0,0,0)
  # BAB<- c(0,0,0,0,f,f,0,0)
  # ABAB<-c(0,0,f,f,f,f,f,f)
  # 
  # GTmatrix <- matrix(rep(c(A,B,AB,BA,AAB,BAB,ABAB),750), ncol=8,byrow = T)
  # sampledRows <- sample(2000:40000, 5250)
  # emptyMatrix <- matrix(0, nrow=nrow(rdf1), ncol=8)
  # emptyMatrix[sampledRows,] <- GTmatrix * sample(c(-1,1), 5250, replace=TRUE)

A <-  c(0,0,1,1,0,0,1,1)
B <-  c(0,0,0,0,1,1,1,1)
AB<-  c(0,0,0,0,0,0,1,1)
BA<-  c(0,0,1,1,-1,-1,0,0)
AAB<- c(0,0,1,1,0,0,0,0)
BAB<- c(0,0,0,0,1,1,0,0)
ABAB<-c(0,0,1,1,1,1,1,1)

gt_df <- generateGroundTruthGenes(n_ground_truth_genes = 5250, # 7 regulation types, each 750 times simulated
                         n_total_genes = 40000,
                         range=2000:40000, # counts are ordered, so by picking a range above 2000 you're avoiding genes that have a count of 1-10 or so.
                         fold_change = params$lfc_reg_recovery_figure,
                         regulation = c(A,B,AB,BA,AAB,BAB,ABAB))

reg_rec_datasets <- generate_multiple_datasets(source_df = cts, nrep = 2, size = params$n_genes_simulated)


ground_truth_datasets <- lapply(reg_rec_datasets, function(df) {
  scaled_df <- 2^gt_df$emptyMatrix * df
  rownames(scaled_df) <- as.character(1:nrow(scaled_df))
  return(scaled_df)
})

sampledRows <- gt_df$ground_truth_data$sampled_rows



## TODO: refactor this function because it contains a lot of repetitive parts
fit_intmodel <- function(counts, meanThr=10, pseudo=1, lowess=FALSE, nrep=2,
                         FDR_pval=0.05, sizeRDF=60000, structureDataFrame){
  
  normed <- MedianNorm(counts, count_threshold = meanThr, pseudo = pseudo)
  rats <- compute_ratios(normed, lowess_norm = F, structureDataFrame = structureDataFrame)
  fitstats <- ratios.fit(rats, CompThreshold = 1, n_rep=nrep) # compute model statistics
  
  RandomDF <- generate_synthetic_data(counts, size = 40000, nrep=2)
  
  message("Random Dataset succesfully created for FDR")
  #mediannorm for random dataset
  normedR <- MedianNorm(RandomDF)
  #compute ratios for random dataset
  ratsR <- compute_ratios(normedR, lowess_norm=lowess, structureDataFrame = params$strucDF)
  #fit model for random dataset
  message("computing random stats")
  statsR <- ratios.fit(ratsR, CompThreshold=1,n_rep=nrep)
  
  statsR$SIGA <- statsR[,6]<0.05 & statsR$R2>0.8 # 0.7 for adjusted R2, 0.8 for Multiple Rsquared
  statsR$SIGB <- statsR[,7]<0.05 & statsR$R2>0.8
  statsR$SIGAB <-statsR[,8]<0.05 & statsR$R2>0.8
  
  print(sum(statsR$SIGA))
  print(sum(statsR$SIGB))
  print(sum(statsR$SIGAB))
  
  FDRA <- statsR[statsR$SIGA,]
  FDRB <- statsR[statsR$SIGB,]
  FDRAB<- statsR[statsR$SIGAB,]
  
  FDRab <- FDRb <- FDRa <- c()
  for (i in 1:length(seq(0,4,0.001))){
    FDRab[i] <- mean(abs(FDRAB$X3)>(seq(0,4,0.001)[i]))
  }
  for (i in 1:length(seq(0,4,0.001))){
    FDRb[i] <- mean(abs(FDRB$X2)>(seq(0,4,0.001)[i]))
  }
  for (i in 1:length(seq(0,4,0.001))){
    FDRa[i] <- mean(abs(FDRA$X1)>(seq(0,4,0.001)[i]))
  }
  xs <- seq(0,4,0.001)
  FDRDF <- data.frame("Xs"=xs, "FDRa"=FDRa,"FDRb"=FDRb, "FDRab"=FDRab)
  print(paste("5% FDR B_threshold A:", seq(0,4,0.001)[which.min(abs(FDRa-0.05))]))
  print(paste("5% FDR B_threshold B:", seq(0,4,0.001)[which.min(abs(FDRb-0.05))]))
  print(paste("5% FDR B_threshold AB:",seq(0,4,0.001)[which.min(abs(FDRab-0.05))]))
  
  NN <- function(x,y,i){ #x:betas, y:FDR, i:numbertotest
    loc <- which.min(abs(x-abs(i)))
    return(y[loc])
  } #nearest neighbor function
  
  #compute FDR chance statistic
  fitstats$FDRA <- sapply(fitstats$X1, NN, x=xs, y=FDRa)
  fitstats$FDRB <- sapply(fitstats$X2, NN, x=xs, y=FDRb)
  fitstats$FDRAB <- sapply(fitstats$X3, NN, x=xs, y=FDRab)
  
  #output dataset
  return(fitstats)
}



##### LIMMA #####
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



for (i in ground_truth_datasets){
  print(head(i))
  
  # LinSig
  vregStats <- fit_intmodel(i, size=40000, structureDataFrame = params$strucDF)
  
  
  onlyRegulatedGenes <- data.frame(vregStats[as.character(sampledRows),])
  regs<-c("A","B", "A:B", "A+B", "B+AB", "A+AB", "A+B+AB")
  ground_truth_labels <- rep(regs, 750)
  onlyRegulatedGenes$true_class <- ground_truth_labels
  
  
  # only true discovery genes
  otruereg <- c(sum(onlyRegulatedGenes$FDRA<0.05 & onlyRegulatedGenes$FDRB>0.05 & onlyRegulatedGenes$FDRAB>0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1<0.05 & onlyRegulatedGenes$cov_x2>0.00 & onlyRegulatedGenes$cov_x3>0.00 & onlyRegulatedGenes$true_class=="B", na.rm=T), #
                #just B rep(1:7, 750)==2
                sum(onlyRegulatedGenes$FDRA>0.05 & onlyRegulatedGenes$FDRB<0.05 & onlyRegulatedGenes$FDRAB>0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1>0.00 & onlyRegulatedGenes$cov_x2<0.05 & onlyRegulatedGenes$cov_x3>0.00 & onlyRegulatedGenes$true_class=="A", na.rm=T),
                #just A+B rep(1:7, 750)==3
                sum(onlyRegulatedGenes$FDRA<0.05 & onlyRegulatedGenes$FDRB<0.05 & onlyRegulatedGenes$FDRAB>0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1<0.05 & onlyRegulatedGenes$cov_x2<0.05 & onlyRegulatedGenes$cov_x3>0 & onlyRegulatedGenes$true_class=="A+B", na.rm=T),
                #just AB rep(1:7, 750)==4
                sum(onlyRegulatedGenes$FDRA>0.05 & onlyRegulatedGenes$FDRB>0.05 & onlyRegulatedGenes$FDRAB<0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1>0.00 & onlyRegulatedGenes$cov_x2>0.00 & onlyRegulatedGenes$cov_x3<0.05 & onlyRegulatedGenes$true_class=="A:B", na.rm=T),
                #just A+AB rep(1:7, 750)==5
                sum(onlyRegulatedGenes$FDRA<0.05 & onlyRegulatedGenes$FDRB>0.05 & onlyRegulatedGenes$FDRAB<0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1<0.05 & onlyRegulatedGenes$cov_x2>0.00 & onlyRegulatedGenes$cov_x3<0.05 & onlyRegulatedGenes$true_class=="A+AB", na.rm=T),
                #just B+AB rep(1:7, 750)==6
                sum(onlyRegulatedGenes$FDRA>0.05 & onlyRegulatedGenes$FDRB<0.05 & onlyRegulatedGenes$FDRAB<0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1>0.00 & onlyRegulatedGenes$cov_x2<0.05 & onlyRegulatedGenes$cov_x3<0.05 & onlyRegulatedGenes$true_class=="B+AB", na.rm=T),
                #just A+B+AB rep(1:7, 750)==7
                sum(onlyRegulatedGenes$FDRA<0.05 & onlyRegulatedGenes$FDRB<0.05 & onlyRegulatedGenes$FDRAB<0.05 & onlyRegulatedGenes$R2>0.8 &
                      onlyRegulatedGenes$cov_x1<0.05 & onlyRegulatedGenes$cov_x2<0.05 & onlyRegulatedGenes$cov_x3<0.05 & onlyRegulatedGenes$true_class=="A+B+AB", na.rm=T))


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
  
  # limma
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
  
  litruereg<- c(sum(l2s$adj.P.Val<0.05 & l3s$adj.P.Val>0.05 & l4s$adj.P.Val>0.05 & onlyRegulatedGenes$true_class=='A'),
                sum(l2s$adj.P.Val>0.05 & l3s$adj.P.Val<0.05 & l4s$adj.P.Val>0.05 & onlyRegulatedGenes$true_class=='B'),
                sum(l2s$adj.P.Val<0.05 & l3s$adj.P.Val<0.05 & l4s$adj.P.Val>0.05 & onlyRegulatedGenes$true_class=='A+B'),
                sum(l2s$adj.P.Val>0.05 & l3s$adj.P.Val>0.05 & l4s$adj.P.Val<0.05 & onlyRegulatedGenes$true_class=='A:B'),
                sum(l2s$adj.P.Val<0.05 & l3s$adj.P.Val>0.05 & l4s$adj.P.Val<0.05 & onlyRegulatedGenes$true_class=='B+AB'),
                sum(l2s$adj.P.Val>0.05 & l3s$adj.P.Val<0.05 & l4s$adj.P.Val<0.05 & onlyRegulatedGenes$true_class=='A+AB'),
                sum(l2s$adj.P.Val<0.05 & l3s$adj.P.Val<0.05 & l4s$adj.P.Val<0.05 & onlyRegulatedGenes$true_class=='A+B+AB'))
  
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
  
  d2truereg <-c(sum(resrA[sampledRows,]$padj<0.05 & resrB[sampledRows,]$padj>0.05 & resrAB[sampledRows,]$padj>0.05 & onlyRegulatedGenes$true_class=='A'),
                sum(resrA[sampledRows,]$padj>0.05 & resrB[sampledRows,]$padj<0.05 & resrAB[sampledRows,]$padj>0.05 & onlyRegulatedGenes$true_class=='B'),
                sum(resrA[sampledRows,]$padj<0.05 & resrB[sampledRows,]$padj<0.05 & resrAB[sampledRows,]$padj>0.05 & onlyRegulatedGenes$true_class=='A+B'),
                sum(resrA[sampledRows,]$padj>0.05 & resrB[sampledRows,]$padj>0.05 & resrAB[sampledRows,]$padj<0.05 & onlyRegulatedGenes$true_class=='A:B'),
                sum(resrA[sampledRows,]$padj<0.05 & resrB[sampledRows,]$padj>0.05 & resrAB[sampledRows,]$padj<0.05 & onlyRegulatedGenes$true_class=='B+AB'),
                sum(resrA[sampledRows,]$padj>0.05 & resrB[sampledRows,]$padj<0.05 & resrAB[sampledRows,]$padj<0.05 & onlyRegulatedGenes$true_class=='A+AB'),
                sum(resrA[sampledRows,]$padj<0.05 & resrB[sampledRows,]$padj<0.05 & resrAB[sampledRows,]$padj<0.05 & onlyRegulatedGenes$true_class=='A+B+AB'))
  
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
                      "replicate"=rep(i, 63)
  )
  
  
  cTMFdf <- rbind(cTMFdf, TMFdf)
}


write.csv(cTMFdf, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/2sig_2rep_RegulationRecoveryData.csv")


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

