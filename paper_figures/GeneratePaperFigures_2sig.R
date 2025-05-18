# figure generation for 2 replicates


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
library(ggplot2)
source("~/Boston Internship/Github/Rsyn/paper_figures/data_standardization.R") # for MedianNorm function
source("~/Boston Internship/Github/Rsyn/paper_figures/gen_synthetic_data_helpers.R")
source("~/Boston Internship/cellsigsyn/R/compute_ratios.R") # for RNAseqLowess function
source("~/Boston Internship/Github/Rsyn/paper_figures/compute_lfc_thresholds.R")
source("~/Boston Internship/cellsigsyn/R/fit_model.R")

cts <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/IL6IL10combDF.csv", row.names=1)[,1:8]


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


###################################################################
# Generate multiple synthetic datasets
generate_multiple_datasets <- function(source_df, nrep = 2, size = 40000, n_datasets = 10) {
  replicate(n_datasets, generate_synthetic_data(source_df = source_df, nrep = nrep, size = size), simplify = FALSE)
}

# Create ground truth scaling matrix
create_ground_truth_scaling <- function(logFCs, n_genes_per_group = 3000) {
  logFCs_all <- rep(logFCs, 2)  # length = 2 × length(logFCs)
  # Repeat the AB pattern 3000 times for positive and negative effects
  scaling_AB_pos <- matrix(rep(c(0, 0, 1, 1, 0, 0, 0, 0), n_genes_per_group), ncol = 8, byrow = TRUE)
  scaling_AB_neg <- matrix(rep(c(0, 0, -1, -1, 0, 0, 0, 0), n_genes_per_group), ncol = 8, byrow = TRUE)
    design_matrix <- rbind(scaling_AB_pos, scaling_AB_neg)  # 6000 × 8
    ground_truth_scaling_matrix <- design_matrix * logFCs_all
  
  return(ground_truth_scaling_matrix)
}

# Apply ground truth scaling to datasets
apply_ground_truth <- function(datasets, ground_truth_scaling, total_genes = 40000, signal_genes = 6000) {
  empty_matrix <- matrix(0, nrow = total_genes, ncol = 8)
  sampled_rows <- sample(1:total_genes, signal_genes)
  print(dim(ground_truth_scaling))         # Should be 6000 x 8
  print(length(sampled_rows)    )           # Should be 6000
  empty_matrix[sampled_rows, ] <- ground_truth_scaling
  
  # Apply scaling to each dataset
  scaled_datasets <- lapply(datasets, function(df) {
    scaled_df <- 2^empty_matrix * df
    rownames(scaled_df) <- as.character(1:nrow(scaled_df))
    return(scaled_df)
  })
  
  list(datasets = scaled_datasets, sampled_rows = as.character(sampled_rows))
}

# Main execution
logFCs <- log2(seq(1.5, 3.25, by = 0.125))
logFCs_all <- rep(logFCs, 2)

datasets <- generate_multiple_datasets(source_df = cts, nrep = 2, size = params$n_genes_simulated)
ground_truth <- create_ground_truth_scaling(logFCs_all)
ground_truth_datasets <- apply_ground_truth(datasets, ground_truth)

# Access outputs
sampledRows <- ground_truth_datasets$sampled_rows
all10DFs <- ground_truth_datasets$datasets



###########################################
########## S T A T S ######################
###########################################


##### A N A L Y S E    L I N S I G #####
LinSigStats <- list()

# compute recommended thresholds -> later replace with FDR<0.05 calculation. Takes a long time though for each df...
compute_lfc_thresholds(sig2rep2v3,nrep=2, size=100000, nclus=7, lowessn=0)



for (i in 1:params$n_synth_dfs){
  randomDataFrame <- all10DFs[[i]]
  GTnorm <- MedianNorm(randomDataFrame)
  GTratios <- compute_ratios(GTnorm, lowess = 0, structureDataFrame = params$strucDF)
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

##### A N A L Y S E     L I M M A  #####
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
for (i in 1:params$n_synth_dfs){
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



##### A N A L Y S E    D E S E Q 2 #####
library(DESeq2)
metadata <- data.frame("condition"=c("Ctrl", "Ctrl", "Ctrl", "Ctrl", "Trt", "Trt", "Trt", "Trt"), 
                       "genotype"=c("WT", "WT", "MU", "MU", "WT","WT", "MU", "MU"))
rownames(metadata) <- colnames(sig2rep2v0)

deseq2Stats <- list()
for (i in 1:params$n_synth_dfs){
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
  coord_cartesian(ylim=c(0.0, 1))

#write.csv(allSensVarLFC, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/SENS_LFC_plotdata_variableLFC1.5_4.csv")
allSensVarLFC <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/2sig_data/SENS_LFC_plotdata_variableLFC1.5_4.csv")

################################################################################
################### R E G U L A T I O N    R E C O V E R Y #####################
################################################################################

#7 different types of regulations

#dataframe with all types of regulation.
#For each regulation compute the number of wrong logic/no logic/right logic?

f=params$lfc_reg_recovery_figure

A <-  c(0,0,f,f,0,0,f,f)
B <-  c(0,0,0,0,f,f,f,f)
AB<-  c(0,0,0,0,0,0,f,f)
BA<-  c(0,0,f,f,-f,-f,0,0)
AAB<- c(0,0,f,f,0,0,0,0)
BAB<- c(0,0,0,0,f,f,0,0)
ABAB<-c(0,0,f,f,f,f,f,f)


GTmatrix <- matrix(rep(c(A,B,AB,BA,AAB,BAB,ABAB),750), ncol=8,byrow = T)
sampledRows <- sample(2000:40000, 5250)

emptyMatrix <- matrix(0, nrow=nrow(rdf1), ncol=8)
emptyMatrix[sampledRows,] <- GTmatrix * sample(c(-1,1), 5250, replace=TRUE)


rdf1 <- generate_synthetic_data(source_df = cts, nrep=2, size=40000)
rdf2 <- generate_synthetic_data(source_df = cts,nrep=2, size=40000)
rdf3 <- generate_synthetic_data(source_df = cts,nrep=2, size=40000)
rdf4 <- generate_synthetic_data(source_df = cts,nrep=2, size=40000)
rdf5 <- generate_synthetic_data(source_df = cts,nrep=2, size=40000)
rdf6 <- generate_synthetic_data(source_df = cts,nrep=2, size=40000)
rdf7 <- generate_synthetic_data(source_df = cts,nrep=2, size=40000)
rdf8 <- generate_synthetic_data(source_df = cts,nrep=2, size=40000)
rdf9 <- generate_synthetic_data(source_df = cts,nrep=2, size=40000)
rdf0 <- generate_synthetic_data(source_df = cts,nrep=2, size=40000)



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

