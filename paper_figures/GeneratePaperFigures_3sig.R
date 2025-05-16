library(dplyr)
library(MASS)
library(cellsigsyn)
library(parallel)
library(DESeq2)
library(limma)
library(edgeR)
setwd("C:/Users/HB/OneDrive/Documents/Boston Internship/")
cts <- read.csv("IL6IL10combDF.csv", row.names=1)[,1:8]

rdf3_0 <- gen_randomstats(cts,nrep=4, size=40000, onlyDF=1)
rdf3_1 <- gen_randomstats(cts,nrep=4, size=40000, onlyDF=1)
rdf3_2 <- gen_randomstats(cts,nrep=4, size=40000, onlyDF=1)
rdf3_3 <- gen_randomstats(cts,nrep=4, size=40000, onlyDF=1)
rdf3_4 <- gen_randomstats(cts,nrep=4, size=40000, onlyDF=1)
rdf3_5 <- gen_randomstats(cts,nrep=4, size=40000, onlyDF=1)
rdf3_6 <- gen_randomstats(cts,nrep=4, size=40000, onlyDF=1)
rdf3_7 <- gen_randomstats(cts,nrep=4, size=40000, onlyDF=1)
rdf3_8 <- gen_randomstats(cts,nrep=4, size=40000, onlyDF=1)
rdf3_9 <- gen_randomstats(cts,nrep=4, size=40000, onlyDF=1)


#   B   A   C   AB   BC   AC   ABC
c / p / l / t / pl / tp / tl / tpl

## Generate 2 replicate, 3 signal dataset with 6 different types of regulation
f<-1 # 2^{1} fold change simulation
signalLogic3 <- c( c(0,0,0,0,f,f,0,0,f,f,0,0,f,f,f,f),    # A 
                c(0,0,f,f,f,f,-f,-f,f+f,f+f,0,0,0,0,f,f), # A+B+C
                c(0,0,0,0,0,0,0,0,f,f,0,0,0,0,f,f),       # A:B
                c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,f,f),       # A:B:C
                c(0,0,0,0,f,f,0,0,0,0,0,0,f,f,0,0),       # A:B+A
                c(0,0,0,0,f,f,0,0,0,0,0,0,f,f,f,f) )      # A:B:C+A:B+A

groundTruth3sig <- matrix(rep(c(signalLogic3, -signalLogic3), 450), ncol=16, byrow=T) # 5400 GT genes. All ground truth logic is generated with positive and negative regulation.
sampledGenes3sig <- sample(2000:40000, 5400)

#vreg3sig <- rdf30
#vreg3sig[sampledGenes3sig,] <- rdf30[sampledGenes3sig,] * 2^groundTruth3sig
#vreg3sig[sampledGenes3sig,]
#rownames(vreg3sig) <-as.character(1:nrow(vreg3sig))

emptyMatrix <- matrix(0, nrow=nrow(rdf3_0), ncol=16)
emptyMatrix[sampledGenes3sig,] <- groundTruth3sig


sig3rep2v0<-(2^(emptyMatrix)*rdf3_0)
sig3rep2v1<-(2^(emptyMatrix)*rdf3_1)
sig3rep2v2<-(2^(emptyMatrix)*rdf3_2)
sig3rep2v3<-(2^(emptyMatrix)*rdf3_3)
sig3rep2v4<-(2^(emptyMatrix)*rdf3_4)
sig3rep2v5<-(2^(emptyMatrix)*rdf3_5)
sig3rep2v6<-(2^(emptyMatrix)*rdf3_6)
sig3rep2v7<-(2^(emptyMatrix)*rdf3_7)
sig3rep2v8<-(2^(emptyMatrix)*rdf3_8)
sig3rep2v9<-(2^(emptyMatrix)*rdf3_9)

rownames(sig3rep2v0) <- as.character(1:nrow(sig3rep2v0))
rownames(sig3rep2v1) <- as.character(1:nrow(sig3rep2v1))
rownames(sig3rep2v2) <- as.character(1:nrow(sig3rep2v2))
rownames(sig3rep2v3) <- as.character(1:nrow(sig3rep2v3))
rownames(sig3rep2v4) <- as.character(1:nrow(sig3rep2v4))
rownames(sig3rep2v5) <- as.character(1:nrow(sig3rep2v5))
rownames(sig3rep2v6) <- as.character(1:nrow(sig3rep2v6))
rownames(sig3rep2v7) <- as.character(1:nrow(sig3rep2v7))
rownames(sig3rep2v8) <- as.character(1:nrow(sig3rep2v8))
rownames(sig3rep2v9) <- as.character(1:nrow(sig3rep2v9))

write.csv(sig3rep2v0, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf0_3sig_2rep.csv")
write.csv(sig3rep2v1, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf1_3sig_2rep.csv")
write.csv(sig3rep2v2, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf2_3sig_2rep.csv")
write.csv(sig3rep2v3, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf3_3sig_2rep.csv")
write.csv(sig3rep2v4, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf4_3sig_2rep.csv")
write.csv(sig3rep2v5, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf5_3sig_2rep.csv")
write.csv(sig3rep2v6, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf6_3sig_2rep.csv")
write.csv(sig3rep2v7, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf7_3sig_2rep.csv")
write.csv(sig3rep2v8, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf8_3sig_2rep.csv")
write.csv(sig3rep2v9, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf9_3sig_2rep.csv")

sampledGenes3sig <- as.character(sampledGenes3sig)
write(sampledGenes3sig, "./LinSigPaper/3sig_data/sampledGenes3sig.txt") # vector with row names of ground truth genes

# sig3rep2v0 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf0_3sig_2rep.csv", row.names = 1)
# sig3rep2v1 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf1_3sig_2rep.csv", row.names = 1)
# sig3rep2v2 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf2_3sig_2rep.csv", row.names = 1)
# sig3rep2v3 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf3_3sig_2rep.csv", row.names = 1)
# sig3rep2v4 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf4_3sig_2rep.csv", row.names = 1)
# sig3rep2v5 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf5_3sig_2rep.csv", row.names = 1)
# sig3rep2v6 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf6_3sig_2rep.csv", row.names = 1)
# sig3rep2v7 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf7_3sig_2rep.csv", row.names = 1)
# sig3rep2v8 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf8_3sig_2rep.csv", row.names = 1)
# sig3rep2v9 <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/rdf9_3sig_2rep.csv", row.names = 1)
# 
# sampledGenes3sig <- read.table("./LinSigPaper/3sig_data/sampledGenes3sig.txt") # vector with row names of ground truth genes
# sampledGenes3sig <- sampledGenes3sig$V1
# sampledGenes3sig <- as.character(sampledGenes3sig)



### FUNCTIONS FOR DECONVOLUTING
deseq2method3 <- function(vreg3sig, sampledGenes3sig){ #inputdf vreg3sig / index Ground TRuth: sampledGenes3sig
  
  metadata <- data.frame("condition"=rep(c("ctrl", "pH","ctrl","ctrl", "pH", "pH","ctrl","pH"),each=2), # made up random signals
                         "genotype"=rep(c("WT","WT","MU","WT","MU","WT","MU","MU"),each=2),
                         "cyto"=rep(c("N", "N", "N", "Y", "N", "Y", "Y", "Y"), each=2) )
  rownames(metadata) <- colnames(vreg3sig)
  
  ddsr <-  DESeqDataSetFromMatrix(round(vreg3sig), colData=metadata, design =~ genotype + condition + cyto + genotype:condition+ genotype:cyto + condition:cyto + condition:cyto:genotype)
  
  ddsr$genotype = relevel(ddsr$genotype, "WT")
  ddsr <- DESeq(ddsr)
  
  # edit for 3 signals
  # B
  resrpH = results(ddsr, contrast=c("condition","pH","ctrl"), independentFiltering = F)
  resrpH$pvalue[is.na(resrpH$pvalue)] <- 1
  resrpH$padj[is.na(resrpH$padj)] <- 1
  sum(resrpH$padj<0.05, na.rm=T)
  # A
  resrLPS = results(ddsr, contrast=c("genotype","MU","WT"), independentFiltering = F)
  resrLPS$pvalue[is.na(resrLPS$pvalue)] <- 1
  resrLPS$padj[is.na(resrLPS$padj)] <- 1
  sum(resrLPS$padj<0.05, na.rm=T)
  # C
  resrTNF = results(ddsr, contrast=c("cyto","Y","N"), independentFiltering = F)
  resrTNF$pvalue[is.na(resrTNF$pvalue)] <- 1
  resrTNF$padj[is.na(resrTNF$padj)] <- 1
  sum(resrTNF$padj<0.05, na.rm=T)
  # A:B
  resrpHLPS = results(ddsr, name="genotypeMU.conditionpH", independentFiltering = F)
  resrpHLPS$pvalue[is.na(resrpHLPS$pvalue)] <- 1
  resrpHLPS$padj[is.na(resrpHLPS$padj)] <- 1
  sum(resrpHLPS$padj<0.05, na.rm=T)
  # B:C
  resrLPSTNF = results(ddsr, name="conditionpH.cytoY", independentFiltering = F)
  resrLPSTNF$pvalue[is.na(resrLPSTNF$pvalue)] <- 1
  resrLPSTNF$padj[is.na(resrLPSTNF$padj)] <- 1
  sum(resrLPSTNF$padj<0.05, na.rm=T)
  # A:C
  resrpHTNF = results(ddsr, name="genotypeMU.cytoY", independentFiltering = F)
  resrpHTNF$pvalue[is.na(resrpHTNF$pvalue)] <- 1
  resrpHTNF$padj[is.na(resrpHTNF$padj)] <- 1
  sum(resrpHTNF$padj<0.05, na.rm=T)
  # A:B:C
  resrpHLPSTNF = results(ddsr, name="genotypeMU.conditionpH.cytoY", independentFiltering = F)
  resrpHLPSTNF$pvalue[is.na(resrpHLPSTNF$pvalue)] <- 1
  resrpHLPSTNF$padj[is.na(resrpHLPSTNF$padj)] <- 1
  sum(resrpHLPSTNF$pvalue<0.05, na.rm=T)
  
  return(list(resrpH, resrLPS, resrTNF, resrpHLPS, resrLPSTNF, resrpHTNF, resrpHLPSTNF))
}

limmaMethod3 <- function(vreg3sig, sampledGenes3sig){
  
  #cond <- c("ctrl", "pH", "LPS", "TNF","pHLPS", "TNFpH", "TNFLPS", "LPSpHTNF")
  stim <- factor(rep(c("ctrl", "LPS", "ctrl", "ctrl", "LPS", "LPS", "ctrl", "LPS"),each=2), levels=c("ctrl", "LPS"))
  pH <- factor(rep(c("H","H","L", "H",  "L", "H", "L", "L"),each=2), levels=c("H", "L")) # high pH / low pH
  gen <- factor(rep(c("N","N","N","Y", "N","Y","Y", "Y"),each=2), levels=c("N","Y")) # Yes TNF , No TNF
  
  design <- model.matrix(~stim*pH*gen)
  
  dge <- DGEList(counts=vreg3sig)
  dge<- calcNormFactors(dge)
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  
  l2 <- (topTable(fit, number=Inf,coef=2, sort.by = "none")) # LPS 
  #l2s<- (topTable(fit, number=Inf,coef=2)[sampledGenes3sig,]$adj.P.Val<0.05)
  l3 <- (topTable(fit, number=Inf,coef=3, sort.by = "none")) # pH 
  #l3s<- (topTable(fit, number=Inf,coef=3)[sampledGenes3sig,]$adj.P.Val<0.05)
  l4 <- (topTable(fit, number=Inf,coef=4, sort.by = "none"))   # TNF over whole dataset
  #l4s<- (topTable(fit, number=Inf,coef=4)[sampledGenes3sig,]$adj.P.Val<0.05) # over just ground truth
  l5 <- (topTable(fit, number=Inf,coef=5, sort.by = "none")) #LPS:pH
  #l5s<- (topTable(fit, number=Inf,coef=5)[sampledGenes3sig,]$adj.P.Val<0.05)
  l6 <- (topTable(fit, number=Inf,coef=6, sort.by = "none")) # LPS:TNF
  #l6s<- (topTable(fit, number=Inf,coef=6)[sampledGenes3sig,]$adj.P.Val<0.05)
  l7 <- (topTable(fit, number=Inf,coef=7, sort.by = "none")) # pH:TNF
  #l7s<- (topTable(fit, number=Inf,coef=7)[sampledGenes3sig,]$adj.P.Val<0.05)
  l8 <- (topTable(fit, number=Inf,coef=8, sort.by = "none")) # LPS:pH:TNF
  #l8s<- (topTable(fit, number=Inf,coef=8)[sampledGenes3sig,]$adj.P.Val<0.05)
  
  return(list(l2,l3,l4,l5,l6,l7,l8))
}
  

parfunct3 <- function(ratios, nclus){
  clust <- makeCluster(nclus)
  decovec <- rep(c(0,1,0,0,0,0,0,
                   0,1,0,1,0,0,0,
                   1,0,0,0,0,0,0,
                   1,0,0,1,0,0,0,
                   1,1,0,1,0,0,0,
                   0,0,1,0,0,0,0,
                   0,1,0,0,0,1,0,
                   1,0,0,0,1,0,0,
                   0,0,1,0,0,1,0,
                   0,0,1,0,1,0,0,
                   1,1,1,1,1,1,1,
                   1,0,0,1,1,0,1,
                   0,1,0,1,0,1,1,
                   0,0,1,0,1,1,1,
                   1,0,1,0,1,0,0,
                   0,1,1,0,0,1,0), 2)
  X <- matrix(decovec, ncol=7, byrow=T)
  clusterExport(clust, "X", envir=environment())
  mstats <- parApply(clust, ratios, 1, function(x){
    r<- summary(lm(x ~ X))$r.squared
    cov<- diag(vcov(lm(x ~ X)))
    return(c(cov,r))
  })
  stopCluster(clust)
  return(t(mstats))
}

ownfunc <- function(vreg3sig){
  real_dat <- MedianNorm(vreg3sig)
  reps <- 2
  dfd <- data.frame("c"=1:reps, "pH"=(1:reps+reps),
                    "LPS"=(1:reps+2*reps), "TNF"=(1:reps+3*reps),
                    "pHLPS"=(1:reps+4*reps), "TNFpH"=(1:reps+5*reps),
                    "TNFLPS"=(1:reps+6*reps),"TNFpHLPS"=(1:reps+7*reps))
  
  LogRatios <- matrix(nrow=length(real_dat[,1]), ncol=16*reps)
  
  for (i in 1:reps){
    LogRatios[,1+(16*(i-1))] <- log2(real_dat[,dfd$LPS[i]] / real_dat[,dfd$c[i]])
    LogRatios[,2+(16*(i-1))] <- log2(real_dat[,dfd$pHLPS[i]] / real_dat[,dfd$pH[i]])
    LogRatios[,3+(16*(i-1))] <- log2(real_dat[,dfd$pH[i]] / real_dat[,dfd$c[i]])
    LogRatios[,4+(16*(i-1))] <- log2(real_dat[,dfd$pHLPS[i]] / real_dat[,dfd$LPS[i]])
    LogRatios[,5+(16*(i-1))] <- log2(real_dat[,dfd$pHLPS[i]] / real_dat[,dfd$c[i]])
    LogRatios[,6+(16*(i-1))] <- log2(real_dat[,dfd$TNF[i]] / real_dat[,dfd$c[i]])
    LogRatios[,7+(16*(i-1))] <- log2(real_dat[,dfd$TNFLPS[i]] / real_dat[,dfd$TNF[i]])
    LogRatios[,8+(16*(i-1))] <- log2(real_dat[,dfd$TNFpH[i]] / real_dat[,dfd$TNF[i]])
    LogRatios[,9+(16*(i-1))] <- log2(real_dat[,dfd$TNFLPS[i]] / real_dat[,dfd$LPS[i]])
    LogRatios[,10+(16*(i-1))]<- log2(real_dat[,dfd$TNFpH[i]] / real_dat[,dfd$pH[i]])
    LogRatios[,11+(16*(i-1))]<- log2(real_dat[,dfd$TNFpHLPS[i]] / real_dat[,dfd$c[i]])
    LogRatios[,12+(16*(i-1))]<- log2(real_dat[,dfd$TNFpHLPS[i]] / real_dat[,dfd$TNFLPS[i]])
    LogRatios[,13+(16*(i-1))]<- log2(real_dat[,dfd$TNFpHLPS[i]] / real_dat[,dfd$TNFpH[i]])
    LogRatios[,14+(16*(i-1))]<- log2(real_dat[,dfd$TNFpHLPS[i]] / real_dat[,dfd$pHLPS[i]])
    LogRatios[,15+(16*(i-1))]<- log2(real_dat[,dfd$TNFpH[i]] / real_dat[,dfd$c[i]])
    LogRatios[,16+(16*(i-1))]<- log2(real_dat[,dfd$TNFLPS[i]] / real_dat[,dfd$c[i]])
  }
  
  decovec <- c(0,1,0,0,0,0,0,
               0,1,0,1,0,0,0,
               1,0,0,0,0,0,0,
               1,0,0,1,0,0,0,
               1,1,0,1,0,0,0,
               0,0,1,0,0,0,0,
               0,1,0,0,0,1,0,
               1,0,0,0,1,0,0,
               0,0,1,0,0,1,0,
               0,0,1,0,1,0,0,
               1,1,1,1,1,1,1,
               1,0,0,1,1,0,1,
               0,1,0,1,0,1,1,
               0,0,1,0,1,1,1,
               1,0,1,0,1,0,0,
               0,1,1,0,0,1,0)
  
  Ratios <- LogRatios
  
  X <- matrix(c(decovec,decovec), byrow = T, ncol=7)
  multiplefit_r <- lm(t(Ratios) ~ X)
  B_r <- t(multiplefit_r$coefficients)
  
  stat43sig <-parfunct3(LogRatios, nclus=7)
  X <- matrix(c(decovec), byrow = T, ncol=7)
  B3 <- lm(t(LogRatios)~ rbind(X,X))
  B3r <- t(B3$coefficients)
  
  GT3sigstat <- data.frame(cbind(B3r, stat43sig))
  colnames(GT3sigstat) <- c("int", "x1", "x2", "x3", "x4", "x5", "x6", "x7",
                            "p_int", "p_x1", "p_x2", "p_x3", "p_x4", "p_x5", "p_x6", "p_x7",
                            "R2")
  rownames(GT3sigstat) <- rownames(real_dat)
  return(list(GT3sigstat, dfd))
}

fdrcalcfunc <- function(rdfrDF, GT3sigstat, reps, dfd){
  #rdrfDF is dataset without any groundtruth genes. so rdf3_0
  rownames(rdfrDF) <- as.character(1:nrow(rdfrDF))
  normed3r <- MedianNorm(rdfrDF)
  #same analysis with rdfr31
  # for each 7 coeffficients distribution. 1 will be 3int, avg of 3 2int, avg of 3 single effect
  
  LogRatiosR <- matrix(nrow=length(normed3r[,1]), ncol=16*reps)
  rand_dat <- normed3r
  for (i in 1:reps){
    LogRatiosR[,1+(16*(i-1))] <- log2(rand_dat[,dfd$LPS[i]] / rand_dat[,dfd$c[i]])
    LogRatiosR[,2+(16*(i-1))] <- log2(rand_dat[,dfd$pHLPS[i]] / rand_dat[,dfd$pH[i]])
    LogRatiosR[,3+(16*(i-1))] <- log2(rand_dat[,dfd$pH[i]] / rand_dat[,dfd$c[i]])
    LogRatiosR[,4+(16*(i-1))] <- log2(rand_dat[,dfd$pHLPS[i]] / rand_dat[,dfd$LPS[i]])
    LogRatiosR[,5+(16*(i-1))] <- log2(rand_dat[,dfd$pHLPS[i]] / rand_dat[,dfd$c[i]])
    LogRatiosR[,6+(16*(i-1))] <- log2(rand_dat[,dfd$TNF[i]] / rand_dat[,dfd$c[i]])
    LogRatiosR[,7+(16*(i-1))] <- log2(rand_dat[,dfd$TNFLPS[i]] / rand_dat[,dfd$TNF[i]])
    LogRatiosR[,8+(16*(i-1))] <- log2(rand_dat[,dfd$TNFpH[i]] / rand_dat[,dfd$TNF[i]])
    LogRatiosR[,9+(16*(i-1))] <- log2(rand_dat[,dfd$TNFLPS[i]] / rand_dat[,dfd$LPS[i]])
    LogRatiosR[,10+(16*(i-1))]<- log2(rand_dat[,dfd$TNFpH[i]] / rand_dat[,dfd$pH[i]])
    LogRatiosR[,11+(16*(i-1))]<- log2(rand_dat[,dfd$TNFpHLPS[i]] / rand_dat[,dfd$c[i]])
    LogRatiosR[,12+(16*(i-1))]<- log2(rand_dat[,dfd$TNFpHLPS[i]] / rand_dat[,dfd$TNFLPS[i]])
    LogRatiosR[,13+(16*(i-1))]<- log2(rand_dat[,dfd$TNFpHLPS[i]] / rand_dat[,dfd$TNFpH[i]])
    LogRatiosR[,14+(16*(i-1))]<- log2(rand_dat[,dfd$TNFpHLPS[i]] / rand_dat[,dfd$pHLPS[i]])
    LogRatiosR[,15+(16*(i-1))]<- log2(rand_dat[,dfd$TNFpH[i]] / rand_dat[,dfd$c[i]])
    LogRatiosR[,16+(16*(i-1))]<- log2(rand_dat[,dfd$TNFLPS[i]] / rand_dat[,dfd$c[i]])
  }
  
  stat43sigR <-parfunct3(LogRatiosR, nclus=7)
  
  decovec <- c(0,1,0,0,0,0,0,
               0,1,0,1,0,0,0,
               1,0,0,0,0,0,0,
               1,0,0,1,0,0,0,
               1,1,0,1,0,0,0,
               0,0,1,0,0,0,0,
               0,1,0,0,0,1,0,
               1,0,0,0,1,0,0,
               0,0,1,0,0,1,0,
               0,0,1,0,1,0,0,
               1,1,1,1,1,1,1,
               1,0,0,1,1,0,1,
               0,1,0,1,0,1,1,
               0,0,1,0,1,1,1,
               1,0,1,0,1,0,0,
               0,1,1,0,0,1,0)
  
  X <- matrix(c(decovec,decovec), byrow = T, ncol=7)
  
  B3R <- lm(t(LogRatiosR)~ X)
  B3rR <- data.frame(t(B3R$coefficients))
  rownames(B3rR) <- rownames(normed3r)
  dim(stat43sigR)
  dim(B3rR)
  
  randomdfdr3 <- data.frame(cbind(B3rR, stat43sigR))
  colnames(randomdfdr3) <- c("int", "x1", "x2", "x3", "x4", "x5", "x6", "x7",
                             "p_int", "p_x1", "p_x2", "p_x3", "p_x4", "p_x5", "p_x6", "p_x7",
                             "R2")
  
  x1fdr <- abs(randomdfdr3[(randomdfdr3$p_x1<0.05 & randomdfdr3$R2>0.8),]$x1)
  x2fdr <- abs(randomdfdr3[(randomdfdr3$p_x2<0.05 & randomdfdr3$R2>0.8),]$x2)
  x3fdr <- abs(randomdfdr3[(randomdfdr3$p_x3<0.05 & randomdfdr3$R2>0.8),]$x3)
  x4fdr <- abs(randomdfdr3[(randomdfdr3$p_x4<0.05 & randomdfdr3$R2>0.8),]$x4)
  x5fdr <- abs(randomdfdr3[(randomdfdr3$p_x5<0.05 & randomdfdr3$R2>0.8),]$x5)
  x6fdr <- abs(randomdfdr3[(randomdfdr3$p_x6<0.05 & randomdfdr3$R2>0.8),]$x6)
  x7fdr <- abs(randomdfdr3[(randomdfdr3$p_x7<0.05 & randomdfdr3$R2>0.8),]$x7)
  
  #assign FDR value to each beta from GT table
  
  #assign FDR value func
  aFDR <- function(Beta, FDRdist){
    return(mean(FDRdist>abs(Beta)))
  }
  
  aFDR(0.48, x1fdr)
  fdr1 <- sapply(GT3sigstat$x1, aFDR, FDRdist=x1fdr)
  fdr2 <- sapply(GT3sigstat$x2, aFDR, FDRdist=x2fdr)
  fdr3 <- sapply(GT3sigstat$x3, aFDR, FDRdist=x3fdr)
  fdr4 <- sapply(GT3sigstat$x4, aFDR, FDRdist=x4fdr)
  fdr5 <- sapply(GT3sigstat$x5, aFDR, FDRdist=x5fdr)
  fdr6 <- sapply(GT3sigstat$x6, aFDR, FDRdist=x6fdr)
  fdr7 <- sapply(GT3sigstat$x7, aFDR, FDRdist=x7fdr)
  
  names(fdr1) <- names(fdr2) <- names(fdr3) <- names(fdr4) <- names(fdr5) <- names(fdr6) <- names(fdr7) <- rownames(GT3sigstat)
  
  return(list(fdr1,fdr2,fdr3,fdr4,fdr5,fdr6,fdr7))
}



############################
#### REPEAT FROM HERE ######
############################

sig3rep2 <- sig3rep2v9

limmalist <- limmaMethod3(sig3rep2, sampledGenes3sig)

limpH <- limmalist[[1]]
limLPS <- limmalist[[2]]
limTNF <- limmalist[[3]]
limpHLPS <- limmalist[[4]]
limLPSTNF <- limmalist[[5]]
limpHTNF <- limmalist[[6]]
limpHLPSTNF <- limmalist[[7]]


LTD3 <- c(sum(limLPS[sampledGenes3sig,]$adj.P.Val<0.05 & limpH[sampledGenes3sig,]$adj.P.Val>0.05 & limTNF[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limpHLPS[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 &  rep(1:6, 900)==1),
          #A+B+C
          sum(limLPS[sampledGenes3sig,]$adj.P.Val<0.05 & limpH[sampledGenes3sig,]$adj.P.Val<0.05 & limTNF[sampledGenes3sig,]$adj.P.Val<0.05 &
                limpHTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limpHLPS[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 &  rep(1:6, 900)==2),
          #A:B
          sum(limLPS[sampledGenes3sig,]$adj.P.Val>0.05 & limpH[sampledGenes3sig,]$adj.P.Val>0.05 & limTNF[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limpHLPS[sampledGenes3sig,]$adj.P.Val<0.05 &
                limpHLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 &  rep(1:6, 900)==3),
          #A:B:C
          sum(limLPS[sampledGenes3sig,]$adj.P.Val>0.05 & limpH[sampledGenes3sig,]$adj.P.Val>0.05 & limTNF[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limpHLPS[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHLPSTNF[sampledGenes3sig,]$adj.P.Val<0.05 &  rep(1:6, 900)==4),
          #A:B + A
          sum(limLPS[sampledGenes3sig,]$adj.P.Val<0.05 & limpH[sampledGenes3sig,]$adj.P.Val>0.05 & limTNF[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limpHLPS[sampledGenes3sig,]$adj.P.Val<0.05 &
                limpHLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 &  rep(1:6, 900)==5),
          #A:B:C + A:B + A
          sum(limLPS[sampledGenes3sig,]$adj.P.Val<0.05 & limpH[sampledGenes3sig,]$adj.P.Val>0.05 & limTNF[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limpHLPS[sampledGenes3sig,]$adj.P.Val<0.05 &
                limpHLPSTNF[sampledGenes3sig,]$adj.P.Val<0.05 &  rep(1:6, 900)==6)
)

#TrueMis Disovery
#A
LTMD3<- c(sum(limLPS[sampledGenes3sig,]$adj.P.Val<0.05 & limpH[sampledGenes3sig,]$adj.P.Val>0.05 & limTNF[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limpHLPS[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 ),
          #A+B+C
          sum(limLPS[sampledGenes3sig,]$adj.P.Val<0.05 & limpH[sampledGenes3sig,]$adj.P.Val<0.05 & limTNF[sampledGenes3sig,]$adj.P.Val<0.05 &
                limpHTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limpHLPS[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05),
          #A:B
          sum(limLPS[sampledGenes3sig,]$adj.P.Val>0.05 & limpH[sampledGenes3sig,]$adj.P.Val>0.05 & limTNF[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limpHLPS[sampledGenes3sig,]$adj.P.Val<0.05 &
                limpHLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05),
          #A:B:C
          sum(limLPS[sampledGenes3sig,]$adj.P.Val>0.05 & limpH[sampledGenes3sig,]$adj.P.Val>0.05 & limTNF[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limpHLPS[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHLPSTNF[sampledGenes3sig,]$adj.P.Val<0.05),
          #A:B + A
          sum(limLPS[sampledGenes3sig,]$adj.P.Val<0.05 & limpH[sampledGenes3sig,]$adj.P.Val>0.05 & limTNF[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limpHLPS[sampledGenes3sig,]$adj.P.Val<0.05 &
                limpHLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05),
          #A:B:C + A:B + A
          sum(limLPS[sampledGenes3sig,]$adj.P.Val<0.05 & limpH[sampledGenes3sig,]$adj.P.Val>0.05 & limTNF[sampledGenes3sig,]$adj.P.Val>0.05 &
                limpHTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limLPSTNF[sampledGenes3sig,]$adj.P.Val>0.05 & limpHLPS[sampledGenes3sig,]$adj.P.Val<0.05 &
                limpHLPSTNF[sampledGenes3sig,]$adj.P.Val<0.05)
)

#All Disovery
#A
LAD3 <- c(sum(limLPS$adj.P.Val<0.05 & limpH$adj.P.Val>0.05 & limTNF$adj.P.Val>0.05 &
                limpHTNF$adj.P.Val>0.05 & limLPSTNF$adj.P.Val>0.05 & limpHLPS$adj.P.Val>0.05 &
                limpHLPSTNF$adj.P.Val>0.05 ),
          #A+B+C
          sum(limLPS$adj.P.Val<0.05 & limpH$adj.P.Val<0.05 & limTNF$adj.P.Val<0.05 &
                limpHTNF$adj.P.Val>0.05 & limLPSTNF$adj.P.Val>0.05 & limpHLPS$adj.P.Val>0.05 &
                limpHLPSTNF$adj.P.Val>0.05),
          #A:B
          sum(limLPS$adj.P.Val>0.05 & limpH$adj.P.Val>0.05 & limTNF$adj.P.Val>0.05 &
                limpHTNF$adj.P.Val>0.05 & limLPSTNF$adj.P.Val>0.05 & limpHLPS$adj.P.Val<0.05 &
                limpHLPSTNF$adj.P.Val>0.05),
          #A:B:C
          sum(limLPS$adj.P.Val>0.05 & limpH$adj.P.Val>0.05 & limTNF$adj.P.Val>0.05 &
                limpHTNF$adj.P.Val>0.05 & limLPSTNF$adj.P.Val>0.05 & limpHLPS$adj.P.Val>0.05 &
                limpHLPSTNF$adj.P.Val<0.05),
          #A:B + A
          sum(limLPS$adj.P.Val<0.05 & limpH$adj.P.Val>0.05 & limTNF$adj.P.Val>0.05 &
                limpHTNF$adj.P.Val>0.05 & limLPSTNF$adj.P.Val>0.05 & limpHLPS$adj.P.Val<0.05 &
                limpHLPSTNF$adj.P.Val>0.05),
          #A:B:C + A:B + A
          sum(limLPS$adj.P.Val<0.05 & limpH$adj.P.Val>0.05 & limTNF$adj.P.Val>0.05 &
                limpHTNF$adj.P.Val>0.05 & limLPSTNF$adj.P.Val>0.05 & limpHLPS$adj.P.Val<0.05 &
                limpHLPSTNF$adj.P.Val<0.05)
)

#######################DESEQ2#########################
deseqlist <- deseq2method3(sig3rep2, sampledGenes3sig)


resrpH <- deseqlist[[1]]
resrLPS <- deseqlist[[2]]
resrTNF <- deseqlist[[3]]
resrpHLPS<-deseqlist[[4]]
resrLPSTNF<- deseqlist[[5]]
resrpHTNF <- deseqlist[[6]]
resrpHLPSTNF<- deseqlist[[7]]

#vecTF <- (resrpHLPSTNF[sampledGenes3sig,]$padj<0.05)

#True Disovery
#A
DTD3 <- c(sum(resrLPS[sampledGenes3sig,]$padj<0.05 & resrpH[sampledGenes3sig,]$padj>0.05 & resrTNF[sampledGenes3sig,]$padj>0.05 &
                resrpHTNF[sampledGenes3sig,]$padj>0.05 & resrLPSTNF[sampledGenes3sig,]$padj>0.05 & resrpHLPS[sampledGenes3sig,]$padj>0.05 &
                resrpHLPSTNF[sampledGenes3sig,]$padj>0.05 &  rep(1:6, 900)==1),
          #A+B+C
          sum(resrLPS[sampledGenes3sig,]$padj<0.05 & resrpH[sampledGenes3sig,]$padj<0.05 & resrTNF[sampledGenes3sig,]$padj<0.05 &
                resrpHTNF[sampledGenes3sig,]$padj>0.05 & resrLPSTNF[sampledGenes3sig,]$padj>0.05 & resrpHLPS[sampledGenes3sig,]$padj>0.05 &
                resrpHLPSTNF[sampledGenes3sig,]$padj>0.05 &  rep(1:6, 900)==2),
          #A:B
          sum(resrLPS[sampledGenes3sig,]$padj>0.05 & resrpH[sampledGenes3sig,]$padj>0.05 & resrTNF[sampledGenes3sig,]$padj>0.05 &
                resrpHTNF[sampledGenes3sig,]$padj>0.05 & resrLPSTNF[sampledGenes3sig,]$padj>0.05 & resrpHLPS[sampledGenes3sig,]$padj<0.05 &
                resrpHLPSTNF[sampledGenes3sig,]$padj>0.05 &  rep(1:6, 900)==3),
          #A:B:C
          sum(resrLPS[sampledGenes3sig,]$padj>0.05 & resrpH[sampledGenes3sig,]$padj>0.05 & resrTNF[sampledGenes3sig,]$padj>0.05 &
                resrpHTNF[sampledGenes3sig,]$padj>0.05 & resrLPSTNF[sampledGenes3sig,]$padj>0.05 & resrpHLPS[sampledGenes3sig,]$padj>0.05 &
                resrpHLPSTNF[sampledGenes3sig,]$padj<0.05 &  rep(1:6, 900)==4),
          #A:B + A
          sum(resrLPS[sampledGenes3sig,]$padj<0.05 & resrpH[sampledGenes3sig,]$padj>0.05 & resrTNF[sampledGenes3sig,]$padj>0.05 &
                resrpHTNF[sampledGenes3sig,]$padj>0.05 & resrLPSTNF[sampledGenes3sig,]$padj>0.05 & resrpHLPS[sampledGenes3sig,]$padj<0.05 &
                resrpHLPSTNF[sampledGenes3sig,]$padj>0.05 &  rep(1:6, 900)==5),
          #A:B:C + A:B + A
          sum(resrLPS[sampledGenes3sig,]$padj<0.05 & resrpH[sampledGenes3sig,]$padj>0.05 & resrTNF[sampledGenes3sig,]$padj>0.05 &
                resrpHTNF[sampledGenes3sig,]$padj>0.05 & resrLPSTNF[sampledGenes3sig,]$padj>0.05 & resrpHLPS[sampledGenes3sig,]$padj<0.05 &
                resrpHLPSTNF[sampledGenes3sig,]$padj<0.05 &  rep(1:6, 900)==6)
)

#TrueMis Disovery
#A
DTMD3<- c(sum(resrLPS[sampledGenes3sig,]$padj<0.05 & resrpH[sampledGenes3sig,]$padj>0.05 & resrTNF[sampledGenes3sig,]$padj>0.05 &
                resrpHTNF[sampledGenes3sig,]$padj>0.05 & resrLPSTNF[sampledGenes3sig,]$padj>0.05 & resrpHLPS[sampledGenes3sig,]$padj>0.05 &
                resrpHLPSTNF[sampledGenes3sig,]$padj>0.05 ),
          #A+B+C
          sum(resrLPS[sampledGenes3sig,]$padj<0.05 & resrpH[sampledGenes3sig,]$padj<0.05 & resrTNF[sampledGenes3sig,]$padj<0.05 &
                resrpHTNF[sampledGenes3sig,]$padj>0.05 & resrLPSTNF[sampledGenes3sig,]$padj>0.05 & resrpHLPS[sampledGenes3sig,]$padj>0.05 &
                resrpHLPSTNF[sampledGenes3sig,]$padj>0.05),
          #A:B
          sum(resrLPS[sampledGenes3sig,]$padj>0.05 & resrpH[sampledGenes3sig,]$padj>0.05 & resrTNF[sampledGenes3sig,]$padj>0.05 &
                resrpHTNF[sampledGenes3sig,]$padj>0.05 & resrLPSTNF[sampledGenes3sig,]$padj>0.05 & resrpHLPS[sampledGenes3sig,]$padj<0.05 &
                resrpHLPSTNF[sampledGenes3sig,]$padj>0.05),
          #A:B:C
          sum(resrLPS[sampledGenes3sig,]$padj>0.05 & resrpH[sampledGenes3sig,]$padj>0.05 & resrTNF[sampledGenes3sig,]$padj>0.05 &
                resrpHTNF[sampledGenes3sig,]$padj>0.05 & resrLPSTNF[sampledGenes3sig,]$padj>0.05 & resrpHLPS[sampledGenes3sig,]$padj>0.05 &
                resrpHLPSTNF[sampledGenes3sig,]$padj<0.05),
          #A:B + A
          sum(resrLPS[sampledGenes3sig,]$padj<0.05 & resrpH[sampledGenes3sig,]$padj>0.05 & resrTNF[sampledGenes3sig,]$padj>0.05 &
                resrpHTNF[sampledGenes3sig,]$padj>0.05 & resrLPSTNF[sampledGenes3sig,]$padj>0.05 & resrpHLPS[sampledGenes3sig,]$padj<0.05 &
                resrpHLPSTNF[sampledGenes3sig,]$padj>0.05),
          #A:B:C + A:B + A
          sum(resrLPS[sampledGenes3sig,]$padj<0.05 & resrpH[sampledGenes3sig,]$padj>0.05 & resrTNF[sampledGenes3sig,]$padj>0.05 &
                resrpHTNF[sampledGenes3sig,]$padj>0.05 & resrLPSTNF[sampledGenes3sig,]$padj>0.05 & resrpHLPS[sampledGenes3sig,]$padj<0.05 &
                resrpHLPSTNF[sampledGenes3sig,]$padj<0.05)
)

#All Disovery
#A
DAD3 <- c(sum(resrLPS$padj<0.05 & resrpH$padj>0.05 & resrTNF$padj>0.05 &
                resrpHTNF$padj>0.05 & resrLPSTNF$padj>0.05 & resrpHLPS$padj>0.05 &
                resrpHLPSTNF$padj>0.05 ),
          #A+B+C
          sum(resrLPS$padj<0.05 & resrpH$padj<0.05 & resrTNF$padj<0.05 &
                resrpHTNF$padj>0.05 & resrLPSTNF$padj>0.05 & resrpHLPS$padj>0.05 &
                resrpHLPSTNF$padj>0.05),
          #A:B
          sum(resrLPS$padj>0.05 & resrpH$padj>0.05 & resrTNF$padj>0.05 &
                resrpHTNF$padj>0.05 & resrLPSTNF$padj>0.05 & resrpHLPS$padj<0.05 &
                resrpHLPSTNF$padj>0.05),
          #A:B:C
          sum(resrLPS$padj>0.05 & resrpH$padj>0.05 & resrTNF$padj>0.05 &
                resrpHTNF$padj>0.05 & resrLPSTNF$padj>0.05 & resrpHLPS$padj>0.05 &
                resrpHLPSTNF$padj<0.05),
          #A:B + A
          sum(resrLPS$padj<0.05 & resrpH$padj>0.05 & resrTNF$padj>0.05 &
                resrpHTNF$padj>0.05 & resrLPSTNF$padj>0.05 & resrpHLPS$padj<0.05 &
                resrpHLPSTNF$padj>0.05),
          #A:B:C + A:B + A
          sum(resrLPS$padj<0.05 & resrpH$padj>0.05 & resrTNF$padj>0.05 &
                resrpHTNF$padj>0.05 & resrLPSTNF$padj>0.05 & resrpHLPS$padj<0.05 &
                resrpHLPSTNF$padj<0.05)
)

###############################################################
#################### O W N    M E T H O D #####################
###############################################################

LinSig_modelresult <- ownfunc(sig3rep2)
GT3sigstat <- LinSig_modelresult[[1]]
dfd <- LinSig_modelresult[[2]]
#Calculate FDR for 3 signals

s3 <- as.character(sampledGenes3sig)
rt <- 0.8 # R2 threshold

# use the original random dataset and the model statistics
fdrlist <- fdrcalcfunc(rdf3_0, GT3sigstat, reps=2, dfd=dfd)

fdr1 <- fdrlist[[1]]
fdr2 <- fdrlist[[2]]
fdr3 <- fdrlist[[3]]
fdr4 <- fdrlist[[4]]
fdr5 <- fdrlist[[5]]
fdr6 <- fdrlist[[6]]
fdr7 <- fdrlist[[7]]

#True Discovery
#A
TD3 <- c(sum(GT3sigstat[s3,]$p_x1<0.05 & fdr2[s3]<0.05 & fdr1[s3]>0.05 & fdr3[s3]>0.05 & rep(1:6, 900)==1 &
               fdr4[s3]>0.05 & fdr5[s3]>0.05 & fdr6[s3]>0.05 & fdr7[s3]>0.05 & GT3sigstat[s3,]$R2>rt),
         #A + B + C
         sum(GT3sigstat[s3,]$p_x1<0.05 & GT3sigstat[s3,]$p_x2<0.05 & GT3sigstat[s3,]$p_x3<0.05 & fdr2[s3]<0.05 & fdr1[s3]<.05 & fdr3[s3]<0.05 &
               fdr4[s3]>0.05 & fdr5[s3]>0.05 & fdr6[s3]>0.05 & fdr7[s3]>0.05 & GT3sigstat[s3,]$R2>rt &rep(1:6, 900)==2),
         #A:B
         sum(GT3sigstat[s3,]$p_x4<0.05 & fdr2[s3]>0.05 & fdr1[s3]>.05 & fdr3[s3]>0.05 & rep(1:6, 900)==3 &
               fdr4[s3]<0.05 & fdr5[s3]>0.05 & fdr6[s3]>0.05 & fdr7[s3]>0.05 & GT3sigstat[s3,]$R2>rt),
         #A:B:C
         sum(GT3sigstat[s3,]$p_x7<0.05 & fdr2[s3]>0.05 & fdr1[s3]>.05 & fdr3[s3]>0.05 & rep(1:6, 900)==4 &
               fdr4[s3]>0.05 & fdr5[s3]>0.05 & fdr6[s3]>0.05 & fdr7[s3]<.05 & GT3sigstat[s3,]$R2>rt),
         #A:B + A
         sum(GT3sigstat[s3,]$p_x1<0.05 & GT3sigstat[s3,]$p_x4<0.05 & fdr2[s3]<0.05 & fdr1[s3]>.05 & fdr3[s3]>0.05 & rep(1:6, 900)==5 &
               fdr4[s3]<0.05 & fdr5[s3]>0.05 & fdr6[s3]>0.05 & fdr7[s3]>0.05 & GT3sigstat[s3,]$R2>rt),
         #A:B:C + A:B + A
         sum(GT3sigstat[s3,]$p_x1<0.05 & GT3sigstat[s3,]$p_x4<0.05 & GT3sigstat[s3,]$p_x7<0.05 & fdr2[s3]<0.05 & fdr1[s3]>.05 & fdr3[s3]>0.05 & rep(1:6, 900)==6 &
               fdr4[s3]<0.05 & fdr5[s3]>0.05 & fdr6[s3]>0.05 & fdr7[s3]<0.05 & GT3sigstat[s3,]$R2>rt)
)

#TrueMis Discovery
TMD3 <- c(sum(GT3sigstat[s3,]$p_x1<0.05 & fdr2[s3]<0.05 & fdr1[s3]>0.05 & fdr3[s3]>0.05 &
                fdr4[s3]>0.05 & fdr5[s3]>0.05 & fdr6[s3]>0.05 & fdr7[s3]>0.05 & GT3sigstat[s3,]$R2>rt),
          #A + B + C
          sum(GT3sigstat[s3,]$p_x1<0.05 & GT3sigstat[s3,]$p_x2<0.05 & GT3sigstat[s3,]$p_x3<0.05 & fdr2[s3]<0.05 & fdr1[s3]<.05 & fdr3[s3]<0.05 &
                fdr4[s3]>0.05 & fdr5[s3]>0.05 & fdr6[s3]>0.05 & fdr7[s3]>0.05 & GT3sigstat[s3,]$R2>rt),
          #A:B
          sum(GT3sigstat[s3,]$p_x4<0.05 & fdr2[s3]>0.05 & fdr1[s3]>.05 & fdr3[s3]>0.05 &
                fdr4[s3]<0.05 & fdr5[s3]>0.05 & fdr6[s3]>0.05 & fdr7[s3]>0.05 & GT3sigstat[s3,]$R2>rt),
          #A:B:C
          sum(GT3sigstat[s3,]$p_x7<0.05 & fdr2[s3]>0.05 & fdr1[s3]>.05 & fdr3[s3]>0.05 &
                fdr4[s3]>0.05 & fdr5[s3]>0.05 & fdr6[s3]>0.05 & fdr7[s3]<.05 & GT3sigstat[s3,]$R2>rt),
          #A:B + A
          sum(GT3sigstat[s3,]$p_x1<0.05 & GT3sigstat[s3,]$p_x4<0.05 & fdr2[s3]<0.05 & fdr1[s3]>.05 & fdr3[s3]>0.05 &
                fdr4[s3]<0.05 & fdr5[s3]>0.05 & fdr6[s3]>0.05 & fdr7[s3]>0.05 & GT3sigstat[s3,]$R2>rt),
          #A:B:C + A:B + A
          sum(GT3sigstat[s3,]$p_x1<0.05 & GT3sigstat[s3,]$p_x4<0.05 & GT3sigstat[s3,]$p_x7<0.05 & fdr2[s3]<0.05 & fdr1[s3]>.05 & fdr3[s3]>0.05 &
                fdr4[s3]<0.05 & fdr5[s3]>0.05 & fdr6[s3]>0.05 & fdr7[s3]<0.05 & GT3sigstat[s3,]$R2>rt)
)

#All Discovery
AD3 <- c(sum(GT3sigstat$p_x1<0.05 & fdr2<0.05 & fdr1>0.05 & fdr3>0.05 &
               fdr4>0.05 & fdr5>0.05 & fdr6>0.05 & fdr7>0.05 & GT3sigstat$R2>rt),
         #A + B + C
         sum(GT3sigstat$p_x1<0.05 & GT3sigstat$p_x2<0.05 & GT3sigstat$p_x3<0.05 & fdr2<0.05 & fdr1<.05 & fdr3<0.05 &
               fdr4>0.05 & fdr5>0.05 & fdr6>0.05 & fdr7>0.05 & GT3sigstat$R2>rt),
         #A:B
         sum(GT3sigstat$p_x4<0.05 & fdr2>0.05 & fdr1>.05 & fdr3>0.05 &
               fdr4<0.05 & fdr5>0.05 & fdr6>0.05 & fdr7>0.05 & GT3sigstat$R2>rt),
         #A:B:C
         sum(GT3sigstat$p_x7<0.05 & fdr2>0.05 & fdr1>.05 & fdr3>0.05 &
               fdr4>0.05 & fdr5>0.05 & fdr6>0.05 & fdr7<.05 & GT3sigstat$R2>rt),
         #A:B + A
         sum(GT3sigstat$p_x1<0.05 & GT3sigstat$p_x4<0.05 & fdr2<0.05 & fdr1>.05 & fdr3>0.05 &
               fdr4<0.05 & fdr5>0.05 & fdr6>0.05 & fdr7>0.05 & GT3sigstat$R2>rt),
         #A:B:C + A:B + A
         sum(GT3sigstat$p_x1<0.05 & GT3sigstat$p_x4<0.05 & GT3sigstat$p_x7<0.05 & fdr2<0.05 & fdr1>.05 & fdr3>0.05 &
               fdr4<0.05 & fdr5>0.05 & fdr6>0.05 & fdr7<0.05 & GT3sigstat$R2>rt)
)


##########
# FIGURE #
##########

FD3 = AD3 - TMD3
DFD3 = DAD3 - DTMD3
LFD3 = LAD3 - LTMD3

TrueDisc <- c(TD3, DTD3,LTD3) # 6 types 2 methods
TrueMisc <- c(TMD3-TD3, DTMD3-DTD3, LTMD3-LTD3)
FalsDisc <- c(AD3-TMD3, DAD3-DTMD3, LAD3-LTMD3)

RegRec3sig <- data.frame("rate"=c(TrueDisc, TrueMisc, FalsDisc),
                         "discovery"= rep(c("True Disc", "Mis Class", "False Disc"), each=18),
                         "method"=rep(rep(c("LinSig", "DEseq2", "limma"),each=6),3),
                         "regulation"=rep(c("A", "A+B+C", "A:B", "A:B:C", "A:B+A", "A:B:C+A:B+A"),9)
)
#cRegRec3sig <- NULL
cRegRec3sig <- rbind(cRegRec3sig, RegRec3sig)
#append dataframe to each other and run this script 10 times or so


ggplot(data=cRegRec3sig, aes(x = method, y = rate/(nrow(cRegRec3sig)/54), fill = discovery)) +
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ factor(regulation, level = c("A", "A+B+C", "A:B", "A:B:C", "A:B+A", "A:B:C+A:B+A"))) +
  theme_bw()+
  ylab("Counts")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_hline(yintercept=900, linetype="dashed", color="red")


cRegRec3sig$replicate <- rep(1:10, each=54)
write.csv(cRegRec3sig, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/cRegRec3sig.csv")
cRegRec3sig <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/cRegRec3sig.csv")

#rdf3 <- rbind(rdf30,rdf31,rdf32,rdf33,rdf34,rdf35,rdf36,rdf37,rdf38,rdf39)
#write.csv(rdf3, "rdf3_df.csv")
#rdfr3 <- rbind(rdfr30,rdfr31,rdfr32,rdfr33,rdfr34,rdfr35,rdfr36,rdfr37,rdfr38,rdfr39)
#write.csv(rdfr3, "rdfr3_df.csv")


########################################################
### Regulation Recovery based on signal size 1.5-10? ###
########################################################
#A / A:B / A:B:C

#################################
############ DESEQ2 #############
#################################

#f<-1 # 2^{1} fold change simulation
#s3mv <- c( c(0,0,0,0,f,f,0,0,f,f,0,0,f,f,f,f),       # A
#          c(0,0,0,0,0,0,0,0,f,f,0,0,0,0,f,f),       # A:B
#          c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,f,f))      # A:B:C
#sig3GTv <- matrix(rep(c(s3mv, -s3mv), 900), ncol=16, byrow=T) # 5400 GT genes

scalingFactors <- log2(rep(c(1.5,1.625,1.75,1.875,2,2.125,2.25,2.375,2.5,2.625,2.75,2.875,3,3.125,3.25),360))


fVGTs3 <- matrix(rbind(matrix(rep(c(0,0,0,0,1,1,0,0,1,1,0,0,1,1,1,1),1800), ncol=16,byrow=T),
                        matrix(rep(c(0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1),1800), ncol=16, byrow=T),
                        matrix(rep(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),1800), ncol=16, byrow=T)),ncol=16)


emptyMatrix <- matrix(0, nrow=nrow(rdf3_0), ncol=16)
emptyMatrix[as.numeric(sampledGenes3sig),] <- fVGTs3 * scalingFactors * sample(c(-1,1), 5400, replace=TRUE)

# CHANGE RANDOM DATAFRAME HERE:
rdf <- rdf3_9

sig3rep2v<-(2^(emptyMatrix)*rdf)


rownames(sig3rep2v) <-as.character(1:nrow(sig3rep2v))


############# DESEQ2 #################
ds2objs <- deseq2method3(sig3rep2v, sampledGenes3sig)

ABCds <- ds2objs[[7]][sampledGenes3sig,][3601:5400,] # A:B:C
ABds <- ds2objs[[4]][sampledGenes3sig,][1801:3600,] # A:B
Ads <- ds2objs[[2]][sampledGenes3sig,][1:1800,] # A

############ LIMMA ###################
limmaObj <- limmaMethod3(sig3rep2v, sampledGenes3sig)

ABCli <- limmaObj[[7]][sampledGenes3sig,][3601:5400,]
ABli <- limmaObj[[4]][sampledGenes3sig,][1801:3600,]
Ali <- limmaObj[[2]][sampledGenes3sig,][1:1800,]
##################################
########### OWN METHOD ###########
##################################

VGTstat3 <- ownfunc(sig3rep2v)
dfd <- VGTstat3[[2]]
VGTstat3 <- VGTstat3[[1]]

#VGTstat3[as.character(sampledGenes3sigv),]

fdrlist <- fdrcalcfunc(rdf, VGTstat3, reps=2, dfd=dfd)
sum(fdrlist[[2]][s3]<0.05)
sum(fdrlist[[4]][s3]<0.05)
sum(fdrlist[[7]][s3]<0.05)

ABCom<- VGTstat3[as.character(sampledGenes3sig),][3601:5400,]
ABCom$ABCFDR <- fdrlist[[7]][s3][3601:5400]
ABom <- VGTstat3[as.character(sampledGenes3sig),][1801:3600,]
ABom$ABFDR <- fdrlist[[4]][s3][1801:3600]
Aom  <- VGTstat3[as.character(sampledGenes3sig),][1:1800,]
Aom$AFDR <- fdrlist[[2]][s3][1:1800]

SR <- c()
#DESEQ2
for (i in 1:15){
  SR <- append(SR, (sum(ABCds[15*0:119+i,]$padj<0.05))/120) # A:B:C
}
for (i in 1:15){
  SR <- append(SR, (sum(ABds[15*0:119+i,]$padj<0.05)/120)) # A:B
}
for (i in 1:15){
  SR <- append(SR, (sum(Ads[15*0:119+i,]$padj<0.05)/120)) # A:B
}

# LIMMA
for (i in 1:15){
  SR <- append(SR, (sum(ABCli[15*0:119+i,]$adj.P.Val<0.05))/120) # A:B:C
}
for (i in 1:15){
  SR <- append(SR, (sum(ABli[15*0:119+i,]$adj.P.Val<0.05)/120)) # A:B
}
for (i in 1:15){
  SR <- append(SR, (sum(Ali[15*0:119+i,]$adj.P.Val<0.05)/120)) # A:B
}

# OWN METHOD
for (i in 1:15){
  SR <- append(SR, sum(ABCom[15*0:119+i,]$p_x7<0.05 & ABCom[15*0:119+i,]$R2>0.8 & ABCom[15*0:119+i,]$ABCFDR<0.05)/120)
}
for (i in 1:15){
  SR <- append(SR, sum(ABom[15*0:119+i,]$p_x7<0.05 & ABom[15*0:119+i,]$R2>0.8 & ABom[15*0:119+i,]$ABFDR<0.05)/120)
}
for (i in 1:15){
  SR <- append(SR, sum(Aom[15*0:119+i,]$p_x7<0.05 & Aom[15*0:119+i,]$R2>0.8 & Aom[15*0:119+i,]$AFDR<0.05)/120)
}


sensDF <- data.frame("rate"=SR,
                     "method"=rep(c("DESeq2", "limma", "LinSig"), each=45),
                     "term" = rep(rep(c("3int","2int", "main"), each=15),3),
                     "beta" = (rep(scalingFactors[1:15], 9))
)
#fullsensDF<-NULL
fullsensDF <- rbind(fullsensDF, sensDF)

# ggplot(data=fullsensDF[fullsensDF$beta>0,], aes(x=beta, y=rate, group=interaction(method, term), color=term, linetype=method))+
#   geom_smooth(method="loess",se=T, span=0.7, aes(fill=term))+
#   #geom_line(linewidth=1, stat="summary", fun=mean)+
#   #geom_point(aes(shape=term))+
#   #ylim(0.1,1)
#   ylab("Sensitivity")+
#   xlab("fold change")+
#   coord_cartesian(ylim=c(0.4, 1))

ggplot(data=fullsensDF[fullsensDF$beta>0,], aes(x=beta, y=rate, group=interaction(method, term), color=method, linetype=term))+
  geom_smooth(method="loess",se=T, span=0.7, aes(fill=method), alpha=0.3)+
  scale_linetype_manual(values=c(1,3,5)) + 
  #geom_line(linewidth=1, stat="summary", fun=mean)+
  #geom_point(aes(shape=term))+
  ylab("Sensitivity")+
  xlab("fold change")+
  coord_cartesian(ylim=c(0.4, 1))

write.csv(fullsensDF, "C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/SEN_10reps_3sig.csv")
fullsensDF <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/LinSigPaper/3sig_data/SEN_10reps_3sig.csv")
#calc FDR for FDR vs SENS plot


avgFDR <- cRegRec3sig[(cRegRec3sig$regulation=="A" | cRegRec3sig$regulation=="A:B" |cRegRec3sig$regulation=="A:B:C") & cRegRec3sig$discovery=="False Disc",] %>%
  group_by(regulation, method, discovery) %>%
  summarise(FDR_count=mean(rate))
avgFDR

avgTRD <- cRegRec3sig[(cRegRec3sig$regulation=="A" | cRegRec3sig$regulation=="A:B" |cRegRec3sig$regulation=="A:B:C") & cRegRec3sig$discovery=="True Disc",] %>%
  group_by(regulation, method, discovery) %>%
  summarise(True_count=mean(rate))
avgTRD

avgSENS <- cRegRec3sig[(cRegRec3sig$regulation=="A" | cRegRec3sig$regulation=="A:B" |cRegRec3sig$regulation=="A:B:C") & cRegRec3sig$discovery!="False Disc",] %>%
  group_by(regulation, method, discovery) %>%
  summarise(mean_rate=mean(rate)) %>%
  group_by(method, regulation) %>%
  summarise(total_count=sum(mean_rate))
avgSENS


combstatFDRSENS <- merge(merge(avgFDR, avgSENS, by=c("method", "regulation")), avgTRD, by=c("method", "regulation"))

FDRs3sig <- cRegRec3sig[(cRegRec3sig$regulation=="A" | cRegRec3sig$regulation=="A:B" |cRegRec3sig$regulation=="A:B:C") & cRegRec3sig$discovery=="False Disc",]
SENs3sig <- cRegRec3sig[(cRegRec3sig$regulation=="A" | cRegRec3sig$regulation=="A:B" |cRegRec3sig$regulation=="A:B:C") & cRegRec3sig$discovery=="True Disc",]

SENFDRdf_3sig <- SENs3sig
SENFDRdf_3sig$fdr_rate <- (FDRs3sig$rate/SENs3sig$rate)
SENFDRdf_3sig$sen_rate <- (SENFDRdf_3sig$rate/900)


ggplot(data=SENFDRdf_3sig, aes(x=sen_rate, y=1-fdr_rate, color=method)) +
  geom_point(aes(shape=regulation), size=2)+
  coord_cartesian(ylim=c(0.80,1), xlim=c(0.45,1))+
  xlab("Sensitivity") +
  ylab("1-FDR")+
  stat_ellipse(geom="polygon", level=0.85, aes(fill=method), alpha=0.25)

ggplot(data=SENFDRdf_3sig, aes(x=sen_rate, y=1-fdr_rate, fill=method)) +
  geom_point(aes(shape=regulation), size=2)+
  scale_shape_manual(values=c(21, 22, 24)) + 
  coord_cartesian(ylim=c(0.80,1), xlim=c(0.45,1))+
  xlab("Sensitivity") +
  ylab("1-FDR")+
  stat_ellipse(geom="polygon", level=0.85, aes(fill=method), alpha=0.25)
