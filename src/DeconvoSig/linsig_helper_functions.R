geomean <- function(x){
  exp(mean(log(x)))
}

MedianNorm <- function(data, countThres=10, pseudo=1){
  
  subdata = data[rowMeans(data) >= countThres,]     # Subset rows with average count higher than 10
  subdata = Subdata + pseudo                       # Add Pseudo count of 1
  
  t = apply(subdata, 1, geomean)         # Calculate Geometric mean
  # Divide expression of every element per row to the geometric mean of that row
  # Every 8 elements of Data need to be divided by an element of t
  
  dataRatio = log2(subdata / t)
  T_med = apply(dataRatio, 2, median) # median per column of array
  T_med = 2^T_med
  C <- t(t(subdata) / T_med)
  return(C)
}

RNAseqLowess <- function(LogIntensity, LogRatios){
  
  Ynorm = loess(LogRatios~LogIntensity,
                span=0.05,degree=1,family="gaussian",
                iterations=1,surface="direct")
  Ratiosnorm = LogRatios - fitted(Ynorm)
  return(Ratiosnorm)
}




###################################
######## F U N C T I O N S ########
###################################


# Function for normalization with LOWESS
normalize <- function(inputdf, cntThres, pseudo, replicates, lowess=FALSE){
  
  print("Applying Normalisation")
  cntMat <- as.matrix(inputdf[,1:8])
  DataPseudo <- MedianNorm(cntMat, CountThres=cntThres, pseudo=pseudo)
  
  #automatically detect replicates?
  #replicate input number
  
  col_ctrl <- 1:replicates
  col_condA <- 1:replicates + replicates
  col_condB <- 1:replicates + replicates*2
  col_condAB <- 1:replicates + replicates*3
  
  LogRatios <- matrix(nrow=length(DataPseudo[,1]), ncol=replicates*5)    #empty matrix
  for (i in 1:replicates){
    LogRatios[,1+ (5*(i-1))] <- log2(DataPseudo[, col_condA[i]] / DataPseudo[, col_ctrl[i]])
    LogRatios[,2+ (5*(i-1))] <- log2(DataPseudo[, col_condAB[i]] / DataPseudo[, col_condB[i]])
    LogRatios[,3+ (5*(i-1))] <- log2(DataPseudo[, col_condB[i]] / DataPseudo[, col_ctrl[i]])
    LogRatios[,4+ (5*(i-1))] <- log2(DataPseudo[, col_condAB[i]] / DataPseudo[, col_condA[i]])
    LogRatios[,5+ (5*(i-1))] <- log2(DataPseudo[, col_condAB[i]] / DataPseudo[, col_ctrl[i]])
  }
  
  if (lowess == FALSE){
    rownames(LogRatios) <- rownames(DataPseudo)
    return(LogRatios)
  }
  
  else if (lowess == TRUE){
    LogIntensity <- matrix(nrow=length(DataPseudo[,1]), ncol=10) #empty matrix
    LogIntensity[,1] <- log2(DataPseudo[,3] * DataPseudo[,1])
    LogIntensity[,2] <- log2(DataPseudo[,7] * DataPseudo[,5]) # LOWESS
    LogIntensity[,3] <- log2(DataPseudo[,5] * DataPseudo[,1])
    LogIntensity[,4] <- log2(DataPseudo[,7] * DataPseudo[,3])
    LogIntensity[,5] <- log2(DataPseudo[,7] * DataPseudo[,1])
    LogIntensity[,6] <- log2(DataPseudo[,4] * DataPseudo[,2])
    LogIntensity[,7] <- log2(DataPseudo[,8] * DataPseudo[,6])
    LogIntensity[,8] <- log2(DataPseudo[,6] * DataPseudo[,2])
    LogIntensity[,9] <- log2(DataPseudo[,8] * DataPseudo[,4])
    LogIntensity[,10] <-log2(DataPseudo[,8] * DataPseudo[,2])
    NumRatios <- length(LogRatios[1,])
    withProgress(message="Computing Normalized Ratios", value=0,{
      Ratios <- matrix(ncol=10, nrow=length(LogRatios[,1]))
      for (i in 1:NumRatios){
        Ratios[,i] <- RNAseqLowess(LogIntensity[,i], LogRatios[,i])
        incProgress(1 / NumRatios)
      }
    })
    rownames(Ratios) <- rownames(DataPseudo)
    return(Ratios)
  }
}

deconvoluteFunction <- function(ratiosDF, countDF){
  Ratios <- ratiosDF
  
  strucvec <- rep(c(0,1,0, # creates a matrix for X in y ~ X. X is dependent on number of replicates
                    0,1,1, # 5 rows in X for each replicate
                    1,0,0,
                    1,0,1,
                    1,1,1),
                  (length(Ratios[1,])/5))
  
  numberOfGenes = length(Ratios[,1])
  
  X <- matrix(strucvec, 
              ncol=3, 
              byrow=T)
  B = matrix(0, nrow=numberOfGenes,ncol=4);
  
  #get column data
  cc <- strsplit(colnames(countDF)[1:8], "_") # change to 4*n_replicates
  print(cc)
  cols <- unique(unlist(cc)[2*(1:length(cc))-1])[-1] # change 2 -> n_replicates
  print(cols)
  colnames(X) <- c(cols[2], cols[1], cols[3])
  
  rsquare <- vector()
  ngen <- length(Ratios[,1])
  covB_l <- NULL
  covB <- data.frame(cov_int=double(ngen), cov_x1=double(ngen), 
                     cov_x2=double(ngen), cov_x3=double(ngen))
  
  #adds progress bar
  withProgress(message="Computing Model Statistics", value=0,{
    
    print(paste("Total number of genes:", ngen))
    for (i in 1:ngen){
      fit <- lm(Ratios[i,] ~ X) # for each gene, fit linear model
      covB[i,] <- diag(vcov(fit))
      rsquare[i] <- summary(fit)$adj.r.squared
      incProgress(1/ngen)
    }
  })
  
  multiplefit <- lm(t(Ratios) ~ X)
  B <- t(multiplefit$coefficients)
  Resid <- multiplefit$residuals
  
  modelStats <- as.data.frame(cbind(B, covB))
  
  colnames(modelStats) <- c("int", cols[2], cols[1], cols[3], "p_int",
                            paste0("cov_", cols[2]), paste0("cov_", cols[1]), paste0("cov_", cols[3]))
  modelStats$R2 <- rsquare
  modelStats <- round(modelStats, digits=7)
  return(modelStats)
}

#compute P value with threshold
calcPvalue <- function(betas, covB, treshold){ # edit P statistic for mulitple reps
  ttest_stat <- (abs(betas) - log2(treshold)) / sqrt(covB)
  ttest_stat <- data.frame(ttest_stat)
  n_samples = 4*2 # replicates * conditions
  n_var = 3 # number of variables B1, B2, B3
  DoF <- n_samples - n_var - 1
  Pvalue = 1 - apply(ttest_stat, 2, pt, df=DoF)
  return(Pvalue)
}

NearestNeighbours <- function(x,y,i){ #x:betas, y:FDR, i:numbertotest
  loc <- which.min(abs(x-abs(i)))
  return(y[loc])
}