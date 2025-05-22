compute_ratios <- function(df, pseudo_count=1, lowess_norm=FALSE, structureDataFrame=NULL,
                           n_replicates=2){
  
  
  if (typeof(structureDataFrame) == "list") {
    col_ctrl <- structureDataFrame$col_ctrl
    col_condA <- structureDataFrame$col_condA
    col_condB <- structureDataFrame$col_condB
    col_condAB<- structureDataFrame$col_condAB
    n_replicates=length(col_ctrl)
  }
  else {
    # else assume column order based on example dataset (n_replicates=2)
    col_ctrl <- 1:n_replicates
    col_condA <- 1:n_replicates + n_replicates
    col_condB <- 1:n_replicates + n_replicates*2
    col_condAB <- 1:n_replicates + n_replicates*3
  }
  
  dfPseudo <- df + pseudo_count # add pseudo count to dataset
  
  
  LogRatios <- matrix(nrow=nrow(dfPseudo), ncol=n_replicates*5) # initialize empty matrix
  LogIntensity <- matrix(nrow=nrow(dfPseudo), ncol=n_replicates*5) # initialize empty matrix
  
  
  # Compute log ratios
  print("computing log ratios")
  for (i in 1:n_replicates){
    LogRatios[,1+ (5*(i-1))] <- log2(dfPseudo[, col_condA[i]] / dfPseudo[, col_ctrl[i]])
    LogRatios[,2+ (5*(i-1))] <- log2(dfPseudo[, col_condAB[i]] / dfPseudo[, col_condB[i]])
    LogRatios[,3+ (5*(i-1))] <- log2(dfPseudo[, col_condB[i]] / dfPseudo[, col_ctrl[i]])
    LogRatios[,4+ (5*(i-1))] <- log2(dfPseudo[, col_condAB[i]] / dfPseudo[, col_condA[i]])
    LogRatios[,5+ (5*(i-1))] <- log2(dfPseudo[, col_condAB[i]] / dfPseudo[, col_ctrl[i]])
    
    LogIntensity[,1+ (5*(i-1))] <- log2(dfPseudo[, col_condA[i]] * dfPseudo[, col_ctrl[i]])
    LogIntensity[,2+ (5*(i-1))] <- log2(dfPseudo[, col_condAB[i]] * dfPseudo[, col_condB[i]])
    LogIntensity[,3+ (5*(i-1))] <- log2(dfPseudo[, col_condB[i]] * dfPseudo[, col_ctrl[i]])
    LogIntensity[,4+ (5*(i-1))] <- log2(dfPseudo[, col_condAB[i]] * dfPseudo[, col_condA[i]])
    LogIntensity[,5+ (5*(i-1))] <- log2(dfPseudo[, col_condAB[i]] * dfPseudo[, col_ctrl[i]])
  }
  n_ratios <- nrow(LogRatios)
  
  
  if (!lowess_norm) { # skip ahead if no Lowess normalization
    rownames(LogRatios) <- rownames(df)
    return(LogRatios)
  }
  
  print("Performing LOWESS normalization...")
  normalizedRatios <- matrix(nrow = nrow(LogRatios), ncol = ncol(LogRatios))
  #colnames(Ratios) <- col_names
  pb <- txtProgressBar(min = 0, max = ncol(LogRatios), style = 3)
  
  for (i in 1:ncol(LogRatios)) {
    normalizedRatios[, i] <- RNAseqLowess(LogIntensity[, i], LogRatios[, i])
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  rownames(normalizedRatios) <- rownames(df)
  return(normalizedRatios)
}



###################################
######## F U N C T I O N S ########
###################################


# Function for normalization with LOWESS
normalize <- function(inputdf, countThres, pseudo, replicates, lowess=FALSE){
  
  print("Applying Normalisation")
  cntMat <- as.matrix(inputdf[,1:8])
  DataPseudo <- MedianNorm(cntMat, countThres=countThres, pseudo=pseudo)
  
  Ratios <- compute_ratios(df=DataPseudo,
                           pseudo_count = pseudo,
                           n_replicates=replicates,
                           lowess_norm = lowess)
  
  return(Ratios)
  
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