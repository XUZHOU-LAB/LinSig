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

deconvoluteFunction <- function(countDF, countThres,
                                n_rep, H0_threshold,
                                pseudo, lowess=F){
  
  ratios <- normalize(countDF, countThres = countThres,
                      pseudo=pseudo, replicates=n_rep,
                      lowess=lowess)
  
  modelStats <- ratios.fit(ratios=ratios, 
                          CompThreshold=H0_threshold,
                          n_rep=n_rep)
  
  colnames(modelStats) <- get_colnames(countDF)
  return(round(modelStats, 5))
}


get_colnames <- function(df){
  #get column data
  cc <- strsplit(colnames(df)[1:8], "_") # change to 4*n_replicates
  cols <- unique(unlist(cc)[2*(1:length(cc))-1])[-1] # change 2 -> n_replicates
  column_names <- c("int", cols[2], cols[1], cols[3], "p_int",
                            paste0("p_", cols[2]), paste0("p_", cols[1]), paste0("p_", cols[3]),
                            "R2")
  return(column_names)
}


NearestNeighbours <- function(x,y,i){ #x:betas, y:FDR, i:numbertotest
  loc <- which.min(abs(x-abs(i)))
  return(y[loc])
}

# Plotting function
plotLFC_R2 <- function(dat, x, y, sigs, xlab) {
  ggplot(dat, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point() +
    geom_point(data = dat[sigs, ], color = "blue") +
    coord_cartesian(xlim = c(-5, 5)) +
    labs(x = xlab, y = y)
}



ratios.fit <- function(ratios, CompThreshold=1.5, n_rep=2){
  
  message("fitting model...")
  # Generate structural vector and design matrix
  n_samples <- n_rep*5 # TODO: n_samples defined twice?
  strucvec <- rep(c(0,1,0, 
                    0,1,1, 
                    1,0,0, 
                    1,0,1, 
                    1,1,1), times = n_rep)
  X <- matrix(strucvec, ncol = 3, byrow = TRUE)
  
  # Fit multivariate linear model
  multiplefit <- lm(t(ratios) ~ X)
  residuals <- t(multiplefit$residuals)
  
  # Precompute model components for covariance calculation
  model_matrix <- cbind(1, X)  # Include intercept
  xtx_inv <- solve(crossprod(model_matrix))
  diag_xtx_inv <- diag(xtx_inv)
  n_coef <- ncol(model_matrix)
  
  # Calculate residual sum of squares and variance estimates
  rss <- rowSums(residuals^2)
  sigma_sq <- rss / (n_samples - n_coef)
  
  # Compute coefficient covariance matrices
  covB <- outer(sigma_sq, diag_xtx_inv, "*")
  colnames(covB) <- c("cov_int", "cov_x1", "cov_x2", "cov_x3")
  
  # Calculate adjusted R-squared efficiently
  row_means <- rowMeans(ratios)
  tss <- rowSums(ratios^2) - n_samples * row_means^2
  rsquared <- 1 - (rss / (n_samples - n_coef)) / (tss / (n_samples - 1))
  
  # Extract coefficients and format results
  coefficients <- t(multiplefit$coefficients)
  covB <- as.data.frame(covB)
  
  n_samples <- 4*n_rep
  n_var <- 3
  DoF <- n_samples - n_var - 1
  ttest_stat <- (abs(coefficients) - log2(CompThreshold)) / sqrt(covB) #CompThreshold is H0 hypothesis
  ttest_stat <- data.frame(ttest_stat)
  
  Pvalue = 1 - apply(ttest_stat, 2, pt, df=DoF)

  outputdf <- data.frame(cbind(coefficients,Pvalue,rsquared))
  colnames(outputdf)[9] <- "R2"
  return(outputdf)
}


## Heatmap functions

# Generate cluster labels (encoded integers)
get_cluster_labels <- function(df, beta_cols, qval_cols) {
  sts <- sign(df[, beta_cols]) * (df[, qval_cols] < 0.05)
  sts[sts == -1] <- 2
  rowSums(t(t(sts) * c(1, 3, 9)))  # 1:LPS, 3:pH, 9:pHLPS
}

# Get readable cluster codes
get_cluster_code_mapping <- function(modelTerms) {
  cluster_order <- c(1,2,3,6,21,22,19,15,17,11,9,12,10,13,18,
                     24,20,26,4,5,7,8,14,23,16,25,0)
  clusterCode <- c("T2+", "T2-", "T1+", "T1-", "T1+T3-", "T1+T2+T3-", "T2+T3-",
                   "T1-T3+", "T1-T2-T3+", "T2-T3+", "T3+", "T1+T3+", "T2+T3+",
                   "T1+T2+T3+", "T3-", "T1-T3-", "T2-T3-", "T1-T2-T3-", "T1+T2+",
                   "T1+T2-", "T1-T2+", "T1-T2-", "T1+T2-T3+", "T1+T2-T3-",
                   "T1-T2+T3+", "T1-T2+T3-", "no regulation")
  clusterCode <- gsub("T1", modelTerms[2], clusterCode)
  clusterCode <- gsub("T2", modelTerms[1], clusterCode)
  clusterCode <- gsub("T3", modelTerms[3], clusterCode)
  data.frame(clusterID = cluster_order, clusterCode = clusterCode)
}

assign_genes <- function(df, B_thr = 0.585, R_thr = 0.8) {
  minimumGenesInClus <- 20
  p <- df[, 5:8]
  
  beta_cols <- 2:4
  qval_cols <- 6:8
  
  SIG_genes <- ((abs(df[,2]) > B_thr & p[,2] < 0.05) |
                  (abs(df[,3]) > B_thr & p[,3] < 0.05) |
                  (abs(df[,4]) > B_thr & p[,4] < 0.05)) & df$R2 > R_thr
  
  df[, qval_cols] <- p[,2:4]
  sigs <- df[SIG_genes, ]
  
  clus <- get_cluster_labels(sigs, beta_cols, qval_cols)
  included <- names(which(table(clus) > minimumGenesInClus))
  
  hmdf <- sigs[clus %in% included, beta_cols]
  modelTerms <- colnames(df)[beta_cols]
  clusterNames <- get_cluster_code_mapping(modelTerms)
  
  c_order <- clusterNames$clusterID
  cc_order <- c_order[c_order %in% included]
  
  clus_split <- factor(clus[clus %in% included], levels = cc_order)
  geneClusters <- clusterNames$clusterCode[match(clus_split, clusterNames$clusterID)]
  geneClusters_ord <- factor(geneClusters, levels = clusterNames$clusterCode)
  
  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  Heatmap(as.matrix(hmdf), 
          split = geneClusters_ord, 
          col = col_fun,
          cluster_row_slices = FALSE,
          cluster_columns = FALSE,
          show_row_dend = FALSE,
          heatmap_legend_param = list(title = "LFC"))
}


