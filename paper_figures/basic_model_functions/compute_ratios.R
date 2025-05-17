#utils::globalVariables(c("col_ctrl", "col_condA", "col_condB", "col_condAB", "global_DataPseudo_obj", "pseudo"))

#' LOWESS normalization of LogRatios
#'
#' `RNAseqLowess` returns the lowess normalized logratios
#' @description The function uses the loess function from R with the following parameters. span=0.05, degree=1, family=gaussian and iter=1
#' @param LogIntensity Uses the LogIntensity computed by the compute_ratios function
#' @param LogRatios Uses the LogRatios computed by the compute_ratios function

RNAseqLowess <- function(LogIntensity, LogRatios){

  Ynorm = loess(LogRatios~LogIntensity,
                span=0.05,degree=1,family="gaussian",
                iterations=1,surface="direct")
  Ratiosnorm = LogRatios - fitted(Ynorm)
  return(Ratiosnorm)
}

#' Compute log2 ratios of 2 signals
#'
#' `compute_ratios` returns the log2 ratios with option of lowess normalization
#' @param data Data should be formatted with every row being a gene and every column being an experimental condition
#' @param pseudo Add a pseudo count in order to deal with zero counts. Default is 1
#' @param lowess_norm Uses lowess normalization in order to improve data quality. It is performed by default
#' @param dataStructureVector if left empty, the user will be asked to provide column numbers in the console. The user can also provide a dataframe with 4 column: col_ctrl, col_condA, col_condB, col_condAB. These columns should contain a vector with the column numbers as integers.
#'
#' @description The function will ask you to enter column numbers per condition of the experiment
#' If you have multiple replicates, enter the column numbers with a comma inbetween.
#' For example: control condition is in columns 1 and 2. Enter: 1,2
#' The output of the function is a dataframe with LogRatios to be used for the fitting of the model

compute_ratios <- function(df, pseudo_count=1, lowess_norm=FALSE, structureDataFrame=NULL){


  if (typeof(structureDataFrame) == "list") {
    col_ctrl <- structureDataFrame$col_ctrl
    col_condA <- structureDataFrame$col_condA
    col_condB <- structureDataFrame$col_condB
    col_condAB<- structureDataFrame$col_condAB
  }
  else {
    message("You will be asked to enter column names for each condition.")
    message("If you have replicates, separate the column numbers with a comma. Example: 1,2")

    ctrl  <- readline("Enter column number(s) for condition ctrl:")
    condA <- readline("Enter column number(s) for condition A:")
    condB <- readline("Enter column number(s) for condition B:")
    condAB<- readline("Enter column number(s) for condition AB:")

    col_ctrl <- as.numeric(unlist(strsplit(ctrl, ",")))
    col_condA <- as.numeric(unlist(strsplit(condA, ",")))
    col_condB <- as.numeric(unlist(strsplit(condB, ",")))
    col_condAB <- as.numeric(unlist(strsplit(condAB, ",")))
  }

  dfPseudo <- df + pseudo_count # add pseudo count to dataset

  nreps=length(col_ctrl)
  LogRatios <- matrix(nrow=nrow(dfPseudo), ncol=nreps*5) # initialize empty matrix
  LogIntensity <- matrix(nrow=nrow(dfPseudo), ncol=nreps*5) # initialize empty matrix


  # Compute log ratios
  print("computing log ratios")
  for (i in 1:nreps){
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
  Ratios <- matrix(nrow = nrow(LogRatios), ncol = ncol(LogRatios))
  #colnames(Ratios) <- col_names
  pb <- txtProgressBar(min = 0, max = ncol(LogRatios), style = 3)

  for (i in 1:ncol(LogRatios)) {
    Ratios[, i] <- RNAseqLowess(LogIntensity[, i], LogRatios[, i])
    setTxtProgressBar(pb, i)
  }
  close(pb)

  rownames(Ratios) <- rownames(df)
  return(Ratios)
}
