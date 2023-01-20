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
#'
#' @description The function will ask you to enter column numbers per condition of the experiment
#' If you have multiple replicates, enter the column numbers with a comma inbetween.
#' For example: control condition is in columns 1 and 2. Enter: 1,2
#' The output of the function is a dataframe with LogRatios to be used for the fitting of the model

compute_ratios <- function(data, pseudo=1, lowess_norm=1){
  message("You will be asked to enter column names for each condition.")
  message("If you have replicates, separate the column numbers with a comma. Example: 1,2")

  ctrl  <- readline("Enter column number(s) for condition ctrl:")
  condA <- readline("Enter column number(s) for condition A:")
  condB <- readline("Enter column number(s) for condition B:")
  condAB<- readline("Enter column number(s) for condition AB:")

  .pkglobenv$col_ctrl <- as.numeric(unlist(strsplit(ctrl, ",")))
  .pkglobenv$col_condA <- as.numeric(unlist(strsplit(condA, ",")))
  .pkglobenv$col_condB <- as.numeric(unlist(strsplit(condB, ",")))
  .pkglobenv$col_condAB <- as.numeric(unlist(strsplit(condAB, ",")))

  DataPseudo <- data + pseudo
  .pkglobenv$pseudo_count <- pseudo
  .pkglobenv$DataPseudo_obj <- DataPseudo
  reps <- length(.pkglobenv$col_ctrl)

  LogRatios <- matrix(nrow=length(DataPseudo[,1]), ncol=reps*5) #empty matrix
  LogIntensity <- matrix(nrow=length(DataPseudo[,1]), ncol=reps*5) #empty matrix

  for (i in 1:reps){
    LogRatios[,1+ (5*(i-1))] <- log2(DataPseudo[, .pkglobenv$col_condA[i]] / DataPseudo[, .pkglobenv$col_ctrl[i]])
    LogRatios[,2+ (5*(i-1))] <- log2(DataPseudo[, .pkglobenv$col_condAB[i]] / DataPseudo[, .pkglobenv$col_condB[i]])
    LogRatios[,3+ (5*(i-1))] <- log2(DataPseudo[, .pkglobenv$col_condB[i]] / DataPseudo[, .pkglobenv$col_ctrl[i]])
    LogRatios[,4+ (5*(i-1))] <- log2(DataPseudo[, .pkglobenv$col_condAB[i]] / DataPseudo[, .pkglobenv$col_condA[i]])
    LogRatios[,5+ (5*(i-1))] <- log2(DataPseudo[, .pkglobenv$col_condAB[i]] / DataPseudo[, .pkglobenv$col_ctrl[i]])

    LogIntensity[,1+ (5*(i-1))] <- log2(DataPseudo[, .pkglobenv$col_condA[i]] * DataPseudo[, .pkglobenv$col_ctrl[i]])
    LogIntensity[,2+ (5*(i-1))] <- log2(DataPseudo[, .pkglobenv$col_condAB[i]] * DataPseudo[, .pkglobenv$col_condB[i]])
    LogIntensity[,3+ (5*(i-1))] <- log2(DataPseudo[, .pkglobenv$col_condB[i]] * DataPseudo[, .pkglobenv$col_ctrl[i]])
    LogIntensity[,4+ (5*(i-1))] <- log2(DataPseudo[, .pkglobenv$col_condAB[i]] * DataPseudo[, .pkglobenv$col_condA[i]])
    LogIntensity[,5+ (5*(i-1))] <- log2(DataPseudo[, .pkglobenv$col_condAB[i]] * DataPseudo[, .pkglobenv$col_ctrl[i]])
  }
  n_ratios <- length(LogRatios[1,])

  if (lowess_norm==0) {
    message("no lowess normalization")
    return(LogRatios)}

  else{
  message("performing lowess normalization")
  Ratios <- matrix(ncol=n_ratios, nrow=length(LogRatios[,1]))
    pb <- txtProgressBar(min = 0, max = n_ratios, initial = 0, char = "=", style = 3)

    for (i in 1:n_ratios){
      Ratios[,i] <- RNAseqLowess(LogIntensity[,i], LogRatios[,i])
      setTxtProgressBar(pb,i)
    }
    return(Ratios)
  }
}
