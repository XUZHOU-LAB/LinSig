
#' Fits ratios to a linear model
#'
#' `ratios.fit` ratios are fitted to a linear model
#' @param ratios input ratios data
#' @param CompThreshold ratio threshold for significance
#' @param n_rep number of replicates in data


ratios.fit <- function(ratios, CompThreshold=1.5, n_rep=2){
  strucvec <- rep(c(0,1,0,0,1,1,1,0,0,1,0,1,1,1,1), (length(ratios[1,])/5))
  X <- matrix(strucvec, ncol=3, byrow=T)

  multiplefit <- lm(t(ratios)~X)
  msums <- summary(multiplefit)
  B <- t(multiplefit$coefficients)
  resid <- t(multiplefit$residuals)




  rsquared <- vector(length = length(msums))
  ngen <- length(msums)
  covB_l <- NULL
  covB <- data.frame(cov_int=double(ngen), cov_x1=double(ngen),
                     cov_x2=double(ngen), cov_x3=double(ngen))

  message("computing covariance of coefficients...")
  pb <- txtProgressBar(min = 0, max = length(msums), initial = 0, char = "=", style = 3)
  for (i in 1:length(msums)){
    rsquared[i] <- msums[[i]]$r.squared
    covB[i,] <- diag(vcov(lm(ratios[i,]~X)))
    #covB[i,] <- msums[[i]]$cov.unscaled*(msums[[i]]$sigma**2) #inaccurate
    setTxtProgressBar(pb,i)
  }
  n_samples <- 4*n_rep
  n_var <- 3
  DoF <- n_samples - n_var - 1
  ttest_stat <- (abs(B) - log2(CompThreshold)) / sqrt(covB)
  ttest_stat <- data.frame(ttest_stat)

  Pvalue = 1 - apply(ttest_stat, 2, pt, df=DoF)

  Datafit = B %*% t(cbind(rep(1,10), X)) #matrix multiplication
  Residual = ratios - Datafit
  Expression_variation = rowMeans(ratios**2)
  Expression_residual = rowMeans(Residual**2)
  Varexplain = 100*(Expression_variation - Expression_residual) / Expression_variation

  Fit <- data.frame(Expression_variation, Varexplain)

  #DataPseudo <- global_DataPseudo_obj

  #F_condA_ctrl <- log2(rowMeans(DataPseudo[,col_condA]) / rowMeans(DataPseudo[,col_ctrl]))
  #F_condAB_condB <- log2(rowMeans(DataPseudo[,col_condAB]) / rowMeans(DataPseudo[,col_condB]))
  #F_condAB_condA <- log2(rowMeans(DataPseudo[,col_condAB]) / rowMeans(DataPseudo[,col_condA]))
  #F_condB_ctrl <- log2(rowMeans(DataPseudo[,col_condB]) / rowMeans(DataPseudo[,col_ctrl]))

  #Fold <- data.frame(F_condA_ctrl, F_condAB_condB, F_condAB_condA, F_condB_ctrl)

  plot_variance <- plot(Expression_variation,Varexplain,
                        main=paste("Variance Explained by fitting - pseudo:",pseudo),
                        xlab="Variance of Measurement",
                        ylab="Variance explained")


  #substats <- c('q_Ph0hLPS', 'q_Ph5hLPS', 'q_LPS5hvs0h', 'q_Ph9h')

  #tt = abs(Fold[,1:4]) >= 1 & subStatData <= 0.05
  #SigIDX <- rowSums(tt) > 0

  #plot_sig_variance <- plot(Expression_variation, Varexplain,
  #     main=paste("Variance Explained by fitting - pseudo:",Pseudo),
  #     xlab="Variance of Measurement",
  #     ylab="Variance explained")

  #points(Expression_variation[SigIDX], Varexplain[SigIDX], col="red")



  return(list(B, Pvalue, Fit, plot_variance))
}
