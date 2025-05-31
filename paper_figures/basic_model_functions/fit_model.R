#' Fits ratios to a linear model
#'
#' `ratios.fit` ratios are fitted to a linear model
#' @param ratios input ratios data
#' @param CompThreshold ratio threshold for significance
#' @param n_rep number of replicates in data

ratios.fit <- function(ratios, CompThreshold=1.5, n_rep=2){
  
  message("fitting model...")
  # Generate structural vector and design matrix
  n_samples <- ncol(ratios) # TODO: n_samples defined twice?
  strucvec <- rep(c(0,1,0, 0,1,1, 1,0,0, 1,0,1, 1,1,1), times = n_samples/5)
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
  
  
  Datafit = coefficients %*% t(cbind(rep(1, 5*n_rep), X)) #matrix multiplication to fit model
  Residual = ratios - Datafit
  Expression_variation = rowMeans(ratios^2)
  Expression_residual = rowMeans(Residual^2)
  Varexplain = 100*(Expression_variation - Expression_residual) / Expression_variation
  
  Fit <- data.frame(Expression_variation, Varexplain)

  
  outputdf <- data.frame(cbind(coefficients,Pvalue,rsquared))
  colnames(outputdf)[9] <- "R2"
  return(outputdf)
}
