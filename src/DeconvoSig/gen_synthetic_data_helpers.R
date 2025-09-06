library(dplyr)
library(MASS)
source("~/bloop/LinSig/src/DeconvoSig/basic_model_functions/data_standardization.R")

# Generate multiple synthetic datasets
generate_multiple_datasets <- function(source_df, nrep = 2, size, n_datasets = 10) {
  replicate(n_datasets, generate_synthetic_data(source_df = source_df, nrep = nrep, size = size), simplify = FALSE)
}

generate_synthetic_data <- function(source_df, size, nrep=2,
                                    count_threshold=10){
  
  sampledMuCoV <- generate_mucov_df(source_df=source_df,size=size,nrep=nrep,
                                    count_threshold=count_threshold)
  
  # draw new counts from normal distribution using mu (row mean) and CoV
  synthetic_dataset <- t(
    apply(sampledMuCoV, 1, FUN=drawreps, nrep=nrep))
  
  
  return(round(synthetic_dataset, 1))
}



generate_mucov_df <- function(source_df, size, nrep=2,
                              count_threshold=10){
  
  normed <- MedianNorm(source_df, count_threshold = count_threshold) # filter out low counts
  
  col_ctrl <- 1:nrep + 0*nrep # col 1,2 if 2 replicates
  col_condA <- 1:nrep + 1*nrep # col 3,4 if 2 replicates
  col_condB <- 1:nrep + 2*nrep  # col 5,6 if 2 replicates
  col_condAB <- 1:nrep + 3*nrep  # col 7,8 if 2 replicates
  
  lmucov <- compute_log_mu_cov(normed, nrep=nrep)
  

  # Estimate 2D kernel density
  kde <- kde2d(lmucov$mu, lmucov$cov, n = 100)
  
  # Sample from the estimated 2D density
  dens_vals <- as.vector(kde$z)
  dens_vals <- dens_vals / sum(dens_vals)  # normalize to make it a probability
  
  # Sample indices according to the density
  sample_indices <- sample(length(dens_vals), size = size, replace = TRUE, prob = dens_vals)
  
  # Convert indices back to x, y values
  grid_x <- rep(kde$x, times = length(kde$y))     # x varies fastest
  grid_y <- rep(kde$y, each = length(kde$x))      # y varies slowest
  
  sampled_x <- grid_x[sample_indices] + runif(size, -diff(kde$x)[1]/2, diff(kde$x)[1]/2)
  sampled_y <- grid_y[sample_indices] + runif(size, -diff(kde$y)[1]/2, diff(kde$y)[1]/2)
  
  sampledMuCoV <- data.frame(mu = sampled_x, cov = sampled_y)
  
  #plot(sampledMuCoV$mu, sampledMuCoV$cov, pch='.', xlim=c(0,12),ylim=c(-14,0))
  return(sampledMuCoV)
}



# compute per‐condition means (between replicates)
compute_condition_means <- function(df, nrep = 2) {
  # Build the column indices for each condition
  col_groups <- list(
    col_ctrl   <- seq_len(nrep),                    # 1:nrep
    col_condA  <- seq_len(nrep) +     nrep,         # (nrep+1):(2*nrep)
    col_condB  <- seq_len(nrep) + 2 * nrep,         # (2*nrep+1):(3*nrep)
    col_condAB <- seq_len(nrep) + 3 * nrep          # (3*nrep+1):(4*nrep)
  )
  
  
  # Compute rowMeans for each group
  mean_list <- lapply(col_groups, function(cols) {
    rowMeans(df[, cols, drop = FALSE])
  })
  
  # Combine into a matrix (rows = features, cols = conditions)
  mean_mat <- do.call(cbind, mean_list)
  
  return(mean_mat)
}


compute_condition_sd <- function(df, nrep = 2) {
  # Define column groups for each condition
  col_groups <- list(
    control      = seq_len(nrep),
    condition_A  = seq_len(nrep) +     nrep,
    condition_B  = seq_len(nrep) + 2 * nrep,
    condition_AB = seq_len(nrep) + 3 * nrep
  )
  
  # Compute row-wise standard deviations for each group
  sd_list <- lapply(col_groups, function(cols) {
    row_means <- rowMeans(df[, cols, drop = FALSE])
    row_sumsq <- rowSums((df[, cols, drop = FALSE] - row_means)^2)
    sqrt(row_sumsq)
  })
  
  # Combine into a matrix
  sd_mat <- do.call(cbind, sd_list)
  colnames(sd_mat) <- names(col_groups)
  return(sd_mat)
}


compute_log_mu_cov <- function(df, nrep = 2) {
  # Compute per-condition means and standard deviations
  mean_per_condition <- compute_condition_means(df, nrep = nrep)
  sd_per_condition   <- compute_condition_sd(df, nrep = nrep)
  
  # Compute coefficient of variation (sd / mean)
  mucov <- data.frame(
    mu  = as.vector(mean_per_condition),
    cov = as.vector(sd_per_condition / mean_per_condition)
  )
  

  # Log-transform mu and cov
  log_mucov <- log2(mucov)
  
  return(log_mucov)
}


drawreps <- function(mu_cov_vec, nrep = 2) {
  mu <- 2^mu_cov_vec["mu"]
  cov <- 2^mu_cov_vec["cov"]#*0.5

  random_counts <- rnorm(
    n = 4 * nrep,
    mean = mu,#exp(mu)
    sd = mu * cov #exp(mu + cov)
  )
  random_counts <- random_counts * rlnorm(1, meanlog = 0, sdlog = 0)#sdlog = 0.10
  
  random_counts[random_counts < 0] <- 1
  return(random_counts)
}
