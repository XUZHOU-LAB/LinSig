library(dplyr)
library(MASS)
source("~/Boston Internship/Github/Rsyn/paper_figures/data_standardization.R")


generate_synthetic_data <- function(sourcedf, size=20000, nrep=2){
  normed <- MedianNorm(sourcedf[rowMeans(sourcedf)>10,]) # filter out low counts
  
  col_ctrl <- 1:nrep + 0*nrep # col 1,2 if 2 replicates
  col_condA <- 1:nrep + 1*nrep # col 3,4 if 2 replicates
  col_condB <- 1:nrep + 2*nrep  # col 5,6 if 2 replicates
  col_condAB <- 1:nrep + 3*nrep  # col 7,8 if 2 replicates
  
  lmucov <- compute_log_mu_cov(normed, nrep=nrep)
  
  
  # randomly draw a coefficient of variation within a bin
  coefs_of_variation <- sample_cov_by_binned_mu(lmucov, size=20000, n_bins=20)
  
  
  # Generate new row means
  mustat          <- MASS::fitdistr(rowMeans(normed), "lognormal")
  new_row_means   <- rlnorm(size, mustat$estimate[1], mustat$estimate[2])
  onrm            <- sort(new_row_means)
  
  
  sampledMuCoV <- data.frame(
    mu  = log(onrm),
    cov = coefs_of_variation
  )
  
  # draw new counts from normal distribution using mu (row mean) and CoV
  synthetic_dataset <- t(
      apply(sampledMuCoV, 1, FUN=drawreps, nrep=nrep))
  
  
  return(round(synthetic_dataset, 1))
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
  log_mucov <- log(mucov)
  
  return(log_mucov)
}


sample_cov_by_binned_mu <- function(log_mucov_df, size = 20000, n_bins = 20) {
  # Ensure mu column is sorted
  
  sorted_df <- log_mucov_df[order(log_mucov_df$mu), ]
  
  
  # Bin by mu into equal-width intervals
  binned_df <- dplyr::mutate(
    sorted_df,
    mu_bin = cut(mu, breaks = n_bins)
  )
  

  # Sample cov values proportionally within each bin
  ncovs <- unlist(lapply(split(binned_df, binned_df$mu_bin), function(bin_df) {
    bin_size <- nrow(bin_df)
    n_samples <- round((bin_size / nrow(binned_df)) * size)
    sample(bin_df$cov, n_samples, replace = TRUE)
  }))
  

  # Ensure exactly `size` values
  if (length(ncovs) > 20000) {
    ncovs <- ncovs[seq_len(20000)]
  } else if (length(ncovs) < 20000) {
    # if too few, sample additional with replacement from itself
    ncovs <- c(ncovs, sample(ncovs, 20000 - length(ncovs), replace = TRUE))
  }
  
  return(unname(ncovs))
}


drawreps <- function(mu_cov_vec, nrep = 2) {
  mu <- mu_cov_vec["mu"]
  cov <- mu_cov_vec["cov"]
  
  random_counts <- rnorm(
    n = 4 * nrep,
    mean = exp(mu),
    sd = exp(mu + cov)
  )
  random_counts[random_counts < 0] <- 1
  return(random_counts)
}


df <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/IL6IL10combDF.csv", row.names=1)[,1:8]
dfs <- df[c(51,45,530,432,825,1034,5440),]
gen_randomstats(df, nrep = 2, size = 20000)
