library(dplyr)
library(MASS)
library(cellsigsyn)
source("~/Boston Internship/Github/Rsyn/paper_figures/data_standardization.R") # for MedianNorm function
source("~/Boston Internship/Github/Rsyn/paper_figures/gen_synthetic_data_helpers.R")
source("~/Boston Internship/cellsigsyn/R/compute_ratios.R") # for RNAseqLowess function
library(parallel)


### PARAMETERS ###
params <- new.env(parent = emptyenv())
params$struc_vec <- c(0,1,0, # 
                      0,1,1, #
                      1,0,0, #
                      1,0,1, #
                      1,1,1) #

params$x_axis_thresholds <- seq(0.001, 4.001, 0.001) # initialize an X-axis for the LFC threshold figure (LFC 0.001-4.001, with a small step size for high resolution curve)


strucDF <- list(
  col_ctrl = c(1, 2),      # Columns for CTRL (e.g., replicates 1 and 2)
  col_condA = c(3, 4),     # Columns for Condition A (e.g., replicates 3 and 4)
  col_condB = c(5, 6),     # Columns for Condition B (e.g., replicates 5 and 6)
  col_condAB = c(7, 8)     # Columns for Condition AB (e.g., replicates 7 and 8)
)

# Helper function to compute FDR for a given column
compute_fdr <- function(df, column, thresholds) {
  sapply(thresholds, function(t) mean(abs(df[[column]]) > t))
}


# Helper function to extract threshold for closest FDR to target
get_threshold <- function(FDR_values, target) {
  params$x_axis_thresholds[which.min(abs(FDR_values - target))]
}


compute_lfc_thresholds <- function(source_df, size=20000, nrep=2, lowessn=0, onlyDF=0){
  
  gtt <- generate_synthetic_data(source_df=source_df, size=size, nrep=2)
  
  #median normalization
  normedR <- MedianNorm(gtt)
  
  #compute ratios
  #rats <- compute_ratiosR(normedR, nreps=nrep, lowess_norm=lowessn)
  strucDF <- list(
    col_ctrl = c(1, 2),      # Columns for CTRL (e.g., replicates 1 and 2)
    col_condA = c(3, 4),     # Columns for Condition A (e.g., replicates 3 and 4)
    col_condB = c(5, 6),     # Columns for Condition B (e.g., replicates 5 and 6)
    col_condAB = c(7, 8)     # Columns for Condition AB (e.g., replicates 7 and 8)
  )
  
  rats<- compute_ratios(normedR, pseudo=1, structureDataFrame=strucDF)
  
  
  #fit model
  stats <- ratios.fit(rats, CompThreshold=1,n_rep=nrep)
  
  # Flag significant hits for A, B, and AB based on p-value and R²
  stats$SIGA  <- stats[, 6] < 0.05 & stats$R2 > 0.8
  stats$SIGB  <- stats[, 7] < 0.05 & stats$R2 > 0.8
  stats$SIGAB <- stats[, 8] < 0.05 & stats$R2 > 0.8
  
  # Subset significant rows
  FDRA  <- stats[stats$SIGA, ]
  FDRB  <- stats[stats$SIGB, ]
  FDRAB <- stats[stats$SIGAB, ]
  
  # Define thresholds

  
  # Compute FDR curves
  FDRa  <- compute_fdr(FDRA,  "X1", params$x_axis_thresholds)
  FDRb  <- compute_fdr(FDRB,  "X2", params$x_axis_thresholds)
  FDRab <- compute_fdr(FDRAB, "X3", params$x_axis_thresholds)
  
  # Plot all
  plot(params$x_axis_thresholds, FDRab, log = "x", type = "l", col = "green", lwd = 1.5,
       xlim = c(0.001, 2), ylim = c(0, 0.99), ylab = "FDR", xlab = "LFC Threshold")
  lines(params$x_axis_thresholds, FDRb,  col = "orange", lwd = 1.5)
  lines(params$x_axis_thresholds, FDRa,  col = "blue",   lwd = 1.5)
  abline(h = 0.05, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = c("A", "B", "AB"), col = c("blue", "orange", "green"), lwd = 2)
  
  
  FDRDF <- data.frame(Xs = params$x_axis_thresholds, 
                      FDRa = FDRa, 
                      FDRb = FDRb, 
                      FDRab = FDRab)
  
  
  # Create summary data
  summary_table <- data.frame(
    Condition          = c("A", "B", "AB"),
    Significant_genes   = c(nrow(FDRA), nrow(FDRB), nrow(FDRAB)),
    Threshold_at_5pct  = c(get_threshold(FDRa, 0.05),
                           get_threshold(FDRb, 0.05),
                           get_threshold(FDRab, 0.05)),
    Threshold_at_10pct = c(get_threshold(FDRa, 0.10),
                           get_threshold(FDRb, 0.10),
                           get_threshold(FDRab, 0.10)),
    row.names = c("A", "B", "AB")
  )
  
  
  # Print summary table
  return(summary_table)
}


#test function
#cts <- read.csv("C:/Users/HB/OneDrive/Documents/Boston Internship/IL6IL10combDF.csv", row.names=1)[,1:8]
#compute_lfc_thresholds(source_df = cts, nrep=2)
