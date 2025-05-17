# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#Do we need 5 fold changes for the model?
# Is it still relevant if there are many low P-values

# Data should contain 4 conditions for 2 signals
# Ctrl, LPS, pH, LPS_pH

#format: Ctrl_A, LPS_A, pH_A, LPS_pH_A, Ctrl_B etc.

#' Geometric mean of a vector
#'
#' `geomean` returns the geometric mean of a vector x
#' @param x should be a numeric vector

geomean <- function(x){
  exp(mean(log(x)))
}

#' Median Normalization of the data
#'
#'  `MedianNorm` normalizes data with the median normalization method
#'
#' @param data Use a count matrix K_ij, with one row for each gene i and one column for each sample j. Matrix entries indicate the number of sequencing reads mapped to a gene in a sample.
#' @param countthres Count threshold for the count data. Default is 10
#' @param pseudo Pseudo count added to deal with zero count problem. Default is 1

MedianNorm <- function(data, countthres=10, pseudo=1){
  Subdata = data[rowMeans(data)>=countthres,]     # Subset rows with average count higher than 10
  Subdata = Subdata + pseudo                       # Add Pseudo count of 1

  t = apply(Subdata,1, geomean)         # Calculate Geometric mean
  # Divide expression of every element per row to the geometric mean of that row

  DataRatio = log2(Subdata / t)
  T_med = apply(DataRatio,2, median) # median per column of array
  T_med = 2^T_med

  C <- t(t(Subdata) / T_med)
  return(C)
}

