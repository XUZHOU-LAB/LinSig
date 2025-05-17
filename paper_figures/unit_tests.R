library(testthat)

source("~/Boston Internship/Github/Rsyn/paper_figures/gen_synthetic_data_helpers.R")

test_that("compute_condition_means returns correct row-means", {
  
  test_df <- matrix(data=1:8, nrow = 1, ncol=8)
  
  out <- compute_condition_means(test_df, nrep = 2)
  
  # Expect a 2×4 matrix:
  expect_equal(dim(out), c(1, 4))

  # Since we filled blocks as constants, each mean should equal that block value
  # For both rows, control = 1, A = 2, B = 3, AB = 4:
  expected <- matrix(ncol=4,
    data=c(mean(c(1,2)), mean(c(3,4)), mean(c(5,6)), mean(c(7,8)))
  )

  expect_equal(out, expected)
})


test_that("compute_condition_sd returns correct row-sds", {
  # 1 row, 8 columns: values 1..8
  test_df <- matrix(data = 1:8, nrow = 1, ncol = 8)
  
  out <- compute_condition_sd(test_df, nrep = 2)
  
  # Expect 1×4 matrix
  expect_equal(dim(out), c(1, 4))
  expect_equal(colnames(out),
               c("control", "condition_A", "condition_B", "condition_AB"))
  
  # Manually compute each block’s SD:
  # control = sd(c(1,2)) = sqrt(( (1-1.5)^2 + (2-1.5)^2 ))
  sd_ctrl  <- sqrt((1 - 1.5)^2 + (2 - 1.5)^2)
  sd_A     <- sqrt((3 - 3.5)^2 + (4 - 3.5)^2)
  sd_B     <- sqrt((5 - 5.5)^2 + (6 - 5.5)^2)
  sd_AB    <- sqrt((7 - 7.5)^2 + (8 - 7.5)^2)
  
  expected <- matrix(
    data = c(sd_ctrl, sd_A, sd_B, sd_AB),
    nrow = 1,
    byrow = TRUE
  )
  colnames(expected) <- colnames(out)
  
  expect_equal(out, expected)
})




test_that("compute_log_mu_cov returns correct log(mu) and log(cov)", {
  
  # Simulated input matrix with 2 replicates per condition, 1 row
  # Columns: ctrl(1,2), condA(3,4), condB(5,6), condAB(7,8)
  test_df <- matrix(c(
    10, 20,   # ctrl
    30, 40,   # A
    50, 60,   # B
    70, 80    # AB
  ), nrow = 1)
  
  # Expected means for each condition
  expected_means <- c(15, 35, 55, 75)
  
  # Since all values in each condition are equal, SD should be 0
  expected_sds <- c(sd(c(10,20)), sd(c(30,40)), sd(c(50,60)), sd(c(70,80)))
  
  # cov = sd / mean → 0 / x = 0
  expected_cov <- expected_sds / expected_means
  
  # log-transformed expected results
  expected_log_mu <- log(expected_means)
  expected_log_cov <- log(expected_cov)
  
  # log(0) is -Inf, so handle that explicitly
  expected <- data.frame(
    mu  = expected_log_mu,
    cov = expected_log_cov
  )
  
  # Run function
  result <- compute_log_mu_cov(test_df, nrep = 2)
  
  # Test structure
  expect_equal(dim(result), c(4, 2))
  expect_equal(colnames(result), c("mu", "cov"))
  
  # Test content
  expect_equal(result$mu, expected$mu)
  expect_equal(result$cov, expected$cov)
})


