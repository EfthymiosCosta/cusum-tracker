# run_tests.R
source("R/cusum_funs.R")

mu_0 <- 0
sigma <- 1

pval_gen <- function(){
  x_sample <- rnorm(n = 100, mean = mu_0, sd = sigma)
  z <- sqrt(length(x_sample)) * (mean(x_sample) - mu_0) / sigma
  pnorm(z, lower.tail = FALSE)
}

result <- expect_pval(
  pval_gen, 
  name = "drift_check", 
  dir = "cusum_logs", 
  num_resims = 10,
  on_signal = "message"
)
