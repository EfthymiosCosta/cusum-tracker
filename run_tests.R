# run_tests.R
source("R/cusum_funs.R")

mu_0 <- 0
sigma <- 1
n <- 100

pval_gen <- function(){
  x_sample <- rnorm(n = n, mean = mu_0, sd = sigma)
  z <- sqrt(length(x_sample)) * (mean(x_sample) - mu_0) / sigma
  pnorm(z, lower.tail = FALSE)
}

pval_gen_var <- function() {
  x <- rnorm(n = n, mean = mu_0, sd = sigma)
  test_stat <- (n - 1) * var(x) / sigma^2
  pchisq(test_stat, df = n - 1, lower.tail = FALSE)
}

expect_pval(
  pval_gen, 
  name = "drift_check", 
  dir = "cusum_logs", 
  num_resims = 10,
  on_signal = "message"
)

expect_pval(pval_gen_var, 
            name = "variance_check",
            dir = "cusum_logs",
            num_resims = 10,
            on_signal = "message")
