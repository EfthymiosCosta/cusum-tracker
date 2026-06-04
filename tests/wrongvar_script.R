# run_tests.R
source("R/cusum_funs.R")

mu_0 <- 0
sigma <- 1
n <- 100


pval_gen_var_wrong <- function() {
  x <- rnorm(n = n, mean = mu_0, sd = sigma*1.5)
  test_stat <- (n - 1) * var(x) / sigma^2
  pchisq(test_stat, df = n - 1, lower.tail = FALSE)
}


expect_pval(pval_gen_var_wrong, name = "wrong_var",
            on_signal = "message")