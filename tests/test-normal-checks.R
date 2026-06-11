mu_0 <- 0
sigma <- 1
n <- 100

pval_gen_mean <- function() {
  x_sample <- rnorm(n = n, mean = mu_0, sd = sigma)
  z <- sqrt(length(x_sample)) * (mean(x_sample) - mu_0) / sigma
  pnorm(z, lower.tail = FALSE)
}

pval_gen_var <- function() {
  x <- rnorm(n = n, mean = mu_0, sd = sigma)
  test_stat <- (n - 1) * var(x) / sigma^2
  pchisq(test_stat, df = n - 1, lower.tail = FALSE)
}

test_that("normal distribution checks", {
  expect_pval_jsonl(pval_gen_mean, name = "mean_check", num_resims = 10)
  expect_pval_jsonl(pval_gen_var, name = "variance_check", num_resims = 10)
})