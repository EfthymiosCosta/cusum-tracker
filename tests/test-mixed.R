mu_0 <- 0
sigma <- 1
n <- 100

pval_gen_mean <- function() {
  x <- rnorm(n = n, mean = mu_0, sd = sigma)
  z <- sqrt(n) * (mean(x) - mu_0) / sigma
  pnorm(z, lower.tail = FALSE)
}

pval_gen_var <- function() {
  x <- rnorm(n = n, mean = mu_0, sd = sigma)
  test_stat <- (n - 1) * var(x) / sigma^2
  pchisq(test_stat, df = n - 1, lower.tail = FALSE)
}

pval_gen_var_wrong <- function() {
  x <- rnorm(n = n, mean = mu_0, sd = sigma * 2)
  test_stat <- (n - 1) * var(x) / sigma^2
  pchisq(test_stat, df = n - 1, lower.tail = FALSE)
}

# Deterministic test
test_that("arithmetic works", {
  expect_equal(2 + 2, 4)
  expect_equal(sqrt(16), 4)
})

# Non-deterministic test
test_that("mean stays in control", {
  expect_pval(pval_gen_mean, name = "mixed_mean_check", num_resims = 10)
})

# Combination
test_that("combined deterministic and non-deterministic checks", {
  expect_equal(length(rnorm(n)), n)
  expect_pval(pval_gen_var_wrong, name = "mixed_var_wrong", num_resims = 10)
  expect_pval(pval_gen_var, name = "mixed_var_check", num_resims = 10)
})