stopifnot("binom.kind" %in% names(formals(RNGkind)))

binom_kind <- "Buggy BTPE"

draw_binom <- function(N, size, prob) {
  suppressWarnings(RNGkind(binom.kind = binom_kind))
  rbinom(N, size, prob)
}

gof_pvalue <- function(y, n, p, min_exp = 5) {
  N <- length(y)
  ok <- which(N * dbinom(0:n, n, p) >= min_exp) - 1
  lo <- min(ok)
  hi <- max(ok)
  probs <- c(pbinom(lo, n, p),
             dbinom((lo + 1):(hi - 1), n, p),
             pbinom(hi - 1, n, p, lower.tail = FALSE))
  y <- pmin(pmax(y, lo), hi)
  chisq.test(tabulate(y - lo + 1, hi - lo + 1), p = probs)$p.value
}

size_binom <- 400
prob_binom <- 0.5
n_draws <- 5e7

test_that("rbinom matches the binomial distribution", {
  expect_pval(function() gof_pvalue(draw_binom(n_draws, size_binom, prob_binom), size_binom, prob_binom),
              name = "rbinom_gof", power = 0.1, U0 = 1e6, amax = 10, cbar0max = 1.1, p_leave = 0.05)
})
