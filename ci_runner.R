source("R/cusum_funs.R")
library(testthat)
options(cusum.dir = normalizePath("cusum_logs", mustWork = FALSE))
testthat::test_dir("tests", stop_on_failure = TRUE)