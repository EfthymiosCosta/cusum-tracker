source("R/cusum_funs.R")
library(testthat)
options(cusum.dir = "cusum_logs")
testthat::test_dir("tests", stop_on_failure = TRUE)