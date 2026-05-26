# ci_runner.R
source("run_tests.R")

# GitHub Actions communication
if (test_result$signal) {
  env_file <- Sys.getenv("GITHUB_ENV")
  
  if (nzchar(env_file)) {
    cat("CUSUM_SIGNAL=true\n", file = env_file, append = TRUE)
    
    msg <- sprintf("Drift detected on test `%s`. Final S_t = %.4f exceeds threshold h = %.4f.",
                   test_result$name, test_result$S_t, test_result$h)
    cat(sprintf("CUSUM_MSG=%s\n", msg), file = env_file, append = TRUE)
  }
}