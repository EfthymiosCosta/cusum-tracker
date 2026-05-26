# ci_runner.R
source("R/cusum_funs.R")
source("run_tests.R")

log_files <- list.files("cusum_logs", pattern = "\\.rds$", full.names = TRUE)

signals <- list()
for (lf in log_files) {
  log <- tryCatch(read_pval_log(lf), error = function(e) NULL)
  if (is.null(log) || nrow(log$data) == 0L) next
  
  last <- log$data[nrow(log$data), ]
  if (isTRUE(last$signal)) {
    nm <- if (!is.null(log$config$test_name) && nzchar(log$config$test_name))
      log$config$test_name
    else tools::file_path_sans_ext(basename(lf))
    signals[[nm]] <- list(S_t = last$S_t,
                          h   = log$config$h,
                          p   = last$p_value)
  }
}

if (length(signals) > 0L) {
  parts <- vapply(names(signals), function(nm) {
    s <- signals[[nm]]
    sprintf("`%s` (S_t = %.4f, h = %.4f, last p = %.3g)",
            nm, s$S_t, s$h, s$p)
  }, character(1))
  
  env_file <- Sys.getenv("GITHUB_ENV")
  if (nzchar(env_file)) {
    cat("CUSUM_SIGNAL=true\n", file = env_file, append = TRUE)
    cat(sprintf("CUSUM_MSG=Drift detected on: %s\n",
                paste(parts, collapse = "; ")),
        file = env_file, append = TRUE)
  }
}