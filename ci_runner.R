# ci_runner.R
source("R/cusum_funs.R")
library(testthat)

log_dir <- "cusum_logs"
test_dir <- "tests"

# Snapshot the logs dir
files_before <- list.files(log_dir, pattern = "\\.rds$",
                           full.names = TRUE, recursive = TRUE)
mtimes_before <- setNames(file.mtime(files_before), files_before)

# Run tests without crashing
test_scripts <- list.files(test_dir, pattern = "\\.[Rr]$", full.names = TRUE)
if (length(test_scripts) == 0L)
  stop("No test scripts found in '", test_dir, "/'.", call. = FALSE)

source_errors <- list()
for (script in test_scripts) {
  tryCatch(
    source(script),
    error = function(e) {
      source_errors[[script]] <<- conditionMessage(e)
    }
  )
}

# Snapshot after
files_after <- list.files(log_dir, pattern = "\\.rds$",
                           full.names = TRUE, recursive = TRUE)
mtimes_after <- setNames(file.mtime(files_after), files_after)

is_new <- !(files_after %in% files_before)
is_updated <- !is_new & (mtimes_after > mtimes_before[files_after])
touched <- files_after[is_new | is_updated]

# Report charts that signal
all_passed <- TRUE

for (script in names(source_errors)) {
  passed <- test_that(sprintf("script '%s' ran without error", basename(script)), {
    fail(source_errors[[script]])
  })
  all_passed <- all_passed && passed
}

for (lf in touched) {
  log <- tryCatch(read_pval_log(lf), error = function(e) NULL)
  if (is.null(log) || nrow(log$data) == 0) next
  
  last <- log$data[nrow(log$data), ]
  test_name <- if (!is.null(log$config$test_name) && nzchar(log$config$test_name)) {
    log$config$test_name
  } else {
    tools::file_path_sans_ext(basename(lf))
  }
  
  desc <- sprintf("%s  (S_t = %.4f, h = %.4f, last p = %.3g)",
                  test_name, last$S_t, log$config$h, last$p_value)
  
  passed <- test_that(desc, {
    expect_lt(last$S_t, log$config$h)
  })
  all_passed <- all_passed && passed
}

if (!all_passed) {
  stop("CUSUM run had failures, see test output above.", call. = FALSE)
}