.append_pval_log <- function(file_path, objs) {
  lines <- vapply(
    objs,
    function(o) as.character(jsonlite::toJSON(o, auto_unbox = TRUE, digits = NA, na = "null")),
    character(1)
  )
  con <- file(file_path, open = "a", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con)
  invisible(file_path)
}

read_pval_header <- function(file_path) {
  if (!file.exists(file_path))
    stop("Log file not found: ", file_path, call. = FALSE)
  con <- file(file_path, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  first <- readLines(con, n = 1, warn = FALSE)
  if (length(first) == 0 || !nzchar(first))
    stop("Log '", file_path, "' is empty or has no header.", call. = FALSE)
  cfg <- jsonlite::fromJSON(first)
  if (!all(c("delta", "h", "H", "v", "amax") %in% names(cfg)))
    stop("First line of '", file_path, "' is not a CUSUM design header.", call. = FALSE)
  cfg
}

read_pval_last <- function(file_path) {
  if (!file.exists(file_path))
    stop("Log file not found: ", file_path, call. = FALSE)
  lines <- readLines(file_path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  if (length(lines) <= 1)
    return(NULL)
  jsonlite::fromJSON(lines[length(lines)])
}

read_pval_log <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  header <- read_pval_header(file_path)
  entries <- lapply(lines[-1], jsonlite::fromJSON)
  data <- if (length(entries)) do.call(rbind, lapply(entries, as.data.frame)) else
    data.frame(t = character(0), S = numeric(0), steps = integer(0), signal = integer(0))
  list(header = header, data = data)
}

.draw_pvalue <- function(fn) {
  p <- fn()
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p < 0 || p > 1)
    stop("`p_value()` must return a single number in [0, 1].", call. = FALSE)
  p
}

expect_pval <- function(p_value,
                        name,
                        dir = getOption("cusum.dir", "."),
                        alpha = 0.05,
                        power = 0.2,
                        U0 = 1e5,
                        cbar0max = 1.1,
                        amax = 10,
                        p_leave = 0.01,
                        n_bins = 151,
                        n_scan = 30,
                        design_dir = file.path(getOption("cusum.dir", "."), "designs"),
                        quiet = FALSE) {
  if (!is.function(p_value))
    stop("`p_value` must be a function that returns a single p-value.", call. = FALSE)
  if (!is.character(name) || length(name) != 1L || !nzchar(name))
    stop("`name` must be a single non-empty string.", call. = FALSE)
  file_name <- if (grepl("\\.jsonl$", name, ignore.case = TRUE)) name else paste0(name, ".jsonl")
  file_path <- file.path(dir, file_name)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    if (!quiet) message("expect_pval(): created directory '", dir, "'.")
  }
  new_log <- !file.exists(file_path)
  if (new_log) {
    design <- cusum_design_cached(alpha, power, U0, cbar0max, amax, p_leave,
                                  n_bins, n_scan, design_dir, quiet)
    header <- c(list(name = name, created = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), design)
    S_prev <- 0
  } else {
    header <- read_pval_header(file_path)
    last <- read_pval_last(file_path)
    S_prev <- if (is.null(last)) 0 else last$S
  }
  S <- S_prev
  steps <- 0L
  repeat {
    steps <- steps + 1L
    p <- .draw_pvalue(p_value)
    S <- min(header$H, max(0, S + qnorm(p, lower.tail = FALSE) - header$delta / 2))
    if (steps >= header$amax || S < header$v) break
  }
  signal <- S >= header$h
  crossed <- signal && S_prev < header$h
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- list(t = timestamp, S = S, steps = steps, signal = as.integer(signal))
  .append_pval_log(file_path, if (new_log) list(header, entry) else list(entry))
  if (new_log && !quiet) message("expect_pval(): created log '", file_path, "'.")
  failure_message <- sprintf(
    "CUSUM chart '%s' %s the signalling threshold: S = %.4f >= h = %.4f after %d computational step%s.",
    name, if (crossed) "has crossed" else "remains above", S, header$h, steps,
    if (steps == 1L) "" else "s")
  testthat::expect(ok = !signal, failure_message = failure_message)
  invisible(list(name = name, file = file_path, timestamp = timestamp,
                 S_prev = S_prev, S = S, steps = steps, signal = signal,
                 crossed = crossed, h = header$h, H = header$H, v = header$v))
}

cusum_group <- function(description, code) {
  group_name <- gsub("[^A-Za-z0-9_-]+", "_", description)
  group_name <- gsub("^_+|_+$", "", group_name)
  if (!nzchar(group_name))
    stop("`description` produced an empty group name.", call. = FALSE)
  group_dir <- file.path(getOption("cusum.dir", "."), group_name)
  dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
  local_env <- new.env(parent = parent.frame())
  local_env$expect_pval <- function(p_value, name, dir = group_dir, ...) {
    expect_pval(p_value = p_value, name = name, dir = dir, ...)
  }
  eval(substitute(code), envir = local_env)
}
