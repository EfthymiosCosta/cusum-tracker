.append_pval_log <- function(file_path, objs) {
  lines <- vapply(
    objs,
    function(o) as.character(
      jsonlite::toJSON(o, auto_unbox = TRUE, digits = NA, na = "null")
    ),
    character(1)
  )
  con <- file(file_path, open = "a", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con)
  invisible(file_path)
}

#' Read the first line of a JSON Lines CUSUM log
#' @export
read_pval_header <- function(file_path) {
  if (!file.exists(file_path))
    stop("Log file not found: ", file_path, call. = FALSE)
  con <- file(file_path, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  first <- readLines(con, n = 1, warn = FALSE)
  if (length(first) == 0 || !nzchar(first))
    stop("Log '", file_path, "' is empty or has no header.", call. = FALSE)
  cfg <- jsonlite::fromJSON(first)
  if (is.null(cfg$alpha) || is.null(cfg$h))
    stop("First line of '", file_path, "' is not a CUSUM config header.",
         call. = FALSE)
  cfg
}

#' Read the latest observation (last line) of a JSON Lines CUSUM log
#'
#' Uses `tail` so the cost does not grow with the number of records. Returns
#' `NULL` when the log holds only its configuration header.
#' @export
read_pval_last <- function(file_path) {
  if (!file.exists(file_path))
    stop("Log file not found: ", file_path, call. = FALSE)
  last2 <- system2("tail", c("-n", "2", shQuote(file_path)),
                   stdout = TRUE, stderr = FALSE)
  if (length(last2) <= 1)
    return(NULL)
  jsonlite::fromJSON(last2[length(last2)])
}

expect_pval_jsonl <- function(p_value,
                              name,
                              dir = getOption("cusum.dir", "."),
                              alpha = 5e-2,
                              power = 2e-1,
                              target_ARL0 = 1000,
                              M = 50,
                              h = NULL,
                              num_resims = 10,
                              quiet = FALSE) {
  
  p_value_fn <- NULL
  if (is.function(p_value)) {
    p_value_fn <- p_value
    p_value <- p_value()
  }
  if (!is.numeric(p_value) || length(p_value) != 1 || !is.finite(p_value))
    stop("The provided or generated p-value must be a single finite numeric value.", call. = FALSE)
  if (p_value < 0 || p_value > 1)
    stop("The provided or generated p-value must lie in [0, 1]; got ", p_value, ".", call. = FALSE)
  if (!is.character(name) || length(name) != 1 || !nzchar(name))
    stop("`name` must be a single non-empty string.", call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be a single numeric value in (0, 1).", call. = FALSE)
  if (!is.numeric(power) || length(power) != 1 || power <= 0 || power >= 1)
    stop("`power` must be a single numeric value in (0, 1).", call. = FALSE)
  if (!is.numeric(num_resims) || length(num_resims) != 1 || num_resims < 0 || num_resims %% 1 != 0)
    stop("`num_resims` must be a single non-negative integer.", call. = FALSE)
  
  file_name <- if (grepl("\\.jsonl$", name, ignore.case = TRUE))
    name else paste0(name, ".jsonl")
  file_path <- file.path(dir, file_name)
  
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    if (!quiet) message("expect_pval(): created directory '", dir, "'.")
  }
  
  append_objs <- list()
  wrote_config <- FALSE
  
  if (!file.exists(file_path)) {
    delta <- qnorm(1 - alpha) - qnorm(1 - power)
    if (is.null(h))
      h <- estimate_h(target_ARL0 = target_ARL0, delta = delta, M = M)
    H <- 2 * h
    ARL0 <- compute_ARL(h = h, delta = delta, M = M, is_H0 = TRUE)
    ARL1 <- compute_ARL(h = h, delta = delta, M = M, is_H0 = FALSE)
    config <- list(
      test_name = name,
      created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      alpha = alpha, power = power, M = M, target_ARL0 = target_ARL0,
      h = h, H = H, ARL0 = ARL0, ARL1 = ARL1
    )
    append_objs <- c(append_objs, list(config))
    if (!quiet) message("expect_pval(): created log '", file_path, "'.")
    wrote_config <- TRUE
    S_prev <- 0
  } else {
    config <- read_pval_header(file_path)
    if (is.null(config$alpha) || is.null(config$power) ||
        is.null(config$h) || is.null(config$H))
      stop("'", file_path, "' is missing CUSUM parameters in its header; ",
           "it does not look like an expect_pval() log.", call. = FALSE)
    last <- read_pval_last(file_path)
    S_prev <- if (is.null(last)) 0 else last$S_t
    alpha <- config$alpha
    power <- config$power
    delta <- qnorm(1 - alpha) - qnorm(1 - power)
    h <- config$h
    H <- config$H
  }
  
  increment <- compute_LAN_increment(p_value, delta)
  S_t <- min(H, max(0, S_prev + increment))
  signal <- S_t >= h
  crossed <- signal && (S_prev < h)
  saturated <- S_t >= H
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  append_objs <- c(append_objs, list(list(
    timestamp = timestamp, p_value = p_value,
    increment = increment, S_t = S_t, signal = signal
  )))
  
  if (signal && !is.null(p_value_fn) && num_resims > 0) {
    for (i in seq_len(num_resims)) {
      S_prev <- S_t
      p_new <- p_value_fn()
      if (!is.numeric(p_new) || length(p_new) != 1L || !is.finite(p_new) || p_new < 0 || p_new > 1)
        stop("Generated p-value during continuation steps is invalid.", call. = FALSE)
      increment <- compute_LAN_increment(p_new, delta)
      S_t <- min(H, max(0, S_prev + increment))
      signal <- S_t >= h
      crossed <- signal && (S_prev < h)
      saturated <- S_t >= H
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      append_objs <- c(append_objs, list(list(
        timestamp = timestamp, p_value = p_new,
        increment = increment, S_t = S_t, signal = signal
      )))
    }
  }
  
  .append_pval_log(file_path, append_objs)
  
  result <- list(
    name = name, file = file_path, timestamp = timestamp,
    p_value = append_objs[[length(append_objs)]]$p_value,
    increment = increment, S_prev = S_prev, S_t = S_t,
    signal = signal, crossed = crossed, saturated = saturated,
    h = h, H = H, alpha = alpha, power = power
  )
  state <- if (crossed) "has crossed" else "remains above"
  extra <- if (saturated)
    " The statistic has reached the upper boundary H and is capped there." else ""
  n_obs <- length(append_objs) - as.integer(wrote_config)
  resim_msg <- if (!is.null(p_value_fn) && num_resims > 0 && n_obs > 1)
    sprintf(" (confirmed after %d additional steps)", num_resims) else ""
  failure_message <- sprintf(
    "CUSUM chart '%s' %s the signalling threshold -- S_t = %.4f exceeds h = %.4f.%s%s",
    name, state, S_t, h, resim_msg, extra)
  testthat::expect(
    ok = !signal,
    failure_message = failure_message
  )
  invisible(result)
}