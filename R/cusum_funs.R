z_pval <- function(x, mu_0, sigma) {
  z <- sqrt(length(x)) * (mean(x) - mu_0) / sigma
  pnorm(z, lower.tail = FALSE)
}

compute_LAN_increment <- function(p_val, delta) {
  z_val <- qnorm(p_val, lower.tail = FALSE)
  return(z_val - (delta / 2))
}

LAN_increment_cdf <- function(w, delta, is_H0 = TRUE) {
  if (is_H0) {
    return(pnorm(w + delta / 2))
  } else {
    return(pnorm(w - delta / 2))
  }
}

build_P_matrix <- function(M, H, delta, is_H0 = TRUE) {
  P <- matrix(0, nrow = M + 1, ncol = M + 1)
  w <- H / M
  state_vals <- (0:M) * w
  for (i in 1:(M + 1)) {
    val_i <- state_vals[i]
    P[i, 1] <- LAN_increment_cdf(w/2 - val_i, delta, is_H0)
    if (M > 1) {
      for (j in 2:M) {
        lower_bound <- (j - 1.5) * w - val_i
        upper_bound <- (j - 0.5) * w - val_i
        P[i, j] <- LAN_increment_cdf(upper_bound, delta, is_H0) - 
          LAN_increment_cdf(lower_bound, delta, is_H0)
      }
    }
    
    P[i, M + 1] <- 1 - LAN_increment_cdf((M - 0.5) * w - val_i, delta, is_H0)
  }
  
  P <- P / rowSums(P)
  return(P)
}

compute_stationary_distribution <- function(P_mat) {
  num_states <- nrow(P_mat)
  A <- t(P_mat) - diag(num_states)
  A[num_states, ] <- rep(1, num_states)
  b <- rep(0, num_states)
  b[num_states] <- 1
  pi_steady <- solve(A, b)
  return(pi_steady)
}

estimate_h <- function(target_ARL0, delta, M = 50) {
  target_func <- function(h) {
    log(compute_ARL(h, delta, M, is_H0 = TRUE)) - log(target_ARL0)
  }
  res <- uniroot(target_func, 
                 interval = c(0.1, 20),
                 extendInt = "yes",
                 tol = 1e-4)
  
  return(res$root)
}

compute_ARL <- function(h, delta, M = 50, is_H0 = TRUE) {
  P_mat <- build_P_matrix(M = M, H = h, delta = delta, is_H0 = is_H0)
  Q <- P_mat[1:M, 1:M] 
  I <- diag(M)
  ARL_vec <- tryCatch({
    solve(I - Q, rep(1, M))
  }, error = function(e) {
    return(rep(1e10, M)) 
  })
  return(ARL_vec[1])
}

CUSUM_nonrestarting_plot <- function(pvals, pi_steady, signaling_threshold, 
                                     upper_boundary, delta, num_regions,
                                     main = "CUSUM Chart"){
  M_bins <- length(pi_steady) - 1
  if (num_regions > M_bins) stop("num_regions cannot be larger than the number of states.")
  if (signaling_threshold > upper_boundary) stop("upper_boundary cannot be lower than signaling_threshold.")
  if (signaling_threshold <= 0) stop("signaling_threshold must be positive.")
  
  bin_width <- upper_boundary/M_bins
  palette_func <- colorRampPalette(c("#66cc66", "#fffbc8", "#fc9272", "#de2d26"))
  band_colors <- palette_func(M_bins)
  t_max <- length(pvals)
  plot(0:t_max, c(0, pvals), type = "n", 
       ylim = c(0, upper_boundary),
       xlab = "Time (t)", ylab = expression("CUSUM Statistic (" * S[t] * ")"),
       main = main)
  for (i in 1:M_bins) {
    y_bottom <- (i - 1) * bin_width
    y_top <- i * bin_width
    rect(xleft = -10, ybottom = y_bottom, xright = t_max + 20, ytop = y_top, 
         col = band_colors[i], border = NA)
  }
  region_width <- upper_boundary / num_regions
  state_vals <- (0:M_bins) * bin_width
  for (r in 1:num_regions) {
    y_lower <- (r - 1) * region_width
    y_upper <- r * region_width
    if (r == num_regions) {
      in_region <- which(state_vals >= y_lower & state_vals <= y_upper)
    } else {
      in_region <- which(state_vals >= y_lower & state_vals < y_upper)
    }
    prob_region <- sum(pi_steady[in_region])
    if (prob_region == 0) {
      label_text <- "p = 0"
    } else if (prob_region >= 0.001) {
      label_text <- sprintf("p = %.1f%%", prob_region * 100)
    } else {
      label_text <- sprintf("p = %.1e", prob_region)
    }
    y_mid <- (y_lower + y_upper) / 2
    if (r < num_regions) {
      abline(h = y_upper, col = adjustcolor("black", alpha.f = 0.2), lty = 1, lwd = 1)
    }
    text(x = t_max * 0.98, y = y_mid, labels = label_text,
         col = adjustcolor("black", alpha.f = 0.6), cex = 0.8, font = 2, pos = 2)
  }
  S_t <- numeric(t_max + 1) 
  for (t in 1:t_max) {
    p_t <- pmax(pmin(pvals[t], 1 - 1e-15), 1e-15)
    inc_t <- compute_LAN_increment(p_t, delta)
    S_t[t + 1] <- min(upper_boundary, max(0, S_t[t] + inc_t))
  }
  abline(h = 0, col = 'gray80', lwd = 2, lty = 3)
  abline(h = signaling_threshold, col = 'gray40', lwd = 2, lty = 2)
  abline(h = upper_boundary, col = 'gray40', lwd = 2, lty = 2)
  lines(0:t_max, S_t, col = "black", lwd = 2)
}

.signal_cusum <- function(action, msg, info) {
  cnd <- structure(
    class = c("cusum_signal", action, "condition"),
    list(
      message = if (action == "message") paste0(msg, "\n") else msg,
      call = NULL,
      name = info$name,
      file = info$file,
      p_value = info$p_value,
      increment = info$increment,
      S_prev = info$S_prev,
      S_t = info$S_t,
      h = info$h,
      H = info$H,
      crossed = info$crossed,
      saturated = info$saturated
    )
  )
  switch(action,
         error = stop(cnd),
         warning = warning(cnd),
         message = message(cnd))
}

.save_pval_log <- function(file_path, log) {
  tmp <- tempfile(tmpdir = dirname(file_path), fileext = ".rds.tmp")
  saveRDS(log, file = tmp)
  if (!file.rename(tmp, file_path)) {
    unlink(tmp)
    stop("Failed to save log to '", file_path, "'.", call. = FALSE)
  }
}


#' Read a CUSUM p-value log written by [expect_pval()]
#'
#' @param file_path Path to a `.rds` log file produced by [expect_pval()].
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{`config`}{Named list of chart parameters (`alpha`, `power`, `M`,
#'       `target_ARL0`, `h`, `H`, `ARL0`, `ARL1`, plus `test_name` and 
#'       `created`).}
#'     \item{`data`}{Data frame of recorded observations with columns
#'       `timestamp`, `p_value`, `increment`, `S_t`, `signal`.}
#'   }
#' @export
read_pval_log <- function(file_path) {
  if (!file.exists(file_path))
    stop("Log file not found: ", file_path, call. = FALSE)
  log <- readRDS(file_path)
  if (!is.list(log) || !all(c("config", "data") %in% names(log)))
    stop("'", file_path, "' does not look like an expect_pval() log.",
         call. = FALSE)
  log
}


#' Record a p-value on a running CUSUM log and assert it has not signalled
#'
#' A stateful testthat expectation: each call advances the CUSUM chart for
#' `name` (persisting the update to a `.rds` log) and then asserts, via
#' [testthat::expect()], that the chart is not above its signalling threshold.
#'
#' @param p_value Single numeric value in \[0, 1\] (the p-value to record), or a
#'   function that generates and returns such a value when called.
#' @param name Character scalar naming the test; also the log file name
#'   (a `.rds` extension is appended if absent).
#' @param dir Directory holding the log file. Created if it does not exist.
#'   Defaults to `getOption("cusum.dir", ".")`.
#' @param alpha Target type I error rate. Used to calculate the drift parameter.
#'   Default is `0.05`. On later calls, any supplied value is safely ignored
#'   and read back from the stored config.
#' @param power Target statistical power. Used to calculate the drift parameter.
#'   Default is `0.2`. On later calls, any supplied value is safely ignored
#'   and read back from the stored config.
#' @param target_ARL0 Target in-control average run length, used to calibrate
#'   `h` via `estimate_h()` when `h` is not given. Used only at file creation.
#' @param M Number of CUSUM discretisation bins. Used only at file creation.
#' @param h Optional signalling threshold. If `NULL`, it is calibrated from
#'   `target_ARL0` and the calculated drift parameter. Used only at file creation.
#' @param num_resims Integer. If the threshold `h` is exceeded and `p_value` is
#'   a function, the CUSUM chart continues for this many additional steps
#'   automatically. The expectation fails only if the final state still exceeds
#'   `h`. All intermediate steps are saved to the log. Default is `10`.
#' @param quiet Logical; if `TRUE`, suppress informational messages about
#'   directory/file creation.
#'
#' @return Invisibly, a list describing this call: `name`, `file`,
#'   `timestamp`, `p_value`, `increment`, `S_prev`, `S_t`, `signal`,
#'   `crossed`, `saturated`, and the active chart parameters `h`, `H`,
#'   `alpha`, `power`. As a side effect, registers a testthat expectation that
#'   passes when `S_t < h` and fails (reporting the breach) when the chart has
#'   signalled.
#'
#' @export
expect_pval <- function(p_value,
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
  if (!is.numeric(p_value) || length(p_value) != 1L || !is.finite(p_value))
    stop("The provided or generated p-value must be a single finite numeric value.", call. = FALSE)
  if (p_value < 0 || p_value > 1)
    stop("The provided or generated p-value must lie in [0, 1]; got ", p_value, ".", call. = FALSE)
  if (!is.character(name) || length(name) != 1L || !nzchar(name))
    stop("`name` must be a single non-empty string.", call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be a single numeric value in (0, 1).", call. = FALSE)
  if (!is.numeric(power) || length(power) != 1L || power <= 0 || power >= 1)
    stop("`power` must be a single numeric value in (0, 1).", call. = FALSE)
  if (!is.numeric(num_resims) || length(num_resims) != 1L || num_resims < 0 || num_resims %% 1 != 0)
    stop("`num_resims` must be a single non-negative integer.", call. = FALSE)
  
  file_name <- if (grepl("\\.rds$", name, ignore.case = TRUE))
    name else paste0(name, ".rds")
  file_path <- file.path(dir, file_name)
  
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    if (!quiet) message("expect_pval(): created directory '", dir, "'.")
  }
  
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
      alpha = alpha,
      power = power,
      M = M,
      target_ARL0 = target_ARL0,
      h = h,
      H = H,
      ARL0 = ARL0,
      ARL1 = ARL1
    )
    log <- list(
      config = config,
      data = data.frame(
        timestamp = character(0),
        p_value = numeric(0),
        increment = numeric(0),
        S_t = numeric(0),
        signal = logical(0),
        stringsAsFactors = FALSE
      )
    )
    if (!quiet) message("expect_pval(): created log '", file_path, "'.")
    S_prev <- 0
  } else {
    log <- read_pval_log(file_path)
    config <- log$config
    if (is.null(config$alpha) || is.null(config$power) || is.null(config$h) || is.null(config$H))
      stop("'", file_path, "' is missing CUSUM parameters in its config; ",
           "it does not look like an expect_pval() log.", call. = FALSE)
    
    S_prev <- if (nrow(log$data) > 0) log$data$S_t[nrow(log$data)] else 0
    alpha <- config$alpha
    power <- config$power
    delta <- qnorm(1 - alpha) - qnorm(1 - power)
    h <- config$h
    H <- config$H
  }
  new_rows <- list()
  increment <- compute_LAN_increment(p_value, delta)
  S_t <- min(H, max(0, S_prev + increment))
  signal <- S_t >= h
  crossed <- signal && (S_prev < h)
  saturated <- S_t >= H
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  new_rows[[1]] <- data.frame(
    timestamp = timestamp,
    p_value = p_value,
    increment = increment,
    S_t = S_t,
    signal = signal,
    stringsAsFactors = FALSE
  )
  
  if (signal && !is.null(p_value_fn) && num_resims > 0) {
    for (i in seq_len(num_resims)) {
      S_prev <- S_t
      p_new <- p_value_fn()
      if (!is.numeric(p_new) || length(p_new) != 1L || !is.finite(p_new) || p_new < 0 || p_new > 1) {
        stop("Generated p-value during continuation steps is invalid.", call. = FALSE)
      }
      increment <- compute_LAN_increment(p_new, delta)
      S_t <- min(H, max(0, S_prev + increment))
      signal <- S_t >= h
      crossed <- signal && (S_prev < h)
      saturated <- S_t >= H
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      
      new_rows[[i + 1]] <- data.frame(
        timestamp = timestamp,
        p_value = p_new,
        increment = increment,
        S_t = S_t,
        signal = signal,
        stringsAsFactors = FALSE
      )
    }
  }
  log$data <- rbind(log$data, do.call(rbind, new_rows))
  .save_pval_log(file_path, log)
  
  result <- list(
    name = name,
    file = file_path,
    timestamp = timestamp,
    p_value = new_rows[[length(new_rows)]]$p_value,
    increment = increment,
    S_prev = S_prev,
    S_t = S_t,
    signal = signal,
    crossed = crossed,
    saturated = saturated,
    h = h,
    H = H,
    alpha = alpha,
    power = power
  )
  
  state <- if (crossed){
    "has crossed"
  } else{
    "remains above"
  }
  extra <- if (saturated){
    " The statistic has reached the upper boundary H and is capped there."
  } else {
    ""
  }
  resim_msg <- if (!is.null(p_value_fn) && num_resims > 0 && length(new_rows) > 1){
    sprintf(" (confirmed after %d additional steps)", num_resims)
  } else {
    ""
  }
  failure_message <- sprintf(
    "CUSUM chart '%s' %s the signalling threshold -- S_t = %.4f exceeds h = %.4f.%s%s",
    name, state, S_t, h, resim_msg, extra)
  testthat::expect(
    ok = !signal,
    failure_message = failure_message
  )
  invisible(result)
}


#' Plot a CUSUM p-value log
#'
#' @param file_path Path to a `.rds` log file produced by [expect_pval()].
#' @param num_regions Number of probability regions to display on the plot. 
#'   If `NULL` (default), it is set to `min(10, M)` where `M` is the number 
#'   of discretisation bins.
#' @param main Optional plot title. If `NULL` (default), a title of the form
#'   `"CUSUM Chart: <test_name>"` is built from the name stored in the config
#'   (falling back to the log file's base name).
#'
#' @return Invisibly, the parsed log (as returned by [read_pval_log()]).
#' @export
plot_pval_log <- function(file_path, num_regions = NULL, main = NULL) {
  log <- read_pval_log(file_path)
  cfg <- log$config
  if (nrow(log$data) == 0L)
    stop("Log '", file_path, "' contains no observations yet.", call. = FALSE)
  
  if (is.null(num_regions)) {
    num_regions <- min(10, cfg$M)
  } else if (num_regions > cfg$M) {
    stop("`num_regions` cannot exceed `M` (", cfg$M, ").", call. = FALSE)
  }
  
  if (is.null(main)) {
    test_nm <- cfg$test_name
    if (is.null(test_nm) || !nzchar(test_nm))
      test_nm <- tools::file_path_sans_ext(basename(file_path))
    main <- paste0("CUSUM Chart: ", test_nm)
  }
  delta_val <- qnorm(1 - cfg$alpha) - qnorm(1 - cfg$power)
  P_mat <- build_P_matrix(M = cfg$M, H = cfg$H,
                          delta = delta_val, is_H0 = TRUE)
  pi_steady <- compute_stationary_distribution(P_mat)
  CUSUM_nonrestarting_plot(
    pvals = log$data$p_value,
    pi_steady = pi_steady,
    signaling_threshold = cfg$h,
    upper_boundary = cfg$H,
    delta = delta_val,
    num_regions = num_regions,
    main = main
  )
  invisible(log)
}

#' Group expect_pval() calls under a shared label
#'
#' Mirrors the testthat::test_that() pattern. Pass a description and a
#' `{}` block of expect_pval() calls; the resulting .rds logs all land in
#' `cusum_logs/<description>/`, so related tests stay together
#' on disk and in the alert message.
#'
#' @param description Human-readable label. Sanitised to form the directory name.
#' @param code A `{}` block containing one or more expect_pval() calls.
#' @return Invisibly, a list with each call's return value.
#' @export
cusum_group <- function(description, code) {
  group_name <- gsub("[^A-Za-z0-9_-]+", "_", description)
  group_name <- gsub("^_+|_+$", "", group_name)
  if (!nzchar(group_name))
    stop("`description` produced an empty group name.", call. = FALSE)
  
  group_dir <- file.path("cusum_logs", group_name)
  dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
  
  local_env <- new.env(parent = parent.frame())
  local_env$expect_pval <- function(p_value, name, dir = group_dir, ...) {
    expect_pval(p_value = p_value, name = name, dir = dir, ...)
  }
  
  eval(substitute(code), envir = local_env)
}