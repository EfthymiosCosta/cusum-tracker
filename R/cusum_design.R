.G <- function(t, xi, c) {
  1 - pnorm(c / sqrt(t) - xi * sqrt(t)) + exp(2 * xi * c) * pnorm(-c / sqrt(t) - xi * sqrt(t))
}

.gap_from_pleave <- function(amax, p_leave, delta, rho = 0.5826) {
  f <- function(g) .G(amax, -delta / 2, g + rho) - p_leave
  if (f(0) <= 0) return(0)
  uniroot(f, c(0, sqrt(amax)), extendInt = "downX", tol = 1e-10)$root
}

.make_mesh <- function(H, n_bins) {
  grid <- seq(0, H, length.out = n_bins)
  dx <- grid[2] - grid[1]
  lo <- grid[2:(n_bins - 1)] - dx / 2
  lo[1] <- 0
  hi <- grid[2:(n_bins - 1)] + dx / 2
  hi[length(hi)] <- H
  list(grid = grid, lo = lo, hi = hi, H = H, ns = n_bins)
}

.refine_mesh <- function(mesh, ts, min_width = 1e-8) {
  H <- mesh$H
  orig <- mesh$grid[2:(mesh$ns - 1)]
  ts <- ts[is.finite(ts) & ts > min_width & ts < H - min_width]
  knots <- sort(unique(c(mesh$lo[1], mesh$hi, ts)))
  knots <- knots[c(TRUE, diff(knots) > min_width)]
  lo <- head(knots, -1)
  hi <- tail(knots, -1)
  pts <- vapply(seq_along(lo), function(j) {
    cand <- orig[orig > lo[j] & orig < hi[j]]
    if (length(cand) == 1L) cand else (lo[j] + hi[j]) / 2
  }, numeric(1))
  list(grid = c(0, pts, H), lo = lo, hi = hi, H = H, ns = length(pts) + 2L)
}

.step_kernel <- function(mesh, mu) {
  ns <- mesh$ns
  m <- mesh$grid + mu
  P <- matrix(0, ns, ns)
  P[, 1] <- pnorm(-m)
  P[, ns] <- 1 - pnorm(mesh$H - m)
  P[, 2:(ns - 1)] <- pnorm(outer(-m, mesh$hi, "+")) - pnorm(outer(-m, mesh$lo, "+"))
  P <- pmax(P, 0)
  P / rowSums(P)
}

.split_kernel <- function(G, mesh, t, mu) {
  ns <- mesh$ns
  t <- min(max(t, 0), mesh$H)
  if (t >= mesh$H) return(list(above = matrix(0, ns, ns), below = G))
  if (t <= 0) return(list(above = G, below = matrix(0, ns, ns)))
  i <- findInterval(t, c(mesh$lo[1], mesh$hi), left.open = TRUE)
  j <- i + 1L
  above <- G
  above[, seq_len(j - 1L)] <- 0
  above[, j] <- pnorm(mesh$hi[i] - mesh$grid - mu) - pnorm(t - mesh$grid - mu)
  list(above = above, below = G - above)
}

.jump_kernel <- function(mesh, v, amax, mu) {
  G <- .step_kernel(mesh, mu)
  sp <- .split_kernel(G, mesh, v, mu)
  ns <- mesh$ns
  Mn <- diag(ns)
  K <- matrix(0, ns, ns)
  cost <- numeric(ns)
  for (n in seq_len(amax) - 1L) {
    S <- if (n < amax - 1L) Mn %*% sp$below else Mn %*% G
    K <- K + S
    cost <- cost + (n + 1) * rowSums(S)
    if (n < amax - 1L) Mn <- Mn %*% sp$above
  }
  list(K = K, cost = cost)
}

.jump_evaluate <- function(v, h, gap, amax, delta, n_bins) {
  mesh <- .refine_mesh(.make_mesh(h + gap, n_bins), v)
  ns <- mesh$ns
  tr <- which(mesh$grid < h - 1e-9 * max(1, h))
  k0 <- .jump_kernel(mesh, v, amax, -delta / 2)
  k1 <- .jump_kernel(mesh, v, amax, delta / 2)
  A <- t(k0$K) - diag(ns)
  A[ns, ] <- 1
  pi0 <- as.numeric(solve(qr(A, LAPACK = TRUE), c(rep(0, ns - 1), 1)))
  pi0 <- pi0 / sum(pi0)
  u0 <- solve(qr(diag(length(tr)) - k0$K[tr, tr, drop = FALSE], LAPACK = TRUE), rep(1, length(tr)))[1]
  u1 <- solve(qr(diag(length(tr)) - k1$K[tr, tr, drop = FALSE], LAPACK = TRUE), rep(1, length(tr)))[1]
  list(h = h, v = v, u0 = u0, u1 = u1, cbar0 = sum(k0$cost * pi0))
}

.solve_h <- function(u0_at, U0, start) {
  g <- function(h) {
    u <- u0_at(h)
    if (!is.finite(u) || u <= 0) Inf else log(u) - log(U0)
  }
  lo <- start
  g_lo <- g(lo)
  while (g_lo > 0) {
    lo <- lo / 2
    g_lo <- g(lo)
  }
  hi <- lo
  g_hi <- g_lo
  while (g_hi < 0) {
    lo <- hi
    g_lo <- g_hi
    hi <- 1.5 * hi
    g_hi <- g(hi)
  }
  if (!is.finite(g_hi))
    stop("u0 could not be computed while bracketing h; try a smaller U0 or a larger n_bins.", call. = FALSE)
  if (g_hi == 0) return(hi)
  uniroot(g, c(lo, hi), f.lower = g_lo, f.upper = g_hi, tol = 1e-8)$root
}

cusum_design <- function(alpha = 0.05, power = 0.2, U0 = 1e5, cbar0max = 1.1,
                         amax = 10, p_leave = 0.01, n_bins = 151, n_scan = 30,
                         tol_h = 1e-3, max_h_iter = 40, cores = NULL) {
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  if (!is.numeric(power) || length(power) != 1L || power <= alpha || power >= 1)
    stop("`power` must be a single number in (alpha, 1).", call. = FALSE)
  if (!is.numeric(U0) || length(U0) != 1L || U0 <= 1)
    stop("`U0` must be a single number greater than 1.", call. = FALSE)
  if (!is.numeric(cbar0max) || length(cbar0max) != 1L || cbar0max <= 1)
    stop("`cbar0max` must be a single number greater than 1.", call. = FALSE)
  if (!is.numeric(amax) || length(amax) != 1L || amax < 2 || amax %% 1 != 0)
    stop("`amax` must be a single integer of at least 2.", call. = FALSE)
  if (!is.numeric(p_leave) || length(p_leave) != 1L || p_leave <= 0 || p_leave >= 1)
    stop("`p_leave` must be a single number in (0, 1).", call. = FALSE)
  if (!is.numeric(n_bins) || length(n_bins) != 1L || n_bins < 3 || n_bins %% 1 != 0)
    stop("`n_bins` must be a single integer of at least 3.", call. = FALSE)
  if (!is.numeric(n_scan) || length(n_scan) != 1L || n_scan < 2 || n_scan %% 1 != 0)
    stop("`n_scan` must be a single integer of at least 2.", call. = FALSE)
  amax <- as.integer(amax)
  n_bins <- as.integer(n_bins)
  if (is.null(cores))
    cores <- if (.Platform$OS.type == "windows") 1L else max(1L, parallel::detectCores(), na.rm = TRUE)
  delta <- qnorm(1 - alpha) - qnorm(1 - power)
  gap <- .gap_from_pleave(amax, p_leave, delta)
  v_lo <- 1e-8
  evaluate <- function(v, h) .jump_evaluate(v, h, gap, amax, delta, n_bins)
  h_lo <- .solve_h(function(h) evaluate(h * (1 - 1e-8), h)$u0, U0, start = 1)
  h_up <- .solve_h(function(h) evaluate(v_lo, h)$u0, U0, start = h_lo)
  infeasible <- list(feasible = FALSE, u1 = Inf)
  profile_at <- function(h) {
    out <- tryCatch({
      v_hi <- h * (1 - 1e-8)
      e_hi <- evaluate(v_hi, h)
      if (!all(is.finite(c(e_hi$u0, e_hi$u1, e_hi$cbar0))) || e_hi$cbar0 > cbar0max || e_hi$u0 < U0) {
        infeasible
      } else {
        e_lo <- evaluate(v_lo, h)
        vc <- if (e_lo$cbar0 <= cbar0max) v_lo else
          uniroot(function(v) cbar0max - evaluate(v, h)$cbar0, c(v_lo, v_hi),
                  f.lower = cbar0max - e_lo$cbar0, f.upper = cbar0max - e_hi$cbar0, tol = 1e-8)$root
        vu <- if (e_lo$u0 >= U0) v_lo else
          uniroot(function(v) evaluate(v, h)$u0 - U0, c(v_lo, v_hi),
                  f.lower = e_lo$u0 - U0, f.upper = e_hi$u0 - U0, tol = 1e-8)$root
        e <- evaluate(max(vc, vu) * (1 + 1e-9), h)
        e$feasible <- e$cbar0 <= cbar0max * (1 + 1e-6) && e$u0 >= U0 * (1 - 1e-6)
        e
      }
    }, error = function(err) infeasible)
    out$h <- h
    out
  }
  key <- function(h) sprintf("%.10f", h)
  scanned <- parallel::mclapply(seq(h_lo, h_up, length.out = n_scan), profile_at, mc.cores = cores)
  cache <- list()
  for (p in scanned) cache[[key(p$h)]] <- p
  at <- function(h) {
    k <- key(h)
    if (is.null(cache[[k]])) cache[[k]] <<- profile_at(h)
    cache[[k]]
  }
  feas <- Filter(function(p) isTRUE(p$feasible), cache)
  if (!length(feas))
    stop("No feasible design found for h in [", signif(h_lo, 4), ", ", signif(h_up, 4), "].", call. = FALSE)
  cur <- feas[[which.min(vapply(feas, function(p) p$u1, numeric(1)))]]
  step <- (h_up - h_lo) / (2 * n_scan)
  dir <- -1
  for (it in seq_len(max_h_iter)) {
    moved <- FALSE
    for (d in unique(c(dir, -dir))) {
      hn <- cur$h + d * step
      if (hn < h_lo || hn > h_up) next
      p <- at(hn)
      if (isTRUE(p$feasible) && p$u1 < cur$u1) {
        cur <- p
        dir <- d
        moved <- TRUE
        break
      }
    }
    if (!moved) {
      step <- step / 2
      if (step < tol_h) break
    }
  }
  feas <- Filter(function(p) isTRUE(p$feasible), cache)
  best <- feas[[which.min(vapply(feas, function(p) p$u1, numeric(1)))]]
  list(alpha = alpha, power = power, U0 = U0, cbar0max = cbar0max, amax = amax,
       p_leave = p_leave, n_bins = n_bins, delta = delta, gap = gap,
       h = best$h, H = best$h + gap, v = best$v,
       u0 = best$u0, u1 = best$u1, cbar0 = best$cbar0,
       h_lo = h_lo, h_up = h_up)
}

.design_key <- function(alpha, power, U0, cbar0max, amax, p_leave, n_bins) {
  sprintf("design_alpha%g_power%g_U0%g_cbar0max%g_amax%d_pleave%g_nbins%d",
          alpha, power, U0, cbar0max, as.integer(amax), p_leave, as.integer(n_bins))
}

cusum_design_cached <- function(alpha, power, U0, cbar0max, amax, p_leave,
                                n_bins, n_scan, design_dir, quiet = FALSE) {
  path <- file.path(design_dir, paste0(.design_key(alpha, power, U0, cbar0max, amax, p_leave, n_bins), ".json"))
  if (file.exists(path)) return(jsonlite::read_json(path, simplifyVector = TRUE))
  if (!quiet) message("cusum_design(): no cached design, running the search.")
  design <- cusum_design(alpha = alpha, power = power, U0 = U0, cbar0max = cbar0max,
                         amax = amax, p_leave = p_leave, n_bins = n_bins, n_scan = n_scan)
  dir.create(design_dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = design_dir, fileext = ".json.tmp")
  jsonlite::write_json(design, tmp, auto_unbox = TRUE, digits = NA, pretty = TRUE)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Failed to write the design cache to '", path, "'.", call. = FALSE)
  }
  if (!quiet) message("cusum_design(): cached design in '", path, "'.")
  design
}
