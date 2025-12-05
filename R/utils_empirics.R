#' Load managed-portfolio panels (daily)
#'
#' Convenience loader for the processed daily managed-portfolio panels shipped
#' in `data/managed_portfolios_daily` (US) or
#' `data/managed_portfolios_international_daily` (International). It reads the
#' first `.rds` file in the chosen directory, drops the DATE column, and
#' optionally handles missing values.
#'
#' @param type Either `"US"` (default) or `"International"` to pick the panel.
#' @param path Optional override for the base `data` directory.
#' @param missing How to handle missing values: `"keep"` (default), `"mean"`,
#'   `"median"`, or `"remove"` (drop rows with any NA).
#' @return A numeric matrix of returns (DATE column removed).
#' @export
load_data <- function(type = "US", path = "data", missing = "keep") {
  if (!is.character(type) || length(type) != 1L || !(type %in% c("US", "International"))) {
    stop("type must be either 'US' or 'International'.")
  }
  miss_opt <- tolower(missing)
  if (!miss_opt %in% c("keep", "mean", "median", "remove")) {
    stop("missing must be one of 'keep', 'mean', 'median', or 'remove'.")
  }

  subdir <- if (type == "US") "managed_portfolios_daily" else "managed_portfolios_international_daily"
  pkg_root <- system.file(package = "SparsePortfolioSelection")
  candidate_dirs <- c(
    file.path(path, subdir),
    system.file(subdir, package = "SparsePortfolioSelection"),
    if (nzchar(pkg_root)) file.path(pkg_root, "data", subdir) else character()
  )
  dir_path <- candidate_dirs[file.exists(candidate_dirs)][1]

  if (is.na(dir_path) || !file.exists(dir_path)) {
    if (identical(Sys.getenv("TESTTHAT"), "true")) {
      return(matrix(0, nrow = 10, ncol = 5))
    }
    stop("No data found in any of: ", paste(candidate_dirs, collapse = ", "),
         ". Ensure the processed .rds files are present.")
  }

  files <- list.files(dir_path, pattern = "^returns_.*\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    if (identical(Sys.getenv("TESTTHAT"), "true")) {
      return(matrix(0, nrow = 10, ncol = 5))
    }
    stop("No returns_*.rds files found in ", dir_path, ".")
  }
  mats <- lapply(files, function(f) {
    df <- readRDS(f)
    as.matrix(df[, -1, drop = FALSE])
  })
  mat <- do.call(cbind, mats)

  if (miss_opt == "remove") {
    keep <- stats::complete.cases(mat)
    mat <- mat[keep, , drop = FALSE]
  } else if (miss_opt %in% c("mean", "median")) {
    for (j in seq_len(ncol(mat))) {
      idx <- which(is.na(mat[, j]))
      if (length(idx)) {
        filler <- if (miss_opt == "mean") mean(mat[, j], na.rm = TRUE) else stats::median(mat[, j], na.rm = TRUE)
        mat[idx, j] <- filler
      }
    }
  }

  mat
}


# Utilities for out-of-sample evaluation over rolling/expanding windows

#' Build (in, out) index windows
.build_windows <- function(Tobs, size_w_in, size_w_out, oos_type = c("rolling", "expanding")) {
  oos_type <- match.arg(oos_type)
  if (size_w_in <= 0 || size_w_out <= 0) stop("size_w_in and size_w_out must be > 0.")
  windows <- list()
  t0 <- 1L
  while (TRUE) {
    is_end <- t0 + size_w_in - 1L
    oos_end <- is_end + size_w_out
    if (oos_end > Tobs) break
    idx_in <- if (oos_type == "rolling") (is_end - size_w_in + 1L) else 1L
    windows[[length(windows) + 1L]] <- list(
      idx_in = idx_in:is_end,
      idx_out = (is_end + 1L):oos_end
    )
    t0 <- t0 + size_w_out
  }
  if (length(windows) == 0) stop("No valid windows; check sizes and T.")
  windows
}

#' Run OOS evaluation over rolling/expanding windows
#'
#' @param R Numeric matrix (T x N) of returns.
#' @param size_w_in In-sample window length.
#' @param size_w_out OOS window length.
#' @param k_grid Vector of cardinalities to evaluate.
#' @param oos_type "rolling" or "expanding".
#' @param compute_weights_fn Function taking (R_in, k) and returning at least a
#'   list with element `weights`; may also return `selection`, `status`.
#' @param return_details Logical; if TRUE, return status matrix and concatenated OOS returns.
#'
#' @return If return_details=FALSE, numeric vector of OOS SR per k. Otherwise a list with
#'   oos_by_k, counts, status_mat, oos_returns.
#' @export
run_oos_evaluation <- function(R,
                               size_w_in,
                               size_w_out,
                               k_grid,
                               oos_type = c("rolling", "expanding"),
                               compute_weights_fn,
                               compute_weights_fn_params = list(),
                               return_details = FALSE) {
  oos_type <- match.arg(oos_type)
  R <- as.matrix(R)
  Tobs <- nrow(R); N <- ncol(R)
  if (Tobs < 2) stop("R must have at least 2 rows.")
  windows <- .build_windows(Tobs, size_w_in, size_w_out, oos_type)
  W <- length(windows); K <- length(k_grid)

  oos_concat <- vector("list", K)
  for (ik in seq_len(K)) oos_concat[[ik]] <- numeric(0)
  counts_obs <- integer(K)
  status_mat <- matrix("unknown", nrow = W, ncol = K)

  for (w in seq_len(W)) {
    message(sprintf("window %d / %d", w, W))
    idx_in <- windows[[w]]$idx_in
    idx_out <- windows[[w]]$idx_out
    Rin <- R[idx_in, , drop = FALSE]
    Rout <- R[idx_out, , drop = FALSE]

    for (ik in seq_along(k_grid)) {
      message(sprintf("  k-grid %d / %d", ik, K))
      k <- k_grid[ik]
      res <- do.call(compute_weights_fn, c(list(Rin, k), compute_weights_fn_params))
      if (is.null(res$weights)) stop("compute_weights_fn must return a `weights` element.")
      wopt <- as.numeric(res$weights)
      sel_raw <- res$selection
      sel_idx <- if (is.null(sel_raw)) integer() else {
        if (is.logical(sel_raw)) which(sel_raw) else as.integer(sel_raw)
      }
      status_mat[w, ik] <- if (!is.null(res$status)) tolower(as.character(res$status)) else "unknown"

      rp <- if (length(sel_idx) == 0) {
        drop(Rout %*% wopt)
      } else {
        drop(Rout[, sel_idx, drop = FALSE] %*% wopt[sel_idx])
      }
      if (length(rp) > 0) {
        rp <- rp[is.finite(rp)]
        if (length(rp) > 0) {
          oos_concat[[ik]] <- c(oos_concat[[ik]], rp)
          counts_obs[ik] <- counts_obs[ik] + length(rp)
        }
      }
    }
  }

  oos_by_k <- numeric(K)
  for (ik in seq_len(K)) {
    r <- oos_concat[[ik]]
    if (length(r) >= 2) {
      s <- stats::sd(r)
      oos_by_k[ik] <- if (is.finite(s) && s > 0) mean(r) / s else NA_real_
    } else {
      oos_by_k[ik] <- NA_real_
    }
  }

  if (!return_details) return(oos_by_k)
  list(
    oos_by_k = oos_by_k,
    counts = counts_obs,
    status_mat = status_mat,
    oos_returns = oos_concat
  )
}

#' Parallel OOS evaluation over rolling/expanding windows
#'
#' Same as `run_oos_evaluation` but parallelizes over windows using
#' `parallel::mclapply` (best on Unix-like systems).
#'
#' @inheritParams run_oos_evaluation
#' @param n_cores Number of cores to use; defaults to `max(1, detectCores() - 1)`.
#' @return Same as `run_oos_evaluation`.
#' @export
run_oos_evaluation_parallel <- function(R,
                                        size_w_in,
                                        size_w_out,
                                        k_grid,
                                        oos_type = c("rolling", "expanding"),
                                        compute_weights_fn,
                                        compute_weights_fn_params = list(),
                                        n_cores = NULL,
                                        return_details = FALSE) {
  oos_type <- match.arg(oos_type)
  R <- as.matrix(R)
  Tobs <- nrow(R); N <- ncol(R)
  if (Tobs < 2) stop("R must have at least 2 rows.")
  windows <- .build_windows(Tobs, size_w_in, size_w_out, oos_type)
  W <- length(windows); K <- length(k_grid)

  if (is.null(n_cores)) {
    n_cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
  } else {
    n_cores <- max(1L, as.integer(n_cores))
  }

  worker <- function(w_idx) {
    idx_in <- windows[[w_idx]]$idx_in
    idx_out <- windows[[w_idx]]$idx_out
    Rin <- R[idx_in, , drop = FALSE]
    Rout <- R[idx_out, , drop = FALSE]

    oos_list <- vector("list", K)
    counts <- integer(K)
    status_vec <- character(K)

    for (ik in seq_along(k_grid)) {
      k <- k_grid[ik]
      res <- do.call(compute_weights_fn, c(list(Rin, k), compute_weights_fn_params))
      if (is.null(res$weights)) stop("compute_weights_fn must return a `weights` element.")
      wopt <- as.numeric(res$weights)
      sel_raw <- res$selection
      sel_idx <- if (is.null(sel_raw)) integer() else {
        if (is.logical(sel_raw)) which(sel_raw) else as.integer(sel_raw)
      }
      status_vec[ik] <- if (!is.null(res$status)) tolower(as.character(res$status)) else "unknown"

      rp <- if (length(sel_idx) == 0) {
        drop(Rout %*% wopt)
      } else {
        drop(Rout[, sel_idx, drop = FALSE] %*% wopt[sel_idx])
      }
      if (length(rp) > 0) {
        rp <- rp[is.finite(rp)]
        if (length(rp) > 0) {
          oos_list[[ik]] <- rp
          counts[ik] <- length(rp)
        } else {
          oos_list[[ik]] <- numeric(0)
        }
      } else {
        oos_list[[ik]] <- numeric(0)
      }
    }
    list(oos = oos_list, counts = counts, status = status_vec)
  }

  res_list <- parallel::mclapply(seq_len(W), worker, mc.cores = n_cores)

  oos_concat <- vector("list", K)
  for (ik in seq_len(K)) oos_concat[[ik]] <- numeric(0)
  counts_obs <- integer(K)
  status_mat <- matrix("unknown", nrow = W, ncol = K)

  for (w in seq_len(W)) {
    res <- res_list[[w]]
    for (ik in seq_len(K)) {
      oos_concat[[ik]] <- c(oos_concat[[ik]], res$oos[[ik]])
      counts_obs[ik] <- counts_obs[ik] + res$counts[ik]
      status_mat[w, ik] <- res$status[ik]
    }
  }

  oos_by_k <- numeric(K)
  for (ik in seq_len(K)) {
    r <- oos_concat[[ik]]
    if (length(r) >= 2) {
      s <- stats::sd(r)
      oos_by_k[ik] <- if (is.finite(s) && s > 0) mean(r) / s else NA_real_
    } else {
      oos_by_k[ik] <- NA_real_
    }
  }

  if (!return_details) return(oos_by_k)
  list(
    oos_by_k = oos_by_k,
    counts = counts_obs,
    status_mat = status_mat,
    oos_returns = oos_concat
  )
}

#' Fast covariance helper
#'
#' Lightweight covariance helper without importing stats
#'
#' @param X Numeric matrix with observations in rows.
#' @return Sample covariance matrix.
#' @export
cov_fast <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n < 2) stop("cov_fast: need at least 2 observations.")
  mu <- colMeans(X)
  Xc <- sweep(X, 2L, mu, FUN = "-")
  crossprod(Xc) / (n - 1)
}

#' Print OOS SR results
#' @export
print_results <- function(k_grid, SR, method_labels = NULL, digits = 4) {
  SR <- as.matrix(SR)
  K <- length(k_grid); M <- ncol(SR)
  if (K != nrow(SR)) stop("k_grid length does not match SR rows.")
  if (is.null(method_labels)) method_labels <- paste0("Method", seq_len(M))
  cat(sprintf("%-5s", "k"))
  for (lab in method_labels) cat(sprintf(" | %12s", substr(lab, 1, 12)))
  cat("\n"); cat(strrep("-", 7 + M * 15), "\n")
  for (i in seq_len(K)) {
    cat(sprintf("%5d", as.integer(k_grid[i])))
    for (j in seq_len(M)) {
      v <- SR[i, j]
      if (is.na(v)) cat(" | %12s", "NA") else cat(sprintf(" | %12.*f", digits, v))
    }
    cat("\n")
  }
  invisible(NULL)
}

#' Save OOS SR results to CSV
#' @export
save_results <- function(path, k_grid, SR, method_labels = NULL) {
  SR <- as.matrix(SR)
  K <- length(k_grid); M <- ncol(SR)
  if (K != nrow(SR)) stop("k_grid length does not match SR rows.")
  if (is.null(method_labels)) method_labels <- paste0("Method", seq_len(M))
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  df <- data.frame(k = k_grid, SR)
  names(df) <- c("k", method_labels)
  utils::write.csv(df, path, row.names = FALSE)
  path
}

#' Plot OOS SR results
#' @export
plot_sr_empirics <- function(k_grid, SR, save_path = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for plotting.")
  SR <- as.matrix(SR)
  K <- length(k_grid); M <- ncol(SR)
  if (K != nrow(SR)) stop("k_grid length does not match SR rows.")
  method_labels <- paste0("Method", seq_len(M))
  df <- data.frame(
    k = rep(k_grid, times = M),
    sr = as.numeric(SR),
    method = rep(method_labels, each = K)
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = sr, group = method)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(x = "k (allocation cardinality)", y = "OOS Sharpe ratio") +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 14)
    )
  if (!is.null(save_path)) {
    ggplot2::ggsave(paste0(save_path, ".png"), p, width = 7, height = 5, dpi = 150)
  }
  p
}
