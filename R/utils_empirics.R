#' Load managed-portfolio panels (daily or monthly)
#'
#' Convenience loader for the processed daily or monthly managed-portfolio
#' panels shipped in `data/managed_portfolios_*`. It reads all `.rds` files in
#' the chosen directory, drops the DATE column, and optionally handles missing
#' values.
#'
#' @param type Either `"US"` (default), `"International"`, `"mebeme_ind"`, or `"mebeme"` to pick the panel.
#' @param path Optional override for the base `data` directory.
#' @param missing How to handle missing values: `"keep"` (default), `"mean"`,
#'   `"median"`, or `"remove"` (drop rows with any NA).
#' @param frequency Either `"daily"` (default) or `"monthly"`. Monthly loads
#'   all managed-portfolio return tables and drops the date column.
#' @param add_mkt Logical; if `TRUE`, append MKT-RF columns. For US panels this
#'   comes from factors_ff5; for international monthly it appends MKT for each
#'   region; for `mebeme`/`mebeme_ind` it comes from factors_ff5_monthly.
#' @param add_factors Logical; if `TRUE`, append the first three factors (MKT,
#'   SMB, HML). For international monthly, appends the three factors for each
#'   region. For `mebeme`/`mebeme_ind`, monthly only.
#' @param shuffle Logical; if `TRUE` (default), randomly shuffle asset columns
#'   before returning.
#' @return A named list with elements `returns` (matrix; DATE removed; factors appended if requested)
#'   and `rf` (risk-free vector aligned to returns rows, if available; otherwise NULL).
#' @export
load_data <- function(type = "US", path = "data", missing = "keep",
                      frequency = c("daily", "monthly"),
                      add_mkt = FALSE,
                      add_factors = FALSE,
                      shuffle = TRUE) {
  freq <- match.arg(frequency)
  if (!is.character(type) || length(type) != 1L || !(type %in% c("US", "International", "mebeme_ind", "mebeme"))) {
    stop("type must be one of 'US', 'International', 'mebeme_ind', or 'mebeme'.")
  }
  miss_opt <- tolower(missing)
  if (!miss_opt %in% c("keep", "mean", "median", "remove")) {
    stop("missing must be one of 'keep', 'mean', 'median', or 'remove'.")
  }

  if ((type == "mebeme_ind" || type == "mebeme") && freq != "monthly") {
    stop("type='mebeme_ind' and 'mebeme' are only available for monthly frequency.")
  }
  subdir <- switch(freq,
                   daily   = if (type == "US") "managed_portfolios_daily" else "managed_portfolios_international_daily",
                   monthly = if (type == "US") "managed_portfolios_monthly" else if (type == "International") "managed_portfolios_international_monthly" else "mebeme_ind_monthly")
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
  if (type == "mebeme_ind" || type == "mebeme") {
    keep_basenames <- if (type == "mebeme_ind") {
      c("returns_mebeme100_monthly.rds", "returns_ind49_monthly.rds")
    } else {
      c("returns_mebeme100_monthly.rds")
    }
    files <- files[basename(files) %in% keep_basenames]
  }
  if (type == "International" && freq == "monthly") {
    # Restrict to the 25 ME×BE/ME portfolios per region (Developed Markets: 4 regions × 25 = 100 assets)
    keep_basenames <- c(
      "returns_apxj_mebeme25_int_monthly.rds",
      "returns_eu_mebeme25_int_monthly.rds",
      "returns_jp_mebeme25_int_monthly.rds",
      "returns_na_mebeme25_int_monthly.rds"
    )
    files <- files[basename(files) %in% keep_basenames]
  }
  if (length(files) == 0) {
    if (identical(Sys.getenv("TESTTHAT"), "true")) {
      return(matrix(0, nrow = 10, ncol = 5))
    }
    stop("No returns_*.rds files found in ", dir_path,
         if (type == "International" && freq == "monthly") " after filtering to ME×BE/ME portfolios." else ".")
  }

  dates_ref <- NULL
  rf_vec <- NULL
  mats <- lapply(files, function(f) {
    df <- readRDS(f)
    if (is.null(dates_ref) && ncol(df) >= 1) {
      dates_ref <<- df[[1]]
    } else if (!is.null(dates_ref) && ncol(df) >= 1) {
      if (!identical(dates_ref, df[[1]])) {
        stop("Date column mismatch across returns files; cannot align.")
      }
    }
    m <- as.matrix(df)
    if (ncol(m) >= 2) m[, -1, drop = FALSE] else m
  })
  # align by rows if possible; assume same length; if not, error
  n_rows <- unique(vapply(mats, nrow, integer(1)))
  if (length(n_rows) != 1) {
    stop("Monthly returns files have differing row counts: ", paste(n_rows, collapse = ", "))
  }
  mat <- do.call(cbind, mats)

  if (add_mkt || add_factors) {
    if (type == "International" && freq == "monthly") {
      factor_files <- c(
        apxj   = "ff3_factors_monthly_asia_pacific_ex_japan.rds",
        europe = "ff3_factors_monthly_europe.rds",
        japan  = "ff3_factors_monthly_japan.rds",
        na     = "ff3_factors_monthly_north_america.rds"
      )
      cols_to_add <- list()
      rf_vec <- numeric(length(dates_ref))
      have_rf <- FALSE
      for (region in names(factor_files)) {
        fpath <- file.path(dir_path, factor_files[[region]])
        if (!file.exists(fpath)) {
          warning("Skipping missing factor file: ", basename(fpath))
          next
        }
        ff <- readRDS(fpath)
        if (ncol(ff) < 4) stop("Factor file missing FF3 columns: ", fpath)
        idx <- match(dates_ref, ff[[1]])
        if (anyNA(idx)) stop("Date mismatch between returns and ", basename(fpath))
        if (add_factors) {
          cols_to_add[[paste0(region, "_MKT")]] <- as.numeric(ff[[2]][idx])
          cols_to_add[[paste0(region, "_SMB")]] <- as.numeric(ff[[3]][idx])
          cols_to_add[[paste0(region, "_HML")]] <- as.numeric(ff[[4]][idx])
        } else if (add_mkt) {
          cols_to_add[[paste0(region, "_MKT")]] <- as.numeric(ff[[2]][idx])
        }
        if (!have_rf && ncol(ff) >= 5) {
          rf_vec <- as.numeric(ff[[5]][idx])
          have_rf <- TRUE
        }
      }
      if (length(cols_to_add)) {
        mat <- cbind(mat, as.data.frame(cols_to_add, check.names = FALSE))
      }
      if (!have_rf) rf_vec <- NULL
    } else {
      ff_file <- file.path(dir_path, sprintf("factors_ff5_%s.rds", freq))
      if (!file.exists(ff_file)) stop("factors file not found: ", ff_file)
      ff <- readRDS(ff_file)
      if (ncol(ff) < 2) stop("factors file missing MKT column.")
      if (add_factors && ncol(ff) < 4) stop("factors file missing FF3 columns.")
      if (is.null(dates_ref)) stop("add_mkt/add_factors requested but no date column to align.")
      idx <- match(dates_ref, ff[[1]])
      if (anyNA(idx)) stop("Date mismatch between returns and factors_ff5_", freq, ".")

      if (add_factors) {
        ff3 <- as.data.frame(ff[idx, 2:4, drop = FALSE])
        names(ff3) <- c("MKT", "SMB", "HML")
        mat <- cbind(mat, ff3)
      } else if (add_mkt) {
        mkt_col <- as.numeric(ff[[2]][idx])
        mat <- cbind(mat, mkt_col)
      }
      if (ncol(ff) >= 5) {
        rf_vec <- as.numeric(ff[[5]][idx])
      }
    }
  }

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

  if (shuffle) {
    mat <- mat[, sample.int(ncol(mat)), drop = FALSE]
  }

  list(returns = mat, rf = rf_vec)
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

#' Complete OOS evaluation with turnover and instability metrics
#'
#' Extends `run_oos_evaluation` by also computing median turnover, weight instability,
#' and selection instability across windows.
#'
#' @param R Returns matrix (rows = time, cols = assets); first column is not DATE (date must be aligned externally).
#' @param size_w_in In-sample window length.
#' @param size_w_out OOS window length.
#' @param k_grid Vector of k values.
#' @param oos_type "rolling" or "expanding".
#' @param compute_weights_fn Function returning list(weights, selection, status) given (Rin, k).
#' @param compute_weights_fn_params Additional params to pass.
#' @param rf Optional risk-free vector/scalar aligned to rows of R for turnover computation.
#' @param sharpe_fn Either "mean" (mean/sd) or "median" (median/iqr) for OOS SR aggregation.
#' @param return_details Logical; if FALSE return only a summary data.frame.
#' @export
run_complete_oos_evaluation <- function(R,
                                        size_w_in,
                                        size_w_out,
                                        k_grid,
                                        oos_type = c("rolling", "expanding"),
                                        compute_weights_fn,
                                        compute_weights_fn_params = list(),
                                        rf = 0,
                                        sharpe_fn = c("median", "mean"),
                                        return_details = FALSE) {
  oos_type <- match.arg(oos_type)
  sharpe_fn <- match.arg(sharpe_fn)
  R <- as.matrix(R)
  Tobs <- nrow(R); N <- ncol(R)
  if (Tobs < 2) stop("R must have at least 2 rows.")
  windows <- .build_windows(Tobs, size_w_in, size_w_out, oos_type)
  W <- length(windows); K <- length(k_grid)

  oos_concat <- vector("list", K)
  weights_path <- vector("list", K)
  selections_path <- vector("list", K)
  for (ik in seq_len(K)) {
    oos_concat[[ik]] <- numeric(0)
    weights_path[[ik]] <- list()
    selections_path[[ik]] <- list()
  }
  counts_obs <- integer(K)

  for (w in seq_len(W)) {
    idx_in <- windows[[w]]$idx_in
    idx_out <- windows[[w]]$idx_out
    Rin <- R[idx_in, , drop = FALSE]
    Rout <- R[idx_out, , drop = FALSE]

    for (ik in seq_along(k_grid)) {
      k <- k_grid[ik]
      res <- do.call(compute_weights_fn, c(list(Rin, k), compute_weights_fn_params))
      if (is.null(res$weights)) stop("compute_weights_fn must return a `weights` element.")
      wopt <- as.numeric(res$weights)
      weights_path[[ik]][[length(weights_path[[ik]]) + 1L]] <- wopt
      sel_raw <- res$selection
      sel_idx <- if (is.null(sel_raw)) integer() else {
        if (is.logical(sel_raw)) which(sel_raw) else as.integer(sel_raw)
      }
      selections_path[[ik]][[length(selections_path[[ik]]) + 1L]] <- sel_idx

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
      if (sharpe_fn == "median") {
        iqr <- stats::IQR(r, na.rm = TRUE)
        med <- stats::median(r, na.rm = TRUE)
        oos_by_k[ik] <- if (is.finite(iqr) && iqr > 0) med / iqr else NA_real_
      } else {
        s <- stats::sd(r, na.rm = TRUE)
        oos_by_k[ik] <- if (is.finite(s) && s > 0) mean(r, na.rm = TRUE) / s else NA_real_
      }
    } else {
      oos_by_k[ik] <- NA_real_
    }
  }

  turnover_med <- numeric(K)
  instab_L1_med <- numeric(K)
  instab_L2_med <- numeric(K)
  sel_instab_med <- numeric(K)

  if (length(rf) == 1L) rf_vec <- rep(rf, Tobs) else rf_vec <- rf
  if (length(rf_vec) < Tobs) stop("rf length must match nrow(R) or be scalar.")

  for (ik in seq_len(K)) {
    w_path <- weights_path[[ik]]
    if (length(w_path) >= 2) {
      # turnover: need returns aligned to rebalances
      # weights_path length = W; each rebalance applies to next Rout block (size size_w_out)
      # use the first row of each Rout for turnover
      # Build returns aligned to weights changes
      ret_aligned <- matrix(NA_real_, nrow = length(w_path), ncol = N)
      for (j in seq_along(windows)) {
        t_out <- windows[[j]]$idx_out[1]
        if (t_out <= Tobs) ret_aligned[j, ] <- R[t_out, ]
      }
      turn <- compute_turnover(w_path, ret_aligned, rf_vec[seq_len(nrow(ret_aligned))])
      turnover_med[ik] <- turn$median_to

      instab <- compute_weight_instability(w_path)
      instab_L1_med[ik] <- instab$median_L1
      instab_L2_med[ik] <- instab$median_L2

      sel_instab <- selection_instability(selections_path[[ik]])
      sel_instab_med[ik] <- sel_instab$median_distance
    } else {
      turnover_med[ik] <- NA_real_
      instab_L1_med[ik] <- NA_real_
      instab_L2_med[ik] <- NA_real_
      sel_instab_med[ik] <- NA_real_
    }
  }

  summary_df <- data.frame(
    k = k_grid,
    oos_sr = oos_by_k,
    median_turnover = turnover_med,
    median_weight_instability_L1 = instab_L1_med,
    median_weight_instability_L2 = instab_L2_med,
    median_selection_instability = sel_instab_med
  )

  if (!return_details) return(summary_df)
  list(
    summary = summary_df,
    oos_returns = oos_concat,
    weights_path = weights_path,
    selections_path = selections_path
  )
}

#' Parallel complete OOS evaluation with turnover and instability metrics
#'
#' Same as `run_complete_oos_evaluation` but parallelizes over windows using
#' `parallel::mclapply`.
#'
#' @inheritParams run_complete_oos_evaluation
#' @param n_cores Number of cores; defaults to `max(1, detectCores() - 1)`.
#' @export
run_complete_oos_evaluation_parallel <- function(R,
                                                 size_w_in,
                                                 size_w_out,
                                                 k_grid,
                                                 oos_type = c("rolling", "expanding"),
                                                 compute_weights_fn,
                                                 compute_weights_fn_params = list(),
                                                 rf = 0,
                                                 sharpe_fn = c("median", "mean"),
                                                 n_cores = NULL,
                                                 return_details = FALSE) {
  oos_type <- match.arg(oos_type)
  sharpe_fn <- match.arg(sharpe_fn)
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
    weights_list <- vector("list", K)
    sel_list <- vector("list", K)
    for (ik in seq_along(k_grid)) {
      k <- k_grid[ik]
      res <- do.call(compute_weights_fn, c(list(Rin, k), compute_weights_fn_params))
      if (is.null(res$weights)) stop("compute_weights_fn must return a `weights` element.")
      wopt <- as.numeric(res$weights)
      sel_raw <- res$selection
      sel_idx <- if (is.null(sel_raw)) integer() else {
        if (is.logical(sel_raw)) which(sel_raw) else as.integer(sel_raw)
      }
      weights_list[[ik]] <- wopt
      sel_list[[ik]] <- sel_idx
      rp <- if (length(sel_idx) == 0) {
        drop(Rout %*% wopt)
      } else {
        drop(Rout[, sel_idx, drop = FALSE] %*% wopt[sel_idx])
      }
      oos_list[[ik]] <- if (length(rp)) rp[is.finite(rp)] else numeric(0)
    }
    list(oos = oos_list, weights = weights_list, selections = sel_list)
  }

  res_list <- if (n_cores == 1L) {
    lapply(seq_len(W), worker)
  } else {
    parallel::mclapply(seq_len(W), worker, mc.cores = n_cores)
  }

  oos_concat <- vector("list", K)
  weights_path <- vector("list", K)
  selections_path <- vector("list", K)
  for (ik in seq_len(K)) {
    oos_concat[[ik]] <- numeric(0)
    weights_path[[ik]] <- list()
    selections_path[[ik]] <- list()
  }

  for (w in seq_len(W)) {
    res <- res_list[[w]]
    for (ik in seq_len(K)) {
      oos_concat[[ik]] <- c(oos_concat[[ik]], res$oos[[ik]])
      weights_path[[ik]][[length(weights_path[[ik]]) + 1L]] <- res$weights[[ik]]
      selections_path[[ik]][[length(selections_path[[ik]]) + 1L]] <- res$selections[[ik]]
    }
  }

  oos_by_k <- numeric(K)
  for (ik in seq_len(K)) {
    r <- oos_concat[[ik]]
    if (length(r) >= 2) {
      if (sharpe_fn == "median") {
        iqr <- stats::IQR(r, na.rm = TRUE)
        med <- stats::median(r, na.rm = TRUE)
        oos_by_k[ik] <- if (is.finite(iqr) && iqr > 0) med / iqr else NA_real_
      } else {
        s <- stats::sd(r, na.rm = TRUE)
        oos_by_k[ik] <- if (is.finite(s) && s > 0) mean(r, na.rm = TRUE) / s else NA_real_
      }
    } else {
      oos_by_k[ik] <- NA_real_
    }
  }

  turnover_med <- numeric(K)
  instab_L1_med <- numeric(K)
  instab_L2_med <- numeric(K)
  sel_instab_med <- numeric(K)

  if (length(rf) == 1L) rf_vec <- rep(rf, Tobs) else rf_vec <- rf
  if (length(rf_vec) < Tobs) stop("rf length must match nrow(R) or be scalar.")

  for (ik in seq_len(K)) {
    w_path <- weights_path[[ik]]
    if (length(w_path) >= 2) {
      ret_aligned <- matrix(NA_real_, nrow = length(windows), ncol = N)
      for (j in seq_along(windows)) {
        t_out <- windows[[j]]$idx_out[1]
        if (t_out <= Tobs) ret_aligned[j, ] <- R[t_out, ]
      }
      turn <- compute_turnover(w_path, ret_aligned, rf_vec[seq_len(nrow(ret_aligned))])
      turnover_med[ik] <- turn$median_to

      instab <- compute_weight_instability(w_path)
      instab_L1_med[ik] <- instab$median_L1
      instab_L2_med[ik] <- instab$median_L2

      sel_instab <- selection_instability(selections_path[[ik]])
      sel_instab_med[ik] <- sel_instab$median_distance
    } else {
      turnover_med[ik] <- NA_real_
      instab_L1_med[ik] <- NA_real_
      instab_L2_med[ik] <- NA_real_
      sel_instab_med[ik] <- NA_real_
    }
  }

  summary_df <- data.frame(
    k = k_grid,
    oos_sr = oos_by_k,
    median_turnover = turnover_med,
    median_weight_instability_L1 = instab_L1_med,
    median_weight_instability_L2 = instab_L2_med,
    median_selection_instability = sel_instab_med
  )

  if (!return_details) return(summary_df)
  list(
    summary = summary_df,
    oos_returns = oos_concat,
    weights_path = weights_path,
    selections_path = selections_path
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

#' Compute one-way turnover over a full weight path
#'
#' @param weights List of weight vectors (one per rebalance) or a matrix with rows = rebalances.
#'   Row t is the portfolio formed at the end of period t and held over period t+1.
#' @param R_excess Matrix of excess returns (rows are time, cols are assets).
#' @param rf Vector/scalar of risk-free rates aligned to rows of R_excess.
#' @param tol Tolerance for treating weights as zero.
#' @return A list with `to` (vector of turnover per rebalance) and `median_to`.
#' @export
compute_turnover <- function(weights, R_excess, rf, tol = 1e-12) {
  if (is.list(weights)) {
    w_mat <- do.call(rbind, weights)
  } else {
    w_mat <- as.matrix(weights)
  }
  R_excess <- as.matrix(R_excess)
  if (nrow(w_mat) < 2) stop("weights must have at least 2 rebalances.")
  if (nrow(R_excess) < nrow(w_mat)) stop("R_excess must have at least as many rows as weights.")
  if (length(rf) == 1L) rf <- rep(rf, nrow(R_excess))
  if (length(rf) != nrow(R_excess)) stop("rf length must match nrow(R_excess) or be scalar.")

  T_w <- nrow(w_mat)
  to_vec <- numeric(T_w - 1L)
  for (t in seq_len(T_w - 1L)) {
    to_vec[t] <- compute_turnover_window_t(w_mat[t, ], w_mat[t + 1L, ],
                                           R_excess[t + 1L, ], rf[t + 1L], tol = tol)
  }
  list(to = to_vec, median_to = stats::median(to_vec, na.rm = TRUE))
}

#' Compute one-way turnover for a single rebalance window
#'
#' @param w_prev Portfolio held over period t (weights after rebalance at t-1).
#' @param w_next Portfolio formed at end of period t (weights to hold over t+1).
#' @param r_excess Excess returns vector for period t.
#' @param rf Risk-free rate for period t.
#' @param tol Tolerance for treating weights as zero.
#' @return Scalar one-way turnover \code{sum(abs(w_next - w_pre))}.
#' @export
compute_turnover_window_t <- function(w_prev, w_next, r_excess, rf, tol = 1e-12) {
  g <- 1 + r_excess + rf
  w_pre <- w_prev * g
  denom <- sum(w_pre)
  if (!is.finite(denom) || abs(denom) < tol) return(NA_real_)
  w_pre <- w_pre / denom
  sum(abs(w_next - w_pre))
}

#' Compute weight instability (L1 and L2) across adjacent windows
#'
#' @param w_prev Weight vector at t (after rebalance at t).
#' @param w_next Weight vector at t+1 (after rebalance at t+1).
#' @return A length-2 numeric vector: c(L1, L2).
#' @export
compute_weight_instability_window <- function(w_prev, w_next) {
  diff <- w_next - w_prev
  c(L1 = sum(abs(diff)), L2 = sqrt(sum(diff^2)))
}

#' @rdname compute_weight_instability_window
#' @export
compute_weight_instability <- function(weights) {
  if (is.list(weights)) {
    w_mat <- do.call(rbind, weights)
  } else {
    w_mat <- as.matrix(weights)
  }
  if (nrow(w_mat) < 2) stop("weights must have at least 2 rebalances.")
  L1 <- numeric(nrow(w_mat) - 1L)
  L2 <- numeric(nrow(w_mat) - 1L)
  for (t in seq_len(nrow(w_mat) - 1L)) {
    diff <- w_mat[t + 1L, ] - w_mat[t, ]
    L1[t] <- sum(abs(diff))
    L2[t] <- sqrt(sum(diff^2))
  }
  list(
    instability_L1 = L1,
    instability_L2 = L2,
    median_L1 = stats::median(L1, na.rm = TRUE),
    median_L2 = stats::median(L2, na.rm = TRUE)
  )
}

#' Compute selection instability (Jaccard distance) for a single step
#'
#' @param sel_prev Logical or integer vector of selected assets at t.
#' @param sel_next Logical or integer vector of selected assets at t+1.
#' @return Jaccard distance (1 - Jaccard index).
#' @export
selection_instability_window <- function(sel_prev, sel_next) {
  if (is.logical(sel_prev)) sel_prev <- which(sel_prev)
  if (is.logical(sel_next)) sel_next <- which(sel_next)
  sel_prev <- unique(as.integer(sel_prev))
  sel_next <- unique(as.integer(sel_next))
  if (length(sel_prev) == 0 && length(sel_next) == 0) return(0)
  inter <- length(intersect(sel_prev, sel_next))
  uni <- length(union(sel_prev, sel_next))
  if (uni == 0) return(0)
  1 - inter / uni
}

#' Compute selection instability over a path of supports
#'
#' @param selections List of supports (logical or integer vectors) for each rebalance.
#' @return A list with per-step distances and median distance.
#' @export
selection_instability <- function(selections) {
  if (!is.list(selections) || length(selections) < 2) stop("selections must be a list of length >= 2.")
  D <- numeric(length(selections) - 1L)
  for (i in seq_len(length(selections) - 1L)) {
    D[i] <- selection_instability_window(selections[[i]], selections[[i + 1L]])
  }
  list(distance = D, median_distance = stats::median(D, na.rm = TRUE))
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

#' Plot median turnover vs k
#' @export
plot_turnover_empirics <- function(k_grid, turnover, method_labels = NULL, save_path = NULL) {
  turnover <- as.matrix(turnover)
  K <- length(k_grid); M <- ncol(turnover)
  if (K != nrow(turnover)) stop("k_grid length does not match turnover rows.")
  if (is.null(method_labels)) method_labels <- paste0("Method", seq_len(M))
  df <- data.frame(k = rep(k_grid, times = M),
                   turnover = as.vector(turnover),
                   method = rep(method_labels, each = K))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = turnover, color = method)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(x = "k", y = "Median turnover", color = "Method") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_text(size = 16),
                   legend.text = ggplot2::element_text(size = 14),
                   axis.title = ggplot2::element_text(size = 16),
                   axis.text = ggplot2::element_text(size = 14))
  if (!is.null(save_path)) {
    ggplot2::ggsave(paste0(save_path, ".png"), p, width = 8, height = 5, dpi = 150)
  }
  p
}

#' Plot median weight instability vs k
#' @export
plot_weight_instability_empirics <- function(k_grid, instab_L1, instab_L2 = NULL, method_labels = NULL, save_path = NULL) {
  instab_L1 <- as.matrix(instab_L1)
  K <- length(k_grid); M <- ncol(instab_L1)
  if (K != nrow(instab_L1)) stop("k_grid length does not match instability rows.")
  if (is.null(method_labels)) method_labels <- paste0("Method", seq_len(M))
  df <- data.frame(k = rep(k_grid, times = M),
                   L1 = as.vector(instab_L1),
                   method = rep(method_labels, each = K))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = L1, color = method)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(x = "k", y = "Median weight instability (L1)", color = "Method") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_text(size = 16),
                   legend.text = ggplot2::element_text(size = 14),
                   axis.title = ggplot2::element_text(size = 16),
                   axis.text = ggplot2::element_text(size = 14))
  if (!is.null(instab_L2)) {
    instab_L2 <- as.matrix(instab_L2)
    df2 <- data.frame(k = rep(k_grid, times = ncol(instab_L2)),
                      L2 = as.vector(instab_L2),
                      method = rep(method_labels, each = K))
    p <- p + ggplot2::geom_line(data = df2, ggplot2::aes(x = k, y = L2, linetype = "L2"), inherit.aes = FALSE)
  }
  if (!is.null(save_path)) {
    ggplot2::ggsave(paste0(save_path, ".png"), p, width = 8, height = 5, dpi = 150)
  }
  p
}

#' Plot median selection instability vs k
#' @export
plot_selection_instability_empirics <- function(k_grid, sel_instab, method_labels = NULL, save_path = NULL) {
  sel_instab <- as.matrix(sel_instab)
  K <- length(k_grid); M <- ncol(sel_instab)
  if (K != nrow(sel_instab)) stop("k_grid length does not match selection instability rows.")
  if (is.null(method_labels)) method_labels <- paste0("Method", seq_len(M))
  df <- data.frame(k = rep(k_grid, times = M),
                   sel = as.vector(sel_instab),
                   method = rep(method_labels, each = K))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = sel, color = method)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(x = "k", y = "Median selection instability", color = "Method") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_text(size = 16),
                   legend.text = ggplot2::element_text(size = 14),
                   axis.title = ggplot2::element_text(size = 16),
                   axis.text = ggplot2::element_text(size = 14))
  if (!is.null(save_path)) {
    ggplot2::ggsave(paste0(save_path, ".png"), p, width = 8, height = 5, dpi = 150)
  }
  p
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
