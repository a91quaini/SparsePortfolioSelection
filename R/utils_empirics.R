# utils_empirics.R
# Empirics utilities (theory-consistent: h = 1, roll by 1)
#
# This file is intentionally minimal: it keeps only the functions needed by
# scripts of the form emp_analysis_* that call:
#   load_data(), sps_make_output_dirs(), sps_subset_panel(), sps_make_k_grid(),
#   sps_make_stem_base(), run_empirical_suite_h1(), print_results(),
#   plot_empirical_suite_h1().

# =============================================================================
# Data loading
# =============================================================================

#' Load managed-portfolio panels (daily or monthly) + MEBEME panels (monthly)
#'
#' Returns a list with:
#'   - returns: matrix (T0 x N) of excess returns (and optionally appended factors)
#'   - rf: numeric vector length T0 (or NULL if unavailable)
#'
#' @param type One of: "US", "International", "International6", "International32", "MEBEME", "MEBEME_IND".
#' @param path Base `data` directory (default "data").
#' @param missing How to handle missing values: "keep" (default), "mean", "median", "remove".
#' @param frequency Either "daily" (default) or "monthly".
#' @param add_mkt Logical; if TRUE, append MKT column (panel-dependent).
#' @param add_factors Logical; if TRUE, append factor columns (panel-dependent).
#'   - For MEBEME/MEBEME_IND: FF5 = (MKT, SMB, HML, RMW, CMA).
#' @param shuffle Logical; if TRUE, randomly shuffle asset columns before returning.
#' @return A list: `returns` (matrix) and `rf` (numeric vector or NULL).
#' @export
load_data <- function(type = "US", path = "data", missing = "keep",
                      frequency = c("daily", "monthly"),
                      add_mkt = FALSE, add_factors = FALSE,
                      shuffle = TRUE) {

  freq <- match.arg(frequency)
  type <- as.character(type)
  type_key <- toupper(gsub("[^A-Za-z0-9]+", "_", type))
  miss_opt <- tolower(as.character(missing))

  if (!miss_opt %in% c("keep", "mean", "median", "remove")) {
    stop("missing must be one of: 'keep', 'mean', 'median', 'remove'.")
  }

  # -------------------------------------------------------------------
  # Special panels: MEBEME / MEBEME_IND (monthly only)
  # Directory: data/mebeme_ind_monthly
  # Files:
  #   - returns_mebeme100_monthly.rds
  #   - returns_ind49_monthly.rds
  #   - rf_mebeme_ind_monthly.rds
  #   - factors_ff5_monthly.rds
  # -------------------------------------------------------------------
  if (type_key %in% c("MEBEME", "MEBEME_IND", "MEBEMEIND", "MEBEME__IND")) {

    if (freq != "monthly") stop("type='", type, "' is only available at frequency='monthly'.")

    subdir <- "mebeme_ind_monthly"

    pkg_root <- system.file(package = "SparsePortfolioSelection")
    candidate_dirs <- c(
      file.path(path, subdir),
      system.file(subdir, package = "SparsePortfolioSelection"),
      if (nzchar(pkg_root)) file.path(pkg_root, "data", subdir) else character()
    )
    dir_path <- candidate_dirs[file.exists(candidate_dirs)][1]

    if (is.na(dir_path) || !file.exists(dir_path)) {
      if (identical(Sys.getenv("TESTTHAT"), "true")) {
        return(list(returns = matrix(0, nrow = 10, ncol = 5), rf = rep(0, 10)))
      }
      stop("No data found in: ", paste(candidate_dirs, collapse = ", "),
           ". Expected folder '", subdir, "'.")
    }

    # --- helpers for reading ---
    read_returns_rds <- function(fpath) {
      obj <- readRDS(fpath)
      m <- as.matrix(obj)
      dates <- NULL
      if (ncol(m) >= 2L) {
        # assume first column is dates if it is not numeric (or if it is clearly YYYYMM-like)
        first <- m[, 1]
        if (!is.numeric(first) || all(nchar(as.character(first)) %in% c(6L, 8L))) {
          dates <- as.vector(first)
          m <- m[, -1, drop = FALSE]
        }
      }
      list(dates = dates, mat = m)
    }

    read_rf_rds <- function(fpath, dates_ref = NULL) {
      obj <- readRDS(fpath)
      m <- as.matrix(obj)

      # if it's already a vector
      if (is.null(dim(m)) || ncol(m) == 1L) {
        rf <- as.numeric(m)
        if (!is.null(dates_ref) && length(rf) != length(dates_ref)) {
          stop("rf length does not match returns length.")
        }
        return(rf)
      }

      # if it has a date column + rf column
      if (ncol(m) >= 2L) {
        d <- m[, 1]
        rf <- as.numeric(m[, 2])
        if (!is.null(dates_ref)) {
          idx <- match(dates_ref, d)
          if (anyNA(idx)) stop("Date mismatch between returns and rf file: ", basename(fpath))
          rf <- rf[idx]
        }
        return(rf)
      }

      stop("Unsupported rf file format: ", basename(fpath))
    }

    # --- load returns ---
    f_mebeme <- file.path(dir_path, "returns_mebeme100_monthly.rds")
    if (!file.exists(f_mebeme)) stop("Missing file: ", f_mebeme)
    r1 <- read_returns_rds(f_mebeme)
    dates_ref <- r1$dates
    mat <- r1$mat

    if (type_key %in% c("MEBEME_IND", "MEBEMEIND", "MEBEME__IND")) {
      f_ind <- file.path(dir_path, "returns_ind49_monthly.rds")
      if (!file.exists(f_ind)) stop("Missing file: ", f_ind)
      r2 <- read_returns_rds(f_ind)
      if (!is.null(dates_ref) && !is.null(r2$dates) && !identical(dates_ref, r2$dates)) {
        stop("Date mismatch between mebeme100 and ind49 returns.")
      }
      if (is.null(dates_ref)) dates_ref <- r2$dates
      mat <- cbind(mat, r2$mat)
    }

    # --- load rf ---
    rf_file <- file.path(dir_path, "rf_mebeme_ind_monthly.rds")
    rf_vec <- if (file.exists(rf_file)) read_rf_rds(rf_file, dates_ref = dates_ref) else NULL

    # --- auto-rescale if looks like percent units ---
    # (Ken French data are typically in percent; your code expects decimals for turnover.)
    scale_to_decimal <- function(x) {
      if (is.null(x)) return(NULL)
      x <- as.numeric(x)
      if (!length(x)) return(x)
      med <- stats::median(abs(x[is.finite(x)]), na.rm = TRUE)
      if (is.finite(med) && med > 0.5) x <- x / 100
      x
    }
    # scale rf first (more reliable for unit detection)
    rf_vec <- scale_to_decimal(rf_vec)

    # scale returns using same heuristic on its own if rf missing
    mat_num <- as.numeric(mat)
    medR <- stats::median(abs(mat_num[is.finite(mat_num)]), na.rm = TRUE)
    if (is.finite(medR) && medR > 0.5) mat <- mat / 100

    # --- convert portfolios to excess returns if rf is available ---
    if (!is.null(rf_vec)) {
      if (nrow(mat) != length(rf_vec)) stop("rf length does not match returns rows.")
      mat <- sweep(mat, 1L, rf_vec, FUN = "-")
    }

    # --- optionally append FF5 factors (already excess/zero-cost factors) ---
    if (add_mkt || add_factors) {
      ff_file <- file.path(dir_path, "factors_ff5_monthly.rds")
      if (!file.exists(ff_file)) stop("factors file not found: ", ff_file)
      ff <- readRDS(ff_file)
      ffm <- as.matrix(ff)

      # align by date if possible
      if (!is.null(dates_ref) && ncol(ffm) >= 2L) {
        dff <- ffm[, 1]
        idx <- match(dates_ref, dff)
        if (anyNA(idx)) stop("Date mismatch between returns and factors_ff5_monthly.")
        ffm <- ffm[idx, , drop = FALSE]
      }

      # drop date col if present
      if (ncol(ffm) >= 2L) ff_vals <- ffm[, -1, drop = FALSE] else ff_vals <- ffm

      # scale factors if needed
      medF <- stats::median(abs(as.numeric(ff_vals)[is.finite(as.numeric(ff_vals))]), na.rm = TRUE)
      if (is.finite(medF) && medF > 0.5) ff_vals <- ff_vals / 100

      # select columns by name if possible, else assume first 5 are FF5
      cn <- colnames(ff_vals)
      pick_col <- function(candidates, fallback_idx) {
        if (!is.null(cn)) {
          j <- match(candidates, cn)
          j <- j[is.finite(j)][1]
          if (!is.na(j)) return(j)
        }
        fallback_idx
      }

      j_mkt <- pick_col(c("MKT", "Mkt-RF", "MKT-RF", "MKT_RF", "Mkt_RF", "Mkt.RF"), 1L)
      j_smb <- pick_col(c("SMB"), 2L)
      j_hml <- pick_col(c("HML"), 3L)
      j_rmw <- pick_col(c("RMW"), 4L)
      j_cma <- pick_col(c("CMA"), 5L)

      if (add_factors) {
        if (ncol(ff_vals) < max(j_mkt, j_smb, j_hml, j_rmw, j_cma)) {
          stop("factors_ff5_monthly does not have enough columns for FF5.")
        }
        X <- ff_vals[, c(j_mkt, j_smb, j_hml, j_rmw, j_cma), drop = FALSE]
        colnames(X) <- c("MKT", "SMB", "HML", "RMW", "CMA")
        mat <- cbind(mat, X)
      } else if (add_mkt) {
        if (ncol(ff_vals) < j_mkt) stop("factors_ff5_monthly missing MKT column.")
        mat <- cbind(mat, MKT = as.numeric(ff_vals[, j_mkt]))
      }
    }

    # missing handling
    if (miss_opt == "remove") {
      keep <- stats::complete.cases(mat)
      mat <- mat[keep, , drop = FALSE]
      if (!is.null(rf_vec)) rf_vec <- rf_vec[keep]
    } else if (miss_opt %in% c("mean", "median")) {
      for (j in seq_len(ncol(mat))) {
        idx_na <- which(is.na(mat[, j]))
        if (length(idx_na)) {
          filler <- if (miss_opt == "mean") mean(mat[, j], na.rm = TRUE) else stats::median(mat[, j], na.rm = TRUE)
          mat[idx_na, j] <- filler
        }
      }
    }

    if (shuffle) {
      perm <- sample.int(ncol(mat))
      mat <- mat[, perm, drop = FALSE]
    }

    return(list(returns = as.matrix(mat), rf = rf_vec))
  }

  # -------------------------------------------------------------------
  # Existing logic for US / International managed-portfolio panels
  # -------------------------------------------------------------------
  if (!(type %in% c("US", "International", "International6", "International32"))) {
    stop("type must be one of: 'US', 'International', 'International6', 'International32', 'MEBEME', 'MEBEME_IND'.")
  }

  if (type %in% c("International6", "International32") && freq != "monthly") {
    stop("type='", type, "' is only available at frequency='monthly'.")
  }

  subdir <- switch(
    freq,
    daily   = if (type == "US") "managed_portfolios_daily" else "managed_portfolios_international_daily",
    monthly = if (type == "US") "managed_portfolios_monthly" else "managed_portfolios_international_monthly"
  )

  pkg_root <- system.file(package = "SparsePortfolioSelection")
  candidate_dirs <- c(
    file.path(path, subdir),
    system.file(subdir, package = "SparsePortfolioSelection"),
    if (nzchar(pkg_root)) file.path(pkg_root, "data", subdir) else character()
  )
  dir_path <- candidate_dirs[file.exists(candidate_dirs)][1]

  if (is.na(dir_path) || !file.exists(dir_path)) {
    if (identical(Sys.getenv("TESTTHAT"), "true")) {
      return(list(returns = matrix(0, nrow = 10, ncol = 5), rf = NULL))
    }
    stop("No data found in: ", paste(candidate_dirs, collapse = ", "),
         ". Ensure processed .rds files exist.")
  }

  files <- list.files(dir_path, pattern = "^returns_.*\\.rds$", full.names = TRUE)

  if (type == "International" && freq == "monthly") {
    keep_basenames <- c(
      "returns_apxj_mebeme25_int_monthly.rds",
      "returns_eu_mebeme25_int_monthly.rds",
      "returns_jp_mebeme25_int_monthly.rds",
      "returns_na_mebeme25_int_monthly.rds",
      "returns_apxj_meop25_int_monthly.rds",
      "returns_eu_meop25_int_monthly.rds",
      "returns_jp_meop25_int_monthly.rds",
      "returns_na_meop25_int_monthly.rds"
    )
    files <- files[basename(files) %in% keep_basenames]
  } else if (type == "International6" && freq == "monthly") {
    keep_basenames <- c(
      "returns_apxj_mebeme6_int_monthly.rds",
      "returns_eu_mebeme6_int_monthly.rds",
      "returns_jp_mebeme6_int_monthly.rds",
      "returns_na_mebeme6_int_monthly.rds",
      "returns_apxj_meinv6_int_monthly.rds",
      "returns_eu_meinv6_int_monthly.rds",
      "returns_jp_meinv6_int_monthly.rds",
      "returns_na_meinv6_int_monthly.rds",
      "returns_apxj_meop6_int_monthly.rds",
      "returns_eu_meop6_int_monthly.rds",
      "returns_jp_meop6_int_monthly.rds",
      "returns_na_meop6_int_monthly.rds",
      "returns_apxj_meprior6_int_monthly.rds",
      "returns_eu_meprior6_int_monthly.rds",
      "returns_jp_meprior6_int_monthly.rds",
      "returns_na_meprior6_int_monthly.rds"
    )
    files <- files[basename(files) %in% keep_basenames]
  } else if (type == "International32" && freq == "monthly") {
    keep_basenames <- c(
      "returns_apxj_mebemeop32_int_monthly.rds",
      "returns_eu_mebemeop32_int_monthly.rds",
      "returns_jp_mebemeop32_int_monthly.rds",
      "returns_na_mebemeop32_int_monthly.rds"
    )
    files <- files[basename(files) %in% keep_basenames]
  }

  if (!length(files)) {
    if (identical(Sys.getenv("TESTTHAT"), "true")) {
      return(list(returns = matrix(0, nrow = 10, ncol = 5), rf = NULL))
    }
    stop("No returns_*.rds files found in ", dir_path, ".")
  }

  dates_ref <- NULL
  mats <- lapply(files, function(f) {
    df <- readRDS(f)
    if (ncol(df) >= 1) {
      if (is.null(dates_ref)) dates_ref <<- df[[1]]
      if (!is.null(dates_ref) && !identical(dates_ref, df[[1]])) {
        stop("Date mismatch across returns files; cannot align.")
      }
    }
    m <- as.matrix(df)
    if (ncol(m) >= 2) m[, -1, drop = FALSE] else m
  })

  n_rows <- unique(vapply(mats, nrow, integer(1)))
  if (length(n_rows) != 1) stop("Returns files have differing row counts: ", paste(n_rows, collapse = ", "))

  mat <- do.call(cbind, mats)
  rf_vec <- NULL

  # optionally append factors / rf (existing behavior)
  if (add_mkt || add_factors) {
    if (freq == "daily") {

      ff_file <- file.path(dir_path, "factors_ff5_daily.rds")
      if (!file.exists(ff_file)) stop("factors file not found: ", ff_file)
      ff <- readRDS(ff_file)

      if (is.null(dates_ref)) stop("Cannot align factors: missing date column in returns.")
      idx <- match(dates_ref, ff[[1]])
      if (anyNA(idx)) stop("Date mismatch between returns and factors_ff5_daily.")

      if (add_factors) {
        if (ncol(ff) < 4) stop("factors file missing FF3 columns.")
        X <- as.matrix(ff[idx, 2:4, drop = FALSE])
        colnames(X) <- c("MKT", "SMB", "HML")
        mat <- cbind(mat, X)
      } else if (add_mkt) {
        X <- as.numeric(ff[[2]][idx])
        mat <- cbind(mat, MKT = X)
      }

      if (ncol(ff) >= 5) rf_vec <- as.numeric(ff[[5]][idx])

    } else if (type %in% c("International", "International6", "International32") && freq == "monthly") {

      factor_files <- c(
        apxj   = "ff3_factors_monthly_asia_pacific_ex_japan.rds",
        europe = "ff3_factors_monthly_europe.rds",
        japan  = "ff3_factors_monthly_japan.rds",
        na     = "ff3_factors_monthly_north_america.rds"
      )

      cols_to_add <- list()
      have_rf <- FALSE

      for (region in names(factor_files)) {
        fpath <- file.path(dir_path, factor_files[[region]])
        if (!file.exists(fpath)) next

        ff <- readRDS(fpath)
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

      if (length(cols_to_add)) mat <- cbind(mat, as.matrix(as.data.frame(cols_to_add, check.names = FALSE)))
      if (!have_rf) rf_vec <- NULL

    } else {

      ff_file <- file.path(dir_path, sprintf("factors_ff5_%s.rds", freq))
      if (!file.exists(ff_file)) stop("factors file not found: ", ff_file)
      ff <- readRDS(ff_file)

      if (is.null(dates_ref)) stop("Cannot align factors: missing date column in returns.")
      idx <- match(dates_ref, ff[[1]])
      if (anyNA(idx)) stop("Date mismatch between returns and factors_ff5_", freq, ".")

      if (add_factors) {
        if (ncol(ff) < 4) stop("factors file missing FF3 columns.")
        X <- as.matrix(ff[idx, 2:4, drop = FALSE])
        colnames(X) <- c("MKT", "SMB", "HML")
        mat <- cbind(mat, X)
      } else if (add_mkt) {
        X <- as.numeric(ff[[2]][idx])
        mat <- cbind(mat, MKT = X)
      }

      if (ncol(ff) >= 5) rf_vec <- as.numeric(ff[[5]][idx])
    }
  }

  # missing handling
  if (miss_opt == "remove") {
    keep <- stats::complete.cases(mat)
    mat <- mat[keep, , drop = FALSE]
    if (!is.null(rf_vec)) rf_vec <- rf_vec[keep]
  } else if (miss_opt %in% c("mean", "median")) {
    for (j in seq_len(ncol(mat))) {
      idx_na <- which(is.na(mat[, j]))
      if (length(idx_na)) {
        filler <- if (miss_opt == "mean") mean(mat[, j], na.rm = TRUE) else stats::median(mat[, j], na.rm = TRUE)
        mat[idx_na, j] <- filler
      }
    }
  }

  if (shuffle) {
    perm <- sample.int(ncol(mat))
    mat <- mat[, perm, drop = FALSE]
  }

  list(returns = as.matrix(mat), rf = rf_vec)
}


# =============================================================================
# Internal helpers
# =============================================================================

.sps_as_rf_vec <- function(rf, T0) {
  if (is.null(rf)) return(NULL)
  if (length(rf) == 1L) return(rep(as.numeric(rf), T0))
  rf <- as.numeric(rf)
  if (length(rf) != T0) stop("rf must be scalar or length nrow(R_excess).")
  rf
}

.sps_sharpe_oos <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  mean(x) / s
}

.sps_jaccard_dist <- function(a, b) {
  a <- unique(as.integer(a)); b <- unique(as.integer(b))
  if (length(a) == 0L && length(b) == 0L) return(0)
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0L) return(0)
  1 - inter / uni
}

.sps_pre_trade_weights <- function(w_prev, r_excess_next, rf_next, tol = 1e-12) {
  G <- 1 + as.numeric(r_excess_next) + as.numeric(rf_next)
  x <- as.numeric(w_prev) * G
  denom <- sum(x)
  if (!is.finite(denom) || abs(denom) < tol) return(rep(NA_real_, length(w_prev)))
  x / denom
}

# Count, for each k in k_grid, whether the solver produced exactly k nonzero weights
# (within a tolerance). Returns an integer vector of length K with entries 0/1.
.sps_nnz_mismatch_vec <- function(W, k_grid, tol = 1e-10) {
  W <- as.matrix(W)
  k_grid <- as.integer(k_grid)
  if (nrow(W) != length(k_grid)) stop("nrow(W) must equal length(k_grid).")

  finite <- is.finite(W)
  nnz <- rowSums(finite & (abs(W) > tol))
  any_bad <- rowSums(!finite) > 0L

  as.integer(any_bad | (nnz != k_grid))
}

# expects compute_weights_fn to return either:
# - scalar-k: list(weights=vecN, selection=..., status=...)
# - vector-k: list(weights=matKxN, selection=listK, status=vecK)
.sps_parse_solver_output_vectorized <- function(res, N, K) {
  if (is.null(res$weights)) stop("compute_weights_fn must return $weights.")
  W <- res$weights

  if (is.matrix(W)) {
    if (nrow(W) != K || ncol(W) != N) stop("weights matrix must be K x N.")
    Sel <- res$selection
    if (!is.list(Sel) || length(Sel) != K) {
      Sel <- lapply(seq_len(K), function(i) which(abs(W[i, ]) > 0))
    }
    Status <- res$status
    if (is.null(Status)) Status <- rep("OK", K)
    if (length(Status) == 1L) Status <- rep(as.character(Status), K)
    if (length(Status) != K) Status <- rep("OK", K)
    return(list(W = W, Sel = Sel, Status = as.character(Status)))
  }

  # scalar output: treat as K=1
  w <- as.numeric(W)
  if (length(w) != N) stop("weights vector must have length N.")
  Sel <- res$selection
  if (is.null(Sel)) Sel <- which(abs(w) > 0)
  Status <- if (!is.null(res$status)) as.character(res$status) else "OK"
  list(W = matrix(w, nrow = 1L), Sel = list(as.integer(Sel)), Status = Status)

}

# =============================================================================
# Core empirical evaluation (theory-consistent h=1, roll by 1)
# =============================================================================

#' Run theory-consistent empirical design for h=1 (monthly rebalancing)
#'
#' Windows end at t = T_in,...,T0-1, weights formed at end of t and held over t+1.
#' Moments use denominator T_in:
#'   mu_t = (1/T_in) sum r_s
#'   Sigma_t = (1/T_in) sum r_s r_s' - mu_t mu_t'
#'
#' `compute_weights_fn` must be a function of the form:
#'   compute_weights_fn(mu, sigma, k) -> list(weights=..., selection=..., status=...)
#' (You should pre-bind any solver parameters via a closure.)
#'
#' @param R_excess T0 x N matrix of excess returns.
#' @param rf NULL, scalar, or length-T0 vector of risk-free rates (needed only for turnover).
#' @param T_in In-sample window length.
#' @param k_grid Integer vector of k values.
#' @param compute_weights_fn Function(mu, sigma, k) returning at least `$weights`.
#' @param tol_turnover Tolerance for the pre-trade normalization denominator.
#' @param annualize If TRUE, multiply monthly Sharpe by sqrt(12).
#' @param parallel If TRUE, parallelize over rebalancing dates.
#' @param n_cores Number of cores when parallel=TRUE.
#' @param parallel_backend "auto", "fork", or "psock".
#' @return data.frame with k-level summaries.
#' @export
run_empirical_design_h1 <- function(
    R_excess,
    rf = NULL,
    T_in,
    k_grid,
    compute_weights_fn,
    tol_turnover = 1e-12,
    annualize = TRUE,
    on_error = c("stop", "skip"),
    return_paths = FALSE,
    parallel = FALSE,
    n_cores = 1L,
    parallel_backend = c("auto", "fork", "psock")
) {
  parallel_backend <- match.arg(parallel_backend)
  on_error <- match.arg(on_error)

  R_excess <- as.matrix(R_excess)
  T0 <- nrow(R_excess)
  N <- ncol(R_excess)

  if (T0 < 3L || N < 1L) stop("R_excess must be T0 x N with T0>=3, N>=1.")
  T_in <- as.integer(T_in)
  if (T_in < 2L) stop("T_in must be >= 2.")
  if (T0 <= T_in) stop("Need T0 > T_in.")

  k_grid <- as.integer(k_grid)
  if (!length(k_grid)) stop("k_grid must be non-empty.")
  if (any(k_grid < 1L) || any(k_grid > N)) stop("Each k must be between 1 and N.")
  K <- length(k_grid)

  rf_vec <- .sps_as_rf_vec(rf, T0)
  have_rf <- !is.null(rf_vec)

  # rebalancing dates t = T_in,...,T0-1
  t_grid <- T_in:(T0 - 1L)
  n_reb <- length(t_grid)

  # store OOS returns as n_reb x K
  r_oos <- matrix(NA_real_, nrow = n_reb, ncol = K)


  # optionally store paths (weights/supports over time)
  if (isTRUE(return_paths)) {
    W_path <- vector("list", n_reb)
    S_path <- vector("list", n_reb)
    Status_path <- vector("list", n_reb)
  }
  # running previous weights/supports for instability + turnover
  prev_w <- matrix(NA_real_, nrow = K, ncol = N)
  prev_s <- vector("list", K)

  # metrics series (length n_reb-1)
  TO_mat <- matrix(NA_real_, nrow = max(0L, n_reb - 1L), ncol = K)
  dW1_mat <- matrix(NA_real_, nrow = max(0L, n_reb - 1L), ncol = K)
  dW2_mat <- matrix(NA_real_, nrow = max(0L, n_reb - 1L), ncol = K)
  dSel_mat <- matrix(NA_real_, nrow = max(0L, n_reb - 1L), ncol = K)


  # cardinality check: count windows where nnz(weights) != k
  nnz_mismatch_count <- integer(K)

  # worker: compute (mu,Sigma) for window ending at t, then solve for all k_grid
  worker_one <- function(j) {
    t <- t_grid[j]
    Rin <- R_excess[(t - T_in + 1L):t, , drop = FALSE]
    mu_t <- colMeans(Rin)
    Sigma_t <- (crossprod(Rin) / T_in) - tcrossprod(mu_t)

    out <- tryCatch(
      compute_weights_fn(mu = mu_t, sigma = Sigma_t, k = k_grid),
      error = function(e) e
    )

    if (inherits(out, "error")) {
      if (on_error == "stop") stop(out)
      Wbad <- matrix(NA_real_, nrow = K, ncol = N)
      Selbad <- replicate(K, integer(0), simplify = FALSE)
      return(list(W = Wbad, Sel = Selbad, Status = rep("ERROR", K)))
    }

    parsed <- .sps_parse_solver_output_vectorized(out, N = N, K = K)
    list(W = parsed$W, Sel = parsed$Sel, Status = parsed$Status)
  }

  use_parallel <- isTRUE(parallel) && n_cores > 1L && n_reb >= 2L
  if (parallel_backend == "auto") {
    parallel_backend <- if (.Platform$OS.type == "unix") "fork" else "psock"
  }

  if (!use_parallel) {
    res_list <- lapply(seq_len(n_reb), worker_one)
  } else if (parallel_backend == "fork") {
    res_list <- parallel::mclapply(seq_len(n_reb), worker_one, mc.cores = n_cores)
  } else {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      library(SparsePortfolioSelection)
      if (requireNamespace("lars", quietly = TRUE)) library(lars)
      NULL
    })
    parallel::clusterExport(
      cl,
      varlist = c("R_excess", "T_in", "t_grid", "k_grid", "K", "N",
                  "compute_weights_fn", "worker_one", ".sps_parse_solver_output_vectorized"),
      envir = environment()
    )
    res_list <- parallel::parLapply(cl, seq_len(n_reb), worker_one)
  }

  # combine + compute OOS returns + metrics in time order
  for (j in seq_len(n_reb)) {
    t <- t_grid[j]
    r_next <- R_excess[t + 1L, , drop = TRUE]

    res <- res_list[[j]]
    Wj <- res$W           # K x N
    Selj <- res$Sel       # list length K


    # check cardinality for this window (counts solver outputs with nnz != k)
    nnz_mismatch_count <- nnz_mismatch_count + .sps_nnz_mismatch_vec(Wj, k_grid)

    # OOS return: w_t' r_{t+1}

    if (isTRUE(return_paths)) {
      W_path[[j]] <- Wj
      S_path[[j]] <- Selj
      Status_path[[j]] <- res$Status
    }

    r_oos[j, ] <- as.numeric(Wj %*% r_next)

    if (j == 1L) {
      prev_w[,] <- Wj
      prev_s <- Selj
      next
    }

    for (ik in seq_len(K)) {
      w_prev <- prev_w[ik, ]
      w_curr <- Wj[ik, ]

      d <- w_curr - w_prev
      dW1_mat[j - 1L, ik] <- sum(abs(d))
      dW2_mat[j - 1L, ik] <- sqrt(sum(d * d))
      dSel_mat[j - 1L, ik] <- .sps_jaccard_dist(prev_s[[ik]], Selj[[ik]])

      if (have_rf) {
        # turnover at time t (rebalance end of t): uses returns at time t
        w_pre <- .sps_pre_trade_weights(w_prev, R_excess[t, ], rf_vec[t], tol = tol_turnover)
        TO_mat[j - 1L, ik] <- if (any(!is.finite(w_pre))) NA_real_ else sum(abs(w_curr - w_pre))
      }
    }

    prev_w[,] <- Wj
    prev_s <- Selj
  }

  sharpe_monthly <- apply(r_oos, 2L, .sps_sharpe_oos)
  sharpe_oos <- if (isTRUE(annualize)) sharpe_monthly * sqrt(12) else sharpe_monthly

  summary_df <- data.frame(
    T_in = T_in,
    k = k_grid,
    sharpe_oos = as.numeric(sharpe_oos),
    turnover_median = apply(TO_mat, 2L, stats::median, na.rm = TRUE),
    weight_instability_L1_median = apply(dW1_mat, 2L, stats::median, na.rm = TRUE),
    weight_instability_L2_median = apply(dW2_mat, 2L, stats::median, na.rm = TRUE),
    selection_instability_mean = colMeans(dSel_mat, na.rm = TRUE),
    nnz_mismatch_count = as.integer(nnz_mismatch_count)
  )

  if (isTRUE(return_paths)) {
    return(list(
      summary = summary_df,
      paths = list(
        W = W_path,
        S = S_path,
        Status = Status_path,
        t_grid = t_grid
      )
    ))
  }

  summary_df
}

# =============================================================================
# High-level orchestration helpers
# =============================================================================

#' Create output directories for an empirics run
#' @export
sps_make_output_dirs <- function(out_dir, fig_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(list(out_dir = out_dir, fig_dir = fig_dir))
}

#' Subset a panel to the first N assets (and align rf)
#' @export
sps_subset_panel <- function(R_excess, rf = NULL, n_assets = NULL) {
  R_excess <- as.matrix(R_excess)
  T0 <- nrow(R_excess)
  N0 <- ncol(R_excess)
  if (is.null(n_assets)) n_assets <- N0
  n_assets <- as.integer(n_assets)
  if (n_assets < 1L) stop("n_assets must be >= 1.")
  N <- min(n_assets, N0)

  rf_vec <- .sps_as_rf_vec(rf, T0)
  if (!is.null(rf_vec)) rf_vec <- rf_vec[seq_len(T0)]

  list(
    returns = R_excess[, 1:N, drop = FALSE],
    rf = rf_vec,
    T0 = T0,
    N = N
  )
}

#' Build a k-grid with basic checks
#' @export
sps_make_k_grid <- function(N, k_min = 1L, k_step = 3L, k_cap = NULL) {
  N <- as.integer(N)
  k_min <- as.integer(k_min)
  k_step <- as.integer(k_step)
  if (is.null(k_cap)) k_cap <- N - 1L
  k_cap <- as.integer(k_cap)

  if (N < 2L) stop("Need N >= 2 to form a meaningful k-grid.")
  if (k_min < 1L) stop("k_min must be >= 1.")
  if (k_step < 1L) stop("k_step must be >= 1.")
  k_max <- min(k_cap, N - 1L)
  if (k_max < k_min) stop("k_max < k_min; adjust k_min/k_cap or increase N.")
  seq.int(k_min, k_max, by = k_step)
}

#' Standard filename stem for an empirics run
#' @export
sps_make_stem_base <- function(method,
                               refit = FALSE,
                               panel_tag,
                               factors_tag,
                               mkt_tag,
                               N,
                               h = 1L,
                               normalize_weights = NULL) {
  refit_suffix <- if (isTRUE(refit)) "_refit" else ""
  normalize_suffix <- if (is.null(normalize_weights)) {
    ""
  } else if (isTRUE(normalize_weights)) {
    "_norm"
  } else {
    "_nonorm"
  }
  sprintf("oos_%s%s%s_%s_%s_%s_N%s_h%d",
          as.character(method), refit_suffix, normalize_suffix,
          as.character(panel_tag),
          as.character(factors_tag),
          as.character(mkt_tag),
          as.integer(N),
          as.integer(h))
}

#' Stack a metric from a list of summary data.frames into a K x M matrix
#' @export
sps_stack_metric <- function(summaries, metric_col, k_grid) {
  if (!length(summaries)) stop("summaries is empty.")
  k_grid <- as.integer(k_grid)
  mat <- vapply(summaries, function(df) {
    df <- as.data.frame(df)
    if (!("k" %in% names(df))) stop("Each summary must contain column 'k'.")
    if (!(metric_col %in% names(df))) stop("Missing column '", metric_col, "'.")
    idx <- match(k_grid, as.integer(df$k))
    if (anyNA(idx)) stop("k_grid does not match the summary's k column.")
    as.numeric(df[[metric_col]][idx])
  }, FUN.VALUE = numeric(length(k_grid)))
  mat <- as.matrix(mat)
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1L)
  mat
}

#' Run multiple T_in values, save CSVs, and return stacked matrices for plotting
#'
#' @param solver_factory function(T_in) -> compute_weights_fn(mu, sigma, k).
#' @param out_dir if not NULL, saves per-T_in CSVs there.
#' @param stem_base filename stem (required if out_dir is not NULL).
#' @return list(k_grid, labels, summaries, mats)
#' @export
run_empirical_suite_h1 <- function(R_excess,
                                   rf = NULL,
                                   T_in_grid,
                                   k_grid,
                                   solver_factory,
                                   annualize = TRUE,
                                   return_paths = FALSE,
                                   out_dir = NULL,
                                   stem_base = NULL,
                                   verbose = TRUE,
                                   ...) {

  R_excess <- as.matrix(R_excess)
  T0 <- nrow(R_excess)
  N <- ncol(R_excess)

  T_in_grid <- as.integer(T_in_grid)
  if (!length(T_in_grid)) stop("T_in_grid is empty.")
  k_grid <- as.integer(k_grid)
  if (!length(k_grid)) stop("k_grid is empty.")
  if (!is.function(solver_factory)) stop("solver_factory must be a function(T_in) -> compute_weights_fn.")

  if (!is.null(out_dir) && is.null(stem_base)) {
    stop("If out_dir is provided, stem_base must also be provided.")
  }

  summaries <- list()
  labels <- character(0)

  paths_by_Tin <- if (isTRUE(return_paths)) list() else NULL

  for (T_in in T_in_grid) {
    if (T_in >= T0) {
      if (isTRUE(verbose)) warning(sprintf("Skipping T_in=%d (>= T0=%d).", T_in, T0))
      next
    }

    if (isTRUE(verbose)) {
      message(sprintf("Running h=1 design: N=%d, T0=%d, T_in=%d, |k_grid|=%d", N, T0, T_in, length(k_grid)))
    }

    compute_weights_fn <- solver_factory(T_in)

    res <- run_empirical_design_h1(
      R_excess = R_excess,
      rf = rf,
      T_in = T_in,
      k_grid = k_grid,
      compute_weights_fn = compute_weights_fn,
      annualize = annualize,
      return_paths = return_paths,
      ...
    )

    if (isTRUE(return_paths)) {
      summary_df <- res$summary
      paths_by_Tin[[as.character(T_in)]] <- res$paths
    } else {
      summary_df <- res
    }

    summaries[[length(summaries) + 1L]] <- summary_df
    labels <- c(labels, as.character(T_in))

    if (!is.null(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      file <- file.path(out_dir, sprintf("%s_Tin%d.csv", stem_base, as.integer(T_in)))
      utils::write.csv(summary_df, file, row.names = FALSE)
    }
  }

  if (!length(summaries)) stop("No successful runs in run_empirical_suite_h1().")

  mats <- list(
    sharpe_oos = sps_stack_metric(summaries, "sharpe_oos", k_grid),
    turnover_median = sps_stack_metric(summaries, "turnover_median", k_grid),
    weight_instability_L1_median = sps_stack_metric(summaries, "weight_instability_L1_median", k_grid),
    weight_instability_L2_median = sps_stack_metric(summaries, "weight_instability_L2_median", k_grid),
    selection_instability_mean = sps_stack_metric(summaries, "selection_instability_mean", k_grid),
    nnz_mismatch_count = sps_stack_metric(summaries, "nnz_mismatch_count", k_grid)
  )

  checks <- list(
    nnz_mismatch_count = mats$nnz_mismatch_count
  )

  out <- list(
    k_grid = k_grid,
    labels = labels,
    summaries = summaries,
    mats = mats,
    checks = checks
  )
  if (isTRUE(return_paths)) out$paths <- paths_by_Tin
  out
}

# =============================================================================
# Plotting
# =============================================================================

.plot_metric_by_k <- function(k_grid, mat, labels, ylab,
                              log_y = FALSE, add_kopt = FALSE, save_path = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for plotting.")

  k_grid <- as.integer(k_grid)
  mat <- as.matrix(mat)

  if (nrow(mat) != length(k_grid)) stop("nrow(mat) must equal length(k_grid).")
  if (is.null(labels)) labels <- paste0("run", seq_len(ncol(mat)))
  if (length(labels) != ncol(mat)) stop("labels length must match ncol(mat).")

  K <- length(k_grid)
  M <- ncol(mat)

  df <- data.frame(
    k = rep(k_grid, times = M),
    value = as.numeric(mat),
    series = rep(labels, each = K)
  )

  # log-scale plots cannot display non-positive values: drop them quietly
  if (log_y) {
    df$value[!is.finite(df$value) | df$value <= 0] <- NA_real_
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = value, color = series)) +
    ggplot2::geom_line(linewidth = 1.1, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.2, na.rm = TRUE) +
    ggplot2::labs(
      x = "Number of holdings k",
      y = ylab,
      color = if (M > 1) "In-sample window" else NULL
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.position = if (M > 1) "bottom" else "none")

  if (log_y) p <- p + ggplot2::scale_y_log10()

  if (add_kopt) {
    kopt <- vapply(seq_len(M), function(j) {
      v <- mat[, j]
      v[!is.finite(v)] <- -Inf
      if (all(v == -Inf)) return(NA_integer_)
      k_grid[which.max(v)[1]]
    }, integer(1))

    vdf <- data.frame(series = labels, k_opt = kopt)
    vdf <- vdf[is.finite(vdf$k_opt), , drop = FALSE]

    if (nrow(vdf)) {
      vdf$y_opt <- vapply(seq_len(nrow(vdf)), function(i) {
        j <- match(vdf$series[i], labels)
        v <- mat[, j]
        v <- v[match(vdf$k_opt[i], k_grid)]
        if (length(v) && is.finite(v)) v else NA_real_
      }, numeric(1))

      if (log_y) vdf$y_opt[!is.finite(vdf$y_opt) | vdf$y_opt <= 0] <- NA_real_

      vdf <- vdf[is.finite(vdf$y_opt), , drop = FALSE]
      if (nrow(vdf)) {
        # Build once to learn the panel's y-axis minimum chosen by ggplot,
        # then draw a vertical dashed line from that minimum up to y_opt.
        pb <- ggplot2::ggplot_build(p)
        pp <- pb$layout$panel_params[[1]]
        y_min <- if (!is.null(pp$y.range)) pp$y.range[1] else pp$y$range[1]

        vdf$y_min <- y_min

        p <- p + ggplot2::geom_segment(
          data = vdf,
          ggplot2::aes(x = k_opt, xend = k_opt, y = y_min, yend = y_opt, color = series),
          linetype = "dashed", linewidth = 0.9, alpha = 0.85,
          inherit.aes = FALSE
        )
      }
    }
  }

  if (!is.null(save_path)) {
    ggplot2::ggsave(paste0(save_path, ".png"), p, width = 7.5, height = 5.0, dpi = 160)
  }

  p
}

plot_sr_empirics <- function(k_grid, SR, labels = NULL, save_path = NULL) {
  .plot_metric_by_k(k_grid, SR, labels, ylab = "OOS Sharpe ratio", log_y = FALSE, add_kopt = TRUE, save_path = save_path)
}

plot_turnover_empirics <- function(k_grid, turnover, labels = NULL, save_path = NULL) {
  .plot_metric_by_k(k_grid, turnover, labels, ylab = "Median one-way turnover", log_y = TRUE, add_kopt = FALSE, save_path = save_path)
}

plot_weight_instability_empirics <- function(k_grid, instab_L1, instab_L2, labels = NULL, save_path_base = NULL) {
  p1 <- .plot_metric_by_k(k_grid, instab_L1, labels, ylab = "Median weight instability (L1)", log_y = TRUE,
                          add_kopt = FALSE, save_path = if (is.null(save_path_base)) NULL else paste0(save_path_base, "_l1"))
  p2 <- .plot_metric_by_k(k_grid, instab_L2, labels, ylab = "Median weight instability (L2)", log_y = TRUE,
                          add_kopt = FALSE, save_path = if (is.null(save_path_base)) NULL else paste0(save_path_base, "_l2"))
  list(l1 = p1, l2 = p2)
}

plot_selection_instability_empirics <- function(k_grid, sel_instab, labels = NULL, save_path = NULL) {
  .plot_metric_by_k(k_grid, sel_instab, labels, ylab = "Mean selection instability", log_y = FALSE, add_kopt = FALSE, save_path = save_path)
}

#' Plot all standard empirics diagnostics from run_empirical_suite_h1()
#' @export
plot_empirical_suite_h1 <- function(suite, fig_dir, stem_base) {
  sps_make_output_dirs(out_dir = tempdir(), fig_dir = fig_dir) # ensures fig_dir exists

  k_grid <- suite$k_grid
  labels <- suite$labels

  plot_sr_empirics(
    k_grid, suite$mats$sharpe_oos, labels = labels,
    save_path = file.path(fig_dir, paste0(stem_base, "_sr"))
  )
  plot_turnover_empirics(
    k_grid, suite$mats$turnover_median, labels = labels,
    save_path = file.path(fig_dir, paste0(stem_base, "_turnover"))
  )
  plot_weight_instability_empirics(
    k_grid,
    suite$mats$weight_instability_L1_median,
    suite$mats$weight_instability_L2_median,
    labels = labels,
    save_path_base = file.path(fig_dir, paste0(stem_base, "_weight_instability"))
  )
  plot_selection_instability_empirics(
    k_grid, suite$mats$selection_instability_mean, labels = labels,
    save_path = file.path(fig_dir, paste0(stem_base, "_selection_instability"))
  )

  invisible(TRUE)
}

# =============================================================================
# Console printing
# =============================================================================

#' Quick console table for results matrices (e.g., Sharpe profiles)
#' @export
print_results <- function(k_grid, SR, method_labels = "value", digits = 4) {
  SR <- as.matrix(SR)
  if (length(k_grid) != nrow(SR)) stop("k_grid length must match nrow(SR).")

  if (ncol(SR) == 1L) {
    out <- data.frame(k = as.integer(k_grid), value = as.numeric(SR[, 1L]))
    names(out)[2] <- if (length(method_labels) >= 1L) method_labels[1L] else "value"
  } else {
    if (length(method_labels) != ncol(SR)) {
      stop("method_labels must have length = ncol(SR).")
    }
    out <- data.frame(k = as.integer(k_grid), SR)
    names(out)[-1L] <- method_labels
  }

  if (is.null(digits) || !is.finite(digits) || digits < 1) {
    print(out, row.names = FALSE)
  } else {
    print(out, row.names = FALSE, digits = as.integer(digits))
  }

  invisible(out)
}
