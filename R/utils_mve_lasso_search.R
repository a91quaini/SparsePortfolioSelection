# LASSO utilities for MVE search

# Build synthetic design X, y from moments (mu, Sigma, T) as in Julia utils.
.design_from_moments <- function(mu, sigma_s, Tobs) {
  design_from_moments_cpp(mu, sigma_s, Tobs, eps_ridge_cpp(), FALSE)
}

.validate_alpha_grid <- function(alpha) {
  if (is.null(alpha) || length(alpha) == 0) stop("alpha must be provided.")
  if (any(!is.finite(alpha))) stop("alpha must be finite.")
  grid <- as.numeric(alpha)
  list(grid = grid, is_grid = length(grid) > 1L)
}

.cv_alpha_return <- function(R, k, alpha_grid,
                             nlambda, lambda_min_ratio,
                             nadd, nnested,
                             standardize, epsilon, stabilize_sigma,
                             compute_weights, normalize_weights, use_refit,
                             n_folds) {
  # Minimal placeholder: pick the first alpha in the grid.
  # Can be extended to do actual CV if needed.
  alpha_grid[1]
}

.select_lasso_support <- function(X, y, k,
                                  nlambda, lambda_min_ratio,
                                  lambda_override,
                                  alpha,
                                  nadd = 80L,
                                  nnested = 2L,
                                  standardize = FALSE) {
  current_lambda_min_ratio <- lambda_min_ratio
  best <- NULL

  for (attempt in 1:3) {
    glm_args <- list(x = X, y = y, alpha = alpha, intercept = FALSE, standardize = standardize)
    if (!is.null(lambda_override)) {
      glm_args$lambda <- lambda_override
    } else {
      glm_args$nlambda <- as.integer(nlambda)
      glm_args$lambda.min.ratio <- current_lambda_min_ratio
    }
    fit <- do.call(glmnet::glmnet, glm_args)
    beta_mat <- as.matrix(fit$beta)
    lambdas <- fit$lambda
    nnz <- colSums(beta_mat != 0)

    feas <- which(nnz <= k)
    if (length(feas) > 0) {
      best_idx <- feas[which.max(nnz[feas])]
      status <- if (nnz[best_idx] == k) "LASSO_PATH_EXACT_K" else "LASSO_PATH_CLOSEST"
      best <- list(beta = beta_mat[, best_idx], support = which(beta_mat[, best_idx] != 0),
                   lambda = lambdas[best_idx], status = status,
                   beta_mat = beta_mat, lambdas = lambdas, nnz = nnz)
      break
    }

    if (max(nnz) < k && is.null(lambda_override)) {
      current_lambda_min_ratio <- current_lambda_min_ratio * 0.1
    } else {
      best_idx <- which.min(abs(nnz - k))
      best <- list(beta = beta_mat[, best_idx], support = which(beta_mat[, best_idx] != 0),
                   lambda = lambdas[best_idx], status = "LASSO_PATH_OVER_K",
                   beta_mat = beta_mat, lambdas = lambdas, nnz = nnz)
      break
    }
  }

  if (is.null(best)) {
    stop("Failed to fit glmnet path.")
  }

  # Densify between last <=k and first >k
  beta_mat <- best$beta_mat; lambdas <- best$lambdas; nnz <- best$nnz
  for (nest in seq_len(nnested)) {
    idx_lo <- max(which(nnz <= k))
    idx_hi <- min(which(nnz > k))
    if (is.infinite(idx_lo) || is.infinite(idx_hi) || idx_lo >= idx_hi) break
    lam_lo <- lambdas[idx_lo]
    lam_hi <- lambdas[idx_hi]
    lam_grid <- exp(seq(log(lam_hi), log(lam_lo), length.out = as.integer(nadd) + 2L))
    lam_grid <- lam_grid[-c(1, length(lam_grid))]  # drop endpoints
    fit_dense <- glmnet::glmnet(x = X, y = y, alpha = alpha, intercept = FALSE,
                                standardize = standardize, lambda = lam_grid)
    beta_dense <- as.matrix(fit_dense$beta)
    nnz_dense <- colSums(beta_dense != 0)
    feas_d <- which(nnz_dense <= k)
    if (length(feas_d) > 0) {
      bidx <- feas_d[which.max(nnz_dense[feas_d])]
      status <- if (nnz_dense[bidx] == k) "LASSO_PATH_EXACT_K" else "LASSO_PATH_CLOSEST"
      return(list(beta = beta_dense[, bidx],
                  support = which(beta_dense[, bidx] != 0),
                  lambda = lam_grid[bidx],
                  status = status))
    } else {
      dist <- abs(nnz_dense - k)
      bidx <- which.min(dist)
      if (nnz_dense[bidx] < nnz[idx_lo]) {
        best <- list(beta = beta_dense[, bidx],
                     support = which(beta_dense[, bidx] != 0),
                     lambda = lam_grid[bidx],
                     status = "LASSO_PATH_OVER_K")
        nnz <- nnz_dense
        lambdas <- lam_grid
        beta_mat <- beta_dense
      }
    }
  }

  list(beta = best$beta, support = best$support, lambda = best$lambda, status = best$status)
}
