# Helpers for LASSO relaxation search (alpha validation and K-fold CV)

.validate_alpha_grid <- function(alpha) {
  if (is.null(alpha)) stop("alpha must be provided.")
  if (is.numeric(alpha) && length(alpha) == 1L) {
    if (alpha < 0 || alpha > 1) stop("alpha must be in [0,1].")
    list(is_grid = FALSE, grid = alpha)
  } else {
    agrid <- sort(unique(as.numeric(alpha)))
    if (any(!is.finite(agrid)) || any(agrid < 0) || any(agrid > 1)) {
      stop("alpha grid must be finite and within [0,1].")
    }
    list(is_grid = TRUE, grid = agrid)
  }
}

.kfold_split <- function(n, n_folds = 5L) {
  if (n_folds < 2L) stop("n_folds must be >= 2.")
  idx <- sample.int(n)
  folds <- split(idx, rep_len(seq_len(n_folds), n))
  lapply(seq_len(n_folds), function(f) {
    val <- folds[[f]]
    train <- setdiff(idx, val)
    list(train = train, val = val)
  })
}

.fit_glmnet_path <- function(X, y, alpha, nlambda, lambda_min_ratio, lambda_override, standardize) {
  args <- list(x = X, y = y, alpha = alpha, intercept = FALSE, standardize = standardize)
  if (is.null(lambda_override)) {
    args$nlambda <- as.integer(nlambda)
    args$lambda.min.ratio <- lambda_min_ratio
  } else {
    args$lambda <- lambda_override
  }
  do.call(glmnet::glmnet, args)
}

.select_support_closest_k <- function(beta_mat, lambdas, k) {
  nnz <- colSums(beta_mat != 0)
  feas <- which(nnz <= k)
  if (length(feas) > 0) {
    j <- feas[which.max(nnz[feas])]
    status <- if (nnz[j] == k) "LASSO_PATH_EXACT_K" else "LASSO_PATH_CLOSEST"
  } else {
    j <- which.min(abs(nnz - k))
    status <- "LASSO_PATH_OVER_K"
  }
  list(beta = beta_mat[, j], support = which(beta_mat[, j] != 0), lambda = lambdas[j], status = status)
}

.cv_alpha_glmnet <- function(X, y, k, alpha_grid, nlambda, lambda_min_ratio, standardize, n_folds = 5L) {
  folds <- .kfold_split(nrow(X), n_folds)
  best_alpha <- alpha_grid[1]
  best_mean <- -Inf

  for (a in alpha_grid) {
    sr_vals <- c()
    for (fold in folds) {
      fit <- .fit_glmnet_path(X[fold$train, , drop = FALSE], y[fold$train], a,
                              nlambda, lambda_min_ratio, NULL, standardize)
      beta_mat <- as.matrix(fit$beta)
      lambdas <- fit$lambda
      sel <- .select_support_closest_k(beta_mat, lambdas, k)
      w <- numeric(ncol(X))
      w[sel$support] <- sel$beta[sel$support]
      if (all(w == 0)) next
      sr <- as.numeric(t(w) %*% colMeans(X[fold$val, , drop = FALSE])) /
        sqrt(as.numeric(t(w) %*% stats::cov(X[fold$val, , drop = FALSE]) %*% w))
      if (is.finite(sr)) sr_vals <- c(sr_vals, sr)
    }
    if (length(sr_vals) > 0) {
      m <- mean(sr_vals)
      if (m > best_mean) {
        best_mean <- m
        best_alpha <- a
      }
    }
  }
  best_alpha
}

.design_from_moments_R <- function(mu, sigma_s, Tobs) {
  Q <- Tobs * (sigma_s + tcrossprod(mu))
  muQ <- mean(diag(Q))
  tau <- .Machine$double.eps * if (is.finite(muQ) && muQ > 0) muQ else 1
  diag(Q) <- diag(Q) + tau
  U <- tryCatch(chol(Q), error = function(e) NULL)
  if (is.null(U)) {
    U <- safe_chol_cpp(Q, base_bump = 1e-10, max_tries = 8)
  }
  X <- t(U)
  y <- backsolve(U, Tobs * mu)
  list(X = X, y = y)
}

.cv_alpha_return <- function(R, k, alpha_grid,
                             nlambda, lambda_min_ratio,
                             nadd, nnested,
                             standardize, epsilon, stabilize_sigma,
                             compute_weights, normalize_weights, use_refit,
                             n_folds = 5L) {
  folds <- .kfold_split(nrow(R), n_folds)
  scores <- rep(-Inf, length(alpha_grid))

  for (ai in seq_along(alpha_grid)) {
    a <- alpha_grid[ai]
    sr_vals <- c()
    for (fold in folds) {
      mu_tr <- colMeans(R[fold$train, , drop = FALSE])
      sigma_tr <- stats::cov(R[fold$train, , drop = FALSE])
      sigma_s <- prep_covariance_cpp(sigma_tr, epsilon, stabilize_sigma)
      design <- .design_from_moments_R(mu_tr, sigma_s, length(fold$train))
      sel <- .select_lasso_support(design$X, design$y, k,
                                   nlambda, lambda_min_ratio,
                                   lambda_override = NULL,
                                   alpha = a,
                                   nadd = nadd,
                                   nnested = nnested,
                                   standardize = standardize)

      weights <- numeric(ncol(R))
      if (length(sel$support) > 0) {
        if (compute_weights && use_refit) {
          weights <- compute_mve_weights_cpp(mu_tr, sigma_s,
                                             selection = as.integer(sel$support - 1),
                                             normalize_w = normalize_weights,
                                             epsilon = epsilon,
                                             stabilize_sigma = FALSE,
                                             do_checks = FALSE)
        } else {
          weights[sel$support] <- sel$beta[sel$support]
          if (normalize_weights) {
            weights <- normalize_weights_cpp(weights, mode = "relative", tol = 1e-6, do_checks = FALSE)
          }
        }
      }

      mu_val <- colMeans(R[fold$val, , drop = FALSE])
      sigma_val <- stats::cov(R[fold$val, , drop = FALSE])
      sr <- compute_sr_cpp(weights, mu_val, sigma_val,
                           selection = integer(),
                           epsilon = epsilon,
                           stabilize_sigma = FALSE,
                           do_checks = FALSE)
      if (is.finite(sr)) sr_vals <- c(sr_vals, sr)
    }
    if (length(sr_vals) > 0) {
      scores[ai] <- mean(sr_vals)
    }
  }

  alpha_grid[which.max(scores)]
}
