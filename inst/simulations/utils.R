################################################################################

# Calibrate a Factor Model
#
# This function calibrates a factor model given a TxN matrix of asset returns and a TxK matrix of factor returns.
# The function computes the population mean vector
# and covariance matrix of the factors, and uses these to compute regression coefficients via
# \eqn{\beta = (\Sigma_f)^{-1} \, \mathrm{Cov}(factors, returns)^T,}
# where \eqn{\Sigma_f} is the covariance matrix of the factors.
# The model-implied mean return vector is computed as
# \eqn{\mu = \beta \, \mu_f,}
# where \eqn{\mu_f} is the mean vector of the factors.
# Residuals are then computed (from demeaned returns) and their variances are
# used to form a diagonal matrix \eqn{\Sigma_0}.
# Finally, the covariance matrix of asset returns is estimated as
# \eqn{\Sigma = \beta \, \Sigma_f \, \beta^T + \Sigma_0.}
#
# The additional parameter \code{weak_coeff} adjusts the factor betas weakness by dividing them by
# \eqn{N^{(\text{weak\_coeff}/2)}}, where \eqn{N} is the number of assets.
# A value of \code{0} (default) indicates no adjustment,
# while a value of \code{1} indicates full weakness.
#
# The parameter \code{idiosy_vol_type} determines the structure of the idiosyncratic volatility:
# \itemize{
#   \item \code{0} (default): homoskedastic volatility, so \eqn{\Sigma_0} is the average residual variance times the identity.
#   \item \code{1}: heteroskedastic volatility (no correlation), so \eqn{\Sigma_0} is a diagonal matrix of asset-specific residual variances.
# }
#
# @param returns A numeric matrix (T x N) of asset returns.
# @param factors A numeric matrix (T x K) of factor returns.
# @param weak_coeff A numeric scalar between 0 and 1 indicating the weakness of the factors.
#                   A value of 0 (default) implies no weakness adjustment,
#                   while a value of 1 implies full adjustment.
# @param idiosy_vol_type A numeric scalar representing the type of idiosyncratic volatility:
#                   0 (default) for homoskedastic; 1 for heteroskedastic.
# @param do_checks Logical flag indicating whether to perform input validation (default is FALSE).
#
# @return A list with two components:
# \describe{
#   \item{\eqn{\mu}}{The model-implied mean return vector (an N x 1 matrix or a vector).}
#   \item{\eqn{\Sigma}}{The model-implied covariance matrix of asset returns (an N x N matrix).}
# }
#
# @examples
# \dontrun{
#   set.seed(123)
#   T <- 100  # number of time periods
#   N <- 5    # number of assets
#   K <- 3    # number of factors
#
#   # Create a TxK matrix of factor returns
#   factors <- matrix(rnorm(T * K), nrow = T, ncol = K)
#
#   # Create a TxN matrix of asset returns from a linear factor model plus noise
#   beta_true <- matrix(runif(N * K), nrow = N, ncol = K)
#   returns <- factors %*% t(beta_true) + matrix(rnorm(T * N, sd = 0.5), nrow = T, ncol = N)
#
#   # Calibrate the model with no weak factor adjustment (weak_coeff = 0)
#   # under homoskedastic volatility
#   model1 <- calibrate_factor_model(returns,
#                                    factors,
#                                    weak_coeff = 0,
#                                    idiosy_vol_type = 0,
#                                    do_checks = TRUE)
#
#   # Calibrate the model with moderate weak factor adjustment (weak_coeff = 0.5)
#   # under heteroskedastic volatility
#   model2 <- calibrate_factor_model(returns,
#                                    factors,
#                                    weak_coeff = 0.5,
#                                    idiosy_vol_type = 1,
#                                    do_checks = TRUE)
#
#   # Display the calibrated mean return vectors and covariance matrices
#   print(model1$mu)
#   print(model1$sigma)
#   print(model2$mu)
#   print(model2$sigma)
# }
#
calibrate_factor_model <- function(returns, factors, weak_coeff = 0.0, idiosy_vol_type = 0, do_checks = FALSE) {

  # Perform input validation if do_checks is TRUE
  if (do_checks) {
    if (missing(returns) || !is.matrix(returns) || nrow(returns) == 0 || !is.numeric(returns)) {
      stop("returns must be provided and be a non-empty numeric matrix")
    }
    if (missing(factors) || !is.matrix(factors) || nrow(factors) == 0 || !is.numeric(factors)) {
      stop("factors must be provided as a non-empty numeric matrix")
    }
    if (nrow(returns) != nrow(factors)) {
      stop("The number of rows in returns must equal the number of rows in factors")
    }
    if (!is.numeric(weak_coeff) || length(weak_coeff) != 1) {
      stop("weak_coeff must be a single numeric value")
    }
    if (weak_coeff < 0 || weak_coeff > 1) {
      stop("weak_coeff must be between 0 and 1")
    }
    if (!is.numeric(idiosy_vol_type) || length(idiosy_vol_type) != 1 || !(idiosy_vol_type %in% c(0, 1))) {
      stop("idiosy_vol_type must be a numeric scalar equal to 0 or 1")
    }
  }

  # Compute the mean vector and covariance matrix of the factors
  mu_f <- colMeans(factors)
  sigma_f <- stats::cov(factors)

  # Compute the regression coefficients beta (each column corresponds to an asset)
  beta <- t(solve(sigma_f, stats::cov(factors, returns))) / ncol(returns)^(weak_coeff / 2.0)

  # Compute the model-implied mean returns
  mu <- beta %*% mu_f

  # Compute the residuals (subtract the fitted values and the asset means)
  residuals <- returns - matrix(1, nrow(returns), 1) %*% colMeans(returns) - factors %*% t(beta)

  # Compute the residual variance for each asset and create a diagonal matrix
  if (idiosy_vol_type == 0) {
    sigma0 <- mean(apply(residuals, 2, stats::var)) * diag(ncol(returns))
  } else if (idiosy_vol_type == 1) {
    sigma0 <- diag(apply(residuals, 2, stats::var))
  }

  # Construct the covariance matrix of asset returns
  sigma <- beta %*% sigma_f %*% t(beta) + sigma0

  return(list(mu = mu, sigma = sigma))
}

################################################################################

# compute_simulation_results
#
# This function runs simulation experiments in parallel over a grid of parameters.
# For each combination of 'n_obs' and 'max_card', it runs 'n_sim' simulations by
# calling the 'simulate_sr_loss' function. The simulation results are organized into
# a nested list keyed first by observation count and then by maximum cardinality.
#
# Parameters:
#   n_obs           - A vector of observation counts to be used in the simulations.
#   n_sim           - An integer specifying the number of simulations to run per parameter combination.
#   mu              - A numeric vector or parameter required for the simulation.
#   sigma           - A numeric vector or parameter required for the simulation.
#   max_card        - A vector of maximum cardinality values.
#   max_comb        - An integer specifying the maximum number of combinations to consider. Default is 0,
#                     which means all combinations are computed.
#   simulate_mve_sr - A function that performs a simulation run. It should accept the parameters:
#                     mve_sr, mu, sigma, n_obs, and max_card, and return a list containing
#                     'sr_loss', 'sr_loss_selection', and 'sr_loss_estimation'.
#   seed            - An integer seed for random number generation for reproducibility (default is 123).
#   save_results    - A boolean indicating whether to save the simulation results to a file (default is TRUE).
#   file_name       - A string specifying the file name to save the results (default is "results_portfolios_1fm_n20.rds").
#
# Returns:
#   A nested list where the results are organized by observation count and maximum cardinality.
#   The structure is: results[[as.character(n_obs)]][[as.character(max_card)]],
#   with each element being a matrix of simulation results (n_sim x 4).
#   This matrix contains the following columns:
#     - \code{sample_mve_sr} computed as the optimal sample mve Sharpe ratio,
#     - \code{sample_mve_sr_cardk} computed as the optimal sample mve Sharpe ratio
#       with cardinality \code{max_card},
#     - \code{mve_sr_cardk_est_term} computed as \eqn{w^T \mu^T / \sqrt{w^T\sigma w}}
#       where \code{w} are the optimal sample mve weights,
#     - \code{mve_sr_cardk_sel_term} computed as \eqn{\mu_S^T  \sigma_S^{-1} \mu_S}
#       where \code{S} is the set of assets yielding the optimal sample mve Sharpe ratio.
#
# If 'save_results' is TRUE, the results will be saved to the file path:
# "inst/simulations/results/<file_name>".
#
compute_simulation_results <- function(n_obs,
                                       n_sim,
                                       mu,
                                       sigma,
                                       max_card,
                                       max_comb = 0,
                                       simulate_mve_sr,
                                       seed = 123,
                                       save_results = TRUE,
                                       file_name = "results_portfolios_1fm_n20.rds") {
  # Set the seed for reproducibility
  set.seed(seed)

  # Create the parameter grid using the supplied n_obs and max_card vectors
  param_grid <- expand.grid(n_obs = n_obs, max_card = max_card)

  # Choose number of cores to use: all but one (at least one)
  n_cores <- max(parallel::detectCores() - 1, 1)

  # Run the simulation in parallel over each parameter combination.
  # For each parameter combination the following is done:
  #   1. For the given n_obs and max_card, run n_sim simulations.
  #   2. Each simulation calls simulate_mve_sr using:
  #      - mu and sigma,
  #      - the current n_obs and max_card values.
  #      - max_comb (if 0, all combinations are computed).
  #   3. Each simulation returns a vector of two elements.
  results_grid <- parallel::mclapply(1:nrow(param_grid), function(i) {
    n_obs_val <- param_grid$n_obs[i]
    max_card_val <- param_grid$max_card[i]

    # Run n_sim simulations and collect results.
    sim_results <- replicate(n_sim, {
      output <- simulate_mve_sr(
        mu = mu,
        sigma = sigma,
        n_obs = n_obs_val,
        max_card = max_card_val,
        max_comb = max_comb,
        do_checks = FALSE
      )
      c(sample_mve_sr = output$sample_mve_sr,
        sample_mve_sr_cardk = output$sample_mve_sr_cardk,
        mve_sr_cardk_est_term = output$mve_sr_cardk_est_term,
        mve_sr_cardk_sel_term = output$mve_sr_cardk_sel_term)
    })

    # transpose the simulation result matrix from 2 x n_sim to n_sim x 2.
    sim_results <- t(sim_results)

    # Return a list with the parameter values and simulation matrix.
    list(n_obs = n_obs_val, max_card = max_card_val, sim_results = sim_results)
  }, mc.cores = n_cores)

  # Organize the results in a nested list keyed first by n_obs and then by max_card.
  results <- list()
  for (res in results_grid) {
    n_obs_key <- as.character(res$n_obs)
    max_card_key <- as.character(res$max_card)
    if (is.null(results[[n_obs_key]]))
      results[[n_obs_key]] <- list()
    results[[n_obs_key]][[max_card_key]] <- res$sim_results
  }

  # Save the results to the specified file.
  output_path <- file.path("inst", "simulations", "results", file_name)
  # Save the results only if save_results is TRUE.
  if (save_results) {
    saveRDS(results, file = output_path)
  }

  # Return results invisibly (or you can simply return the results)
  return(results)
}



################################################################################

# evaluate_simulation_results
#
# This function loads simulation results from a given file,
# computes summary statistics across different sample sizes (T) and
# cardinality constraints (k), and produces a set of histogram plots
# (one per combination of T and k) as well as summary plots.
#
# In addition, for each sample size T it creates a summary plot comparing:
#   a) Population full MVE Sharpe ratio (\(\theta_{*}\)), which is constant,
#   b) Population sparse MVE Sharpe ratio for each k (\(\theta_{k,*}\)),
#   c) Sample MVE Sharpe ratio (\(\theta(\hat{w}_k^{\text{card}})\)) computed as
#      \(\theta_{k,*}\) minus the mean overall loss,
#   d) Adjusted sample MVE Sharpe ratio (\(\theta_{H(\hat{w}_k^{\text{card}})}\)) computed as
#      \(\theta_{k,*}\) minus \(\text{mean sr loss} \times \big( 1 - \frac{\text{mean \% selection loss}}{100} \big)\).
#
# Moreover, the function computes quantile tables: for each combination of T and k,
# it computes the lower quantile at \(\alpha/2\) and the upper quantile at \(1-\alpha/2\)
# for both the overall Sharpe ratio loss and the percentage selection loss.
#
# The histogram and summary plots are saved as PNG files in "inst/simulations/figures".
#
# Inputs:
#   f_name        - Base file name (e.g., "portfolios_1fm_n20")
#   N             - Number of assets (used for labeling)
#   mve_sr        - Population full MVE Sharpe ratio (scalar)
#   mve_sr_cardk - Named list with MVE Sharpe ratios for each max cardinality
#                   (e.g., list("5" = ..., "10" = ..., "15" = ...))
#   alpha         - Significance level for quantile computation (default 0.05)
#
# Returns:
#   A list with the following matrices:
#     \item{mean_sr_loss_mat}{Matrix of mean Sharpe ratio loss for each sample size T and cardinality k.}
#     \item{mean_pct_est_mat}{Matrix of mean percentage estimation loss for each T and k.}
#     \item{mean_pct_sel_mat}{Matrix of mean percentage selection loss for each T and k.}
#     \item{quant_sr_loss_lower}{Matrix of the lower quantiles (\(\alpha/2\)) of Sharpe ratio loss.}
#     \item{quant_sr_loss_upper}{Matrix of the upper quantiles (\(1-\alpha/2\)) of Sharpe ratio loss.}
#     \item{quant_pct_est_lower}{Matrix of the lower quantiles (\(\alpha/2\)) of percentage estimation loss.}
#     \item{quant_pct_est_upper}{Matrix of the upper quantiles (\(1-\alpha/2\)) of percentage estimation loss.}
evaluate_simulation_results <- function(f_name = "portfolios_1fm_n20",
                                        N = 20,
                                        mve_sr,
                                        mve_sr_cardk,
                                        alpha = 0.05) {
  # Load the simulation results from file
  f_path <- file.path("inst", "simulations", "results", paste0("results_", f_name, ".rds"))
  results <- readRDS(f_path)

  # Extract the unique sample sizes (T) and cardinality levels (k)
  T_vals <- sort(as.numeric(names(results)))
  k_vals <- sort(as.numeric(unique(unlist(lapply(results, names)))))

  # Prepare empty matrices to store summary statistics
  mean_sr_loss_mat <- matrix(NA, nrow = length(T_vals), ncol = length(k_vals),
                             dimnames = list(paste0("T=", T_vals), paste0("k=", k_vals)))
  mean_pct_sel_mat <- matrix(NA, nrow = length(T_vals), ncol = length(k_vals),
                             dimnames = list(paste0("T=", T_vals), paste0("k=", k_vals)))

  # Prepare empty matrices for quantiles (overall sr_loss and pct_est)
  quant_sr_loss_lower <- matrix(NA, nrow = length(T_vals), ncol = length(k_vals),
                                dimnames = list(paste0("T=", T_vals), paste0("k=", k_vals)))
  quant_sr_loss_upper <- matrix(NA, nrow = length(T_vals), ncol = length(k_vals),
                                dimnames = list(paste0("T=", T_vals), paste0("k=", k_vals)))
  quant_pct_sel_lower <- matrix(NA, nrow = length(T_vals), ncol = length(k_vals),
                                dimnames = list(paste0("T=", T_vals), paste0("k=", k_vals)))
  quant_pct_sel_upper <- matrix(NA, nrow = length(T_vals), ncol = length(k_vals),
                                dimnames = list(paste0("T=", T_vals), paste0("k=", k_vals)))

  # Directory to save figures.
  output_dir <- file.path("inst", "simulations", "figures")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # List to store the arranged histogram plots for each (T, k)
  plot_list <- list()

  # Loop over each parameter combination to produce histograms and summary stats.
  for (T_val in T_vals) {
    for (k_val in k_vals) {
      sim_results <- results[[as.character(T_val)]][[as.character(k_val)]]
      if (is.null(sim_results)) next

      current_mve_sr_cardk <- mve_sr_cardk[[as.character(k_val)]]
      # Create a data frame from the simulation matrix.
      # Assumes sim_results has columns "mve_sr_cardk_est_term", "mve_sr_cardk_est_term"
      df <- as.data.frame(sim_results)

      # Compute the Sharpe ratio loss
      df$sr_loss <- with(df, current_mve_sr_cardk - sim_results[, "mve_sr_cardk_est_term"])

      # Compute percentage of selection loss (if sr_loss is 0, set 0)
      df$pct_sel <- with(df, ifelse(sr_loss != 0, (current_mve_sr_cardk - sim_results[, "mve_sr_cardk_sel_term"]) / df$sr_loss * 100, 0))

      # Compute quantiles for sr_loss and pct_est
      T_key <- paste0("T=", T_val)
      k_key <- paste0("k=", k_val)
      quant_sr_loss_lower[T_key, k_key] <- as.numeric(quantile(df$sr_loss, probs = alpha/2, na.rm = TRUE))
      quant_sr_loss_upper[T_key, k_key] <- as.numeric(quantile(df$sr_loss, probs = 1 - alpha/2, na.rm = TRUE))
      quant_pct_sel_lower[T_key, k_key] <- as.numeric(quantile(df$pct_sel, probs = alpha/2, na.rm = TRUE))
      quant_pct_sel_upper[T_key, k_key] <- as.numeric(quantile(df$pct_sel, probs = 1 - alpha/2, na.rm = TRUE))

      # Create ggplot for overall Sharpe ratio loss with a gray fill and enlarged text.
      p1 <- ggplot2::ggplot(df, ggplot2::aes(x = sr_loss)) +
        ggplot2::geom_histogram(fill = "grey70", color = "black", bins = 30) +
        ggplot2::labs(title = bquote(bold("Sharpe Ratio Loss" ~
                                       (theta[italic(k) * "," * "*" ] - theta(hat(w)[italic(k)]^{card})) ~
                                       ": N =" ~ .(N) ~ ", T =" ~ .(T_val) ~ " and k =" ~ .(k_val))),
                      x = "Sharpe Ratio Loss", y = "Frequency") +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 14),
                       axis.title = ggplot2::element_text(size = 16),
                       plot.title = ggplot2::element_text(size = 18, face = "bold"))

      # Create ggplot for percentage selection loss with a gray fill and enlarged text.
      p2 <- ggplot2::ggplot(df, ggplot2::aes(x = pct_sel)) +
        ggplot2::geom_histogram(fill = "grey70", color = "black", bins = 30) +
        ggplot2::labs(title = paste("% Selection Loss: N =", N, ", T =", T_val, "and k =", k_val),
                      x = "Percentage Selection Loss", y = "Frequency") +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 14),
                       axis.title = ggplot2::element_text(size = 16),
                       plot.title = ggplot2::element_text(size = 18, face = "bold"))

      # Arrange the two plots in a single panel using gridExtra.
      arranged <- gridExtra::grid.arrange(p1, p2, nrow = 2)

      # Save the arranged plot to a PNG file.
      filename <- file.path(output_dir, paste0(f_name, "_t", T_val, "_k", k_val, ".png"))
      ggplot2::ggsave(filename, arranged, width = 8, height = 10, dpi = 300)

      # Optionally store the arranged plot in a list.
      plot_list[[paste0("T", T_val, "_k", k_val)]] <- arranged

      # Compute mean overall sr_loss and mean percentages, removing NAs.
      mean_sr_loss_mat[T_key, k_key] <- mean(df$sr_loss, na.rm = TRUE)
      mean_pct_sel_mat[T_key, k_key] <- mean(df$pct_sel, na.rm = TRUE)

    }

    # Define initial factor levels of interest
    initial_types <- c("theta_*",
                       "theta_{k,*}",
                       "theta_{H(hat(w)_k^{card})}",
                       "theta(w_k^card)")

    # Create the initial data.frame with one row per combination of k and the 4 types
    summary_df <- data.frame(
      k = rep(k_vals, each = length(initial_types)),
      type = factor(rep(initial_types, times = length(k_vals)),
                    levels = initial_types),
      sharpe = NA_real_
    )

    # Specify T key (if needed later to index results)
    T_key <- paste0("T=", T_val)

    # Loop over each value of k to fill in the values for the initial types
    for (k_val in k_vals) {
      k_key <- paste0("k=", k_val)
      sim_results <- results[[as.character(T_val)]][[as.character(k_val)]]

      # a) Population full MVE Sharpe ratio (theta_*)
      summary_df$sharpe[summary_df$k == k_val & summary_df$type == "theta_*"] <- mve_sr

      # b) Population MVE Sharpe ratio with max cardinality k (theta_{k,*})
      summary_df$sharpe[summary_df$k == k_val & summary_df$type == "theta_{k,*}"] <-
        mve_sr_cardk[[as.character(k_val)]]

      # c) Population MVE Sharpe ratio with max cardinality k based on estimated index set (theta_{H(hat(w)_k^{card})})
      summary_df$sharpe[summary_df$k == k_val & summary_df$type == "theta_{H(hat(w)_k^{card})}"] <-
        mean(sim_results[, "mve_sr_cardk_sel_term"], na.rm = TRUE)

      # d) Population Sharpe ratio of sample MVE portfolio weights with max cardinality k (theta(w_k^card))
      summary_df$sharpe[summary_df$k == k_val & summary_df$type == "theta(w_k^card)"] <-
        mean(sim_results[, "mve_sr_cardk_est_term"], na.rm = TRUE)
      print(mean(sim_results[, "mve_sr_cardk_est_term"], na.rm = TRUE))
    }


    p_summary <- ggplot2::ggplot(summary_df, ggplot2::aes(x = k, y = sharpe, color = type, shape = type)) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_point(size = 3) +
      ggplot2::labs(title = paste("Sharpe Ratios for N =", N, ", T =", T_val),
                    x = "Cardinality (k)", y = "Sharpe Ratio") +
      ggplot2::scale_color_manual(
        name = "Type",
        values = c("theta_*" = "black",
                   "theta_{k,*}" = "blue",
                   "theta_{H(hat(w)_k^{card})}" = "darkgreen",
                   "theta(w_k^card)" = "red"),
        breaks = c("theta_*", "theta_{k,*}", "theta_{H(hat(w)_k^{card})}", "theta(w_k^card)"),
        labels = c(expression(theta["*"]),
                   bquote(theta[italic(k) * "," * "*" ]),
                   expression(theta[H(hat(w)[k]^"card")]),
                   expression(theta(hat(w)[k]^"card")))
      ) +
      ggplot2::scale_shape_manual(
        name = "Type",
        values = c("theta_*" = 16,
                   "theta_{k,*}" = 17,
                   "theta_{H(hat(w)_k^{card})}" = 18,
                   "theta(w_k^card)" = 15),
        breaks = c("theta_*", "theta_{k,*}", "theta_{H(hat(w)_k^{card})}", "theta(w_k^card)"),
        labels = c(expression(theta["*"]),
                   bquote(theta[italic(k) * "," * "*" ]),
                   expression(theta[H(hat(w)[k]^"card")]),
                   expression(theta(hat(w)[k]^"card")))
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(), shape = ggplot2::guide_legend()) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(axis.text = ggplot2::element_text(size = 14),
                     axis.title = ggplot2::element_text(size = 16),
                     plot.title = ggplot2::element_text(size = 18, face = "bold"),
                     legend.title = ggplot2::element_text(size = 18),
                     legend.text = ggplot2::element_text(size = 18))

    print(p_summary)

    # Save the summary plot.
    summary_filename <- file.path(output_dir, paste0(f_name, "_t", T_val, "_summary.png"))
    ggplot2::ggsave(summary_filename, p_summary, width = 8, height = 6, dpi = 300, bg = "white")

    # Define the new extra types to be added
    new_types <- c("hat(theta)", "hat(theta)_k")

    # Update factor levels to include the new types
    all_types <- c(initial_types, new_types)
    summary_df$type <- factor(summary_df$type, levels = all_types)

    # Create an extra data.frame for the new types for each value of k
    extra_df <- data.frame(
      k = rep(k_vals, each = length(new_types)),
      type = factor(rep(new_types, times = length(k_vals)), levels = all_types),
      sharpe = NA_real_
    )

    # Append the extra rows to the original summary_df
    summary_df <- rbind(summary_df, extra_df)

    # Loop again over each k value to fill in the values for the extra types
    for (k_val in k_vals) {
      k_key <- paste0("k=", k_val)
      sim_results <- results[[as.character(T_val)]][[as.character(k_val)]]

      # e) Sample MVE Sharpe ratio (hat(theta))
      summary_df$sharpe[summary_df$k == k_val & summary_df$type == "hat(theta)"] <-
        mean(sim_results[, "sample_mve_sr"], na.rm = TRUE)

      # f) Sample MVE Sharpe ratio with max cardinality k (hat(theta)_k)
      summary_df$sharpe[summary_df$k == k_val & summary_df$type == "hat(theta)_k"] <-
        mean(sim_results[, "sample_mve_sr_cardk"], na.rm = TRUE)
    }

    q_summary <- ggplot2::ggplot(summary_df, ggplot2::aes(x = k, y = sharpe, color = type, shape = type)) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_point(size = 3) +
      ggplot2::labs(title = paste("Sharpe Ratios for N =", N, ", T =", T_val),
                    x = "Cardinality (k)", y = "Sharpe Ratio") +
      ggplot2::scale_color_manual(
        name = "Type",
        values = c("theta_*" = "black",
                   "theta_{k,*}" = "blue",
                   "theta_{H(hat(w)_k^{card})}" = "darkgreen",
                   "theta(w_k^card)" = "red",
                   "hat(theta)" = "gray",
                   "hat(theta)_k" = "darkgray"),
        breaks = c("theta_*", "theta_{k,*}", "theta_{H(hat(w)_k^{card})}", "theta(w_k^card)", "hat(theta)", "hat(theta)_k"),
        labels = c(expression(theta["*"]),
                   bquote(theta[italic(k) * "," * "*" ]),
                   expression(theta[H(hat(w)[k]^"card")]),
                   expression(theta(hat(w)[k]^"card")),
                   expression(hat(theta)),
                   expression(hat(theta)[k]))
      ) +
      ggplot2::scale_shape_manual(
        name = "Type",
        values = c("theta_*" = 16,
                   "theta_{k,*}" = 17,
                   "theta_{H(hat(w)_k^{card})}" = 18,
                   "theta(w_k^card)" = 15,
                   "hat(theta)" = 13,
                   "hat(theta)_k" = 12),
        breaks = c("theta_*", "theta_{k,*}", "theta_{H(hat(w)_k^{card})}", "theta(w_k^card)", "hat(theta)", "hat(theta)_k"),
        labels = c(expression(theta["*"]),
                   bquote(theta[italic(k) * "," * "*" ]),
                   expression(theta[H(hat(w)[k]^"card")]),
                   expression(theta(hat(w)[k]^"card")),
                   expression(hat(theta)),
                   expression(hat(theta)[k]))
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(), shape = ggplot2::guide_legend()) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(axis.text = ggplot2::element_text(size = 14),
                     axis.title = ggplot2::element_text(size = 16),
                     plot.title = ggplot2::element_text(size = 18, face = "bold"),
                     legend.title = ggplot2::element_text(size = 18),
                     legend.text = ggplot2::element_text(size = 18))

    print(q_summary)

    # Save the summary plot.
    summary_filename <- file.path(output_dir, paste0(f_name, "_t", T_val, "_summary_full.png"))
    ggplot2::ggsave(summary_filename, q_summary, width = 8, height = 6, dpi = 300, bg = "white")
  }

  # Display the summary matrices.
  cat("Mean Sharpe Ratio Loss:\n")
  print(mean_sr_loss_mat)
  cat("\nQuantiles of sr_loss - Lower (alpha/2):\n")
  print(quant_sr_loss_lower)
  cat("\nQuantiles of sr_loss + Upper (1-alpha/2):\n")
  print(quant_sr_loss_upper)
  cat("\nMean % Selection Loss:\n")
  print(mean_pct_sel_mat)
  cat("\nQuantiles of % Selection Loss - Lower (alpha/2):\n")
  print(quant_pct_sel_lower)
  cat("\nQuantiles of % Selection Loss + Upper (1-alpha/2):\n")
  print(quant_pct_sel_upper)

  return(list(mean_sr_loss_mat = mean_sr_loss_mat,
              quant_sr_loss_lower = quant_sr_loss_lower,
              quant_sr_loss_upper = quant_sr_loss_upper,
              mean_pct_sel_mat = mean_pct_sel_mat,
              quant_pct_sel_lower = quant_pct_sel_lower,
              quant_pct_sel_upper = quant_pct_sel_upper))
}
