# Load the simulation results from file
f_name <- "portfolios_1fm_n20"
f_path <- file.path("inst", "simulations", "results", paste0("results_", f_name, ".rds"))
results <- readRDS(f_path)
N <- 20

# Extract the unique sample sizes (T) and cardinality levels (k)
T_vals <- sort(as.numeric(names(results)))
k_vals <- sort(as.numeric(unique(unlist(lapply(results, names)))))

# Prepare empty matrices to store summary statistics
mean_sr_loss_mat <- matrix(NA, nrow = length(T_vals), ncol = length(k_vals),
                           dimnames = list(paste0("T=", T_vals), paste0("k=", k_vals)))
mean_pct_est_mat <- matrix(NA, nrow = length(T_vals), ncol = length(k_vals),
                           dimnames = list(paste0("T=", T_vals), paste0("k=", k_vals)))
mean_pct_sel_mat <- matrix(NA, nrow = length(T_vals), ncol = length(k_vals),
                           dimnames = list(paste0("T=", T_vals), paste0("k=", k_vals)))

# Directory to save figures.
output_dir <- file.path("inst", "simulations", "figures")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each parameter combination and produce the plots and summary stats.
plot_list <- list()  # to store the arranged plots

for (T_val in T_vals) {
  for (k_val in k_vals) {
    sim_matrix <- results[[as.character(T_val)]][[as.character(k_val)]]
    if (is.null(sim_matrix)) next

    # Create a data frame from the simulation matrix.
    # Assumes sim_matrix has columns "sr_loss", "sr_loss_selection", "sr_loss_estimation".
    df <- as.data.frame(sim_matrix)

    # Compute percentage of estimation loss (if sr_loss is 0, set NA)
    df$pct_est <- with(df, ifelse(sr_loss != 0, sr_loss_estimation / sr_loss * 100, NA))

    # Create ggplot for overall Sharpe ratio loss with a gray fill and enlarged text.
    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = sr_loss)) +
      ggplot2::geom_histogram(fill = "grey70", color = "black", bins = 30) +
      ggplot2::labs(title = paste("Sharpe Ratio Loss: N =", N, ", T =", T_val, "and k =", k_val),
                    x = "Sharpe Ratio Loss", y = "Frequency") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(axis.text = ggplot2::element_text(size = 14),
                     axis.title = ggplot2::element_text(size = 16),
                     plot.title = ggplot2::element_text(size = 18, face = "bold"))

    # Create ggplot for percentage estimation loss with a gray fill and enlarged text.
    p2 <- ggplot2::ggplot(df, ggplot2::aes(x = pct_est)) +
      ggplot2::geom_histogram(fill = "grey70", color = "black", bins = 30) +
      ggplot2::labs(title = paste("% Estimation Loss: N =", N, ", T =", T_val, "and k =", k_val),
                    x = "Percentage Estimation Loss", y = "Frequency") +
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
    mean_sr_loss_mat[paste0("T=", T_val), paste0("k=", k_val)] <- mean(df$sr_loss, na.rm = TRUE)
    mean_pct_est_mat[paste0("T=", T_val), paste0("k=", k_val)] <- mean(df$pct_est, na.rm = TRUE)
    # For selection loss, compute percentage if desired.
    df$pct_sel <- with(df, ifelse(sr_loss != 0, sr_loss_selection / sr_loss * 100, NA))
    mean_pct_sel_mat[paste0("T=", T_val), paste0("k=", k_val)] <- mean(df$pct_sel, na.rm = TRUE)
  }
}

# Display the summary matrices.
cat("Mean Sharpe Ratio Loss:\n")
print(mean_sr_loss_mat)

cat("\nMean % Estimation Loss:\n")
print(mean_pct_est_mat)

cat("\nMean % Selection Loss:\n")
print(mean_pct_sel_mat)

