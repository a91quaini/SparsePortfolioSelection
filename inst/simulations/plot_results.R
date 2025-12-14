# Plot saved Monte Carlo results under inst/simulations/results.
# This script loads each results_*_mc_analysis_ff3.rds, plots the MVE SR
# decomposition, and saves PNGs alongside the source RDS files.

library(SparsePortfolioSelection)

results_dir <- file.path("inst", "simulations", "results")
files <- list.files(results_dir, pattern = "^results_ff3_.*_mc_analysis_ff3\\.rds$", full.names = TRUE)

if (length(files) == 0) {
  stop("No results_ff3_*.rds files found in ", results_dir)
}

for (f in files) {
  message("Processing ", basename(f))
  res <- readRDS(f)
  # Ensure required fields exist
  if (is.null(res$config) || is.null(res$population_sr) ||
      is.null(res$est_terms) || is.null(res$sel_terms)) {
    warning("Skipping malformed results file: ", f)
    next
  }
  plot_mve_sr_decomposition(res, save = TRUE)
}

message("Finished plotting all result files.")
