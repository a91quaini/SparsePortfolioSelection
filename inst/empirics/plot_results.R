# Plot selected OOS results together and save figures using the same
# conventions as the analysis scripts (one combined plot per metric).

library(readr)
library(SparsePortfolioSelection)

# ---- User configuration ---------------------------------------------------
# Directory under inst/empirics/results where the CSVs live
RESULTS_SUBDIR <- "managed_portfolios_international_monthly" # "mebeme_ind_monthly"

# Filenames to combine (relative to RESULTS_SUBDIR)
FILES <- c(
  "oos_lars_nonorm_International_ff3_mkt_N212_h1_Tin240.csv",
  "oos_lars_nonorm_International_ff3_mkt_N212_h1_Tin360.csv"
)

# Labels for the legend (same order/length as FILES)
LABELS <- c("240", "360")

# Base stem for saved figures (do NOT include metric or extension);
# each metric appends _sr, _turnover, _weight_instability(_l1/_l2),
# _selection_instability.
FIG_STEM <- "oos_lars_nonorm_International_ff3_mkt_N212_h1"

# ---- No edits below -------------------------------------------------------

if (length(FILES) != length(LABELS)) {
  stop("FILES and LABELS must have the same length.")
}

results_dir <- file.path("inst", "empirics", "results", RESULTS_SUBDIR)
fig_dir <- file.path("inst", "empirics", "figures", RESULTS_SUBDIR)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

sr_list <- list()
turn_list <- list()
instab1_list <- list()
instab2_list <- list()
selinst_list <- list()
used_labels <- character(0)
k_ref <- NULL

for (i in seq_along(FILES)) {
  f <- FILES[i]
  lab <- LABELS[i]
  path <- file.path(results_dir, f)
  if (!file.exists(path)) {
    warning("Skipping missing file: ", path)
    next
  }
  dat <- readr::read_csv(path, show_col_types = FALSE)
  if (!"k" %in% names(dat)) {
    warning("Skipping file missing k column: ", path)
    next
  }
  if (is.null(k_ref)) {
    k_ref <- dat$k
  } else if (!identical(k_ref, dat$k)) {
    stop("All files must share the same k grid.")
  }

  # Sharpe ratio column: prefer sharpe_oos (new), otherwise SharpeRatio (legacy)
  if ("sharpe_oos" %in% names(dat)) {
    sr_col <- "sharpe_oos"
  } else if ("SharpeRatio" %in% names(dat)) {
    sr_col <- "SharpeRatio"
  } else {
    warning("Skipping file missing sharpe_oos/SharpeRatio column: ", path)
    next
  }

  sr_list[[length(sr_list) + 1L]] <- dat[[sr_col]]
  used_labels <- c(used_labels, lab)

  if ("turnover_median" %in% names(dat)) {
    turn_list[[length(turn_list) + 1L]] <- dat$turnover_median
  }
  if ("weight_instability_L1_median" %in% names(dat)) {
    instab1_list[[length(instab1_list) + 1L]] <- dat$weight_instability_L1_median
  }
  if ("weight_instability_L2_median" %in% names(dat)) {
    instab2_list[[length(instab2_list) + 1L]] <- dat$weight_instability_L2_median
  }
  if ("selection_instability_mean" %in% names(dat)) {
    selinst_list[[length(selinst_list) + 1L]] <- dat$selection_instability_mean
  }
}

if (length(sr_list) == 0) stop("No data to plot.")

base_path <- file.path(fig_dir, FIG_STEM)

sr_mat <- do.call(cbind, sr_list)
invisible(SparsePortfolioSelection:::plot_sr_empirics(k_ref, sr_mat, labels = used_labels,
                                                      save_path = paste0(base_path, "_sr")))

if (length(turn_list) == length(sr_list)) {
  turn_mat <- do.call(cbind, turn_list)
  invisible(SparsePortfolioSelection:::plot_turnover_empirics(k_ref, turn_mat, labels = used_labels,
                                                              save_path = paste0(base_path, "_turnover")))
}

if (length(instab1_list) == length(sr_list) && length(instab2_list) == length(sr_list)) {
  instab1_mat <- do.call(cbind, instab1_list)
  instab2_mat <- do.call(cbind, instab2_list)
  invisible(SparsePortfolioSelection:::plot_weight_instability_empirics(k_ref, instab1_mat, instab2_mat,
                                                                        labels = used_labels,
                                                                        save_path_base = paste0(base_path, "_weight_instability")))
}

if (length(selinst_list) == length(sr_list)) {
  sel_mat <- do.call(cbind, selinst_list)
  invisible(SparsePortfolioSelection:::plot_selection_instability_empirics(k_ref, sel_mat, labels = used_labels,
                                                                           save_path = paste0(base_path, "_selection_instability")))
}
