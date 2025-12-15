# Plot selected OOS results together and save figures using the same
# conventions as the analysis scripts (one combined plot per metric).

library(readr)
library(SparsePortfolioSelection)

# ---- User configuration ---------------------------------------------------
# Directory under inst/empirics/results where the CSVs live
RESULTS_SUBDIR <- "managed_portfolios_monthly" # "mebeme_ind_monthly"

# Filenames to combine (relative to RESULTS_SUBDIR)
FILES <- c(
  "combined_oos_lasso_norefit_us_ff3_mkt_Wout1_Win240_sr.csv",
  "combined_oos_lasso_norefit_us_ff3_mkt_Wout1_Win2360_sr.csv",
  "combined_oos_lasso_norefit_us_ff3_mkt_Wout1_Win480_sr.csv"
)

# Labels for the legend (same order/length as FILES)
LABELS <- c("240", "360", "480")

# Base stem for saved figures (do NOT include metric or extension);
# each metric appends _sr, _turnover, _weight_instability(_l1/_l2),
# _selection_instability, _less_than_k.
FIG_STEM <- "combined_oos_lasso_norefit_us_ff3_mkt_Wout1"

# Optional: total number of OOS windows for each file (needed to plot less_than_k)
# Set to NULL to skip the less-than-k plot, or supply a numeric vector
# of the same length as FILES.
TOTAL_WINDOWS <- NULL

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
lessk_list <- list()
lessk_totals <- numeric(0)
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

  # Sharpe ratio column: prefer SharpeRatio, otherwise first non-k column
  if ("SharpeRatio" %in% names(dat)) {
    sr_col <- "SharpeRatio"
  } else {
    sr_candidates <- setdiff(names(dat), "k")
    if (length(sr_candidates) == 1) {
      sr_col <- sr_candidates
    } else {
      warning("Skipping file missing recognizable SR column: ", path)
      next
    }
  }

  sr_list[[length(sr_list) + 1L]] <- dat[[sr_col]]
  used_labels <- c(used_labels, lab)

  if ("turnover" %in% names(dat)) {
    turn_list[[length(turn_list) + 1L]] <- dat$turnover
  }
  if ("weight_instability_l1" %in% names(dat)) {
    instab1_list[[length(instab1_list) + 1L]] <- dat$weight_instability_l1
  }
  if ("weight_instability_l2" %in% names(dat)) {
    instab2_list[[length(instab2_list) + 1L]] <- dat$weight_instability_l2
  }
  if ("selection_instability" %in% names(dat)) {
    selinst_list[[length(selinst_list) + 1L]] <- dat$selection_instability
  }
  if ("less_than_k" %in% names(dat)) {
    lessk_list[[length(lessk_list) + 1L]] <- dat$less_than_k
    if (!is.null(TOTAL_WINDOWS) && length(TOTAL_WINDOWS) >= i && !is.na(TOTAL_WINDOWS[i])) {
      lessk_totals <- c(lessk_totals, TOTAL_WINDOWS[i])
    } else {
      lessk_totals <- c(lessk_totals, NA_real_)
    }
  }
}

if (length(sr_list) == 0) stop("No data to plot.")

base_path <- file.path(fig_dir, FIG_STEM)

sr_mat <- do.call(cbind, sr_list)
invisible(plot_sr_empirics(k_ref, sr_mat, labels = used_labels,
                           save_path = paste0(base_path, "_sr")))

if (length(turn_list) == length(sr_list)) {
  turn_mat <- do.call(cbind, turn_list)
  invisible(plot_turnover_empirics(k_ref, turn_mat, labels = used_labels,
                                   save_path = paste0(base_path, "_turnover")))
}

if (length(instab1_list) == length(sr_list) && length(instab2_list) == length(sr_list)) {
  instab1_mat <- do.call(cbind, instab1_list)
  instab2_mat <- do.call(cbind, instab2_list)
  invisible(plot_weight_instability_empirics(k_ref, instab1_mat, instab2_mat,
                                             labels = used_labels,
                                             save_path_base = paste0(base_path, "_weight_instability")))
}

if (length(selinst_list) == length(sr_list)) {
  sel_mat <- do.call(cbind, selinst_list)
  invisible(plot_selection_instability_empirics(k_ref, sel_mat, labels = used_labels,
                                                save_path = paste0(base_path, "_selection_instability")))
}

# less-than-k needs total windows to compute fractions
if (length(lessk_list) == length(sr_list) &&
    length(lessk_totals) == length(sr_list) &&
    all(!is.na(lessk_totals))) {
  lessk_mat <- do.call(cbind, lessk_list)
  invisible(plot_less_than_k(k_ref, lessk_mat, labels = used_labels,
                             save_path = paste0(base_path, "_less_than_k"),
                             total_windows = lessk_totals))
} else if (length(lessk_list) > 0) {
  warning("Skipping less-than-k plot: TOTAL_WINDOWS not supplied or mismatched.")
}
