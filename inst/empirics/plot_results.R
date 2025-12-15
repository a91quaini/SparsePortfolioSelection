# Plot selected OOS Sharpe ratio results on the same figure.
# Configure the source folder, a set of CSV filenames, and labels below.

library(ggplot2)
library(readr)

# ---- User configuration ---------------------------------------------------
RESULTS_SUBDIR <- "mebeme_ind_monthly"  # folder under inst/empirics/results "mebeme_ind
FILES <- c(
  "oos_sr_lasso_mebeme_monthly_ff3_mkt_N103_Win240_Wout1.csv",
  "oos_sr_lasso_mebeme_monthly_ff3_mkt_N103_Win360_Wout1.csv",
  "oos_sr_lasso_mebeme_monthly_ff3_mkt_N103_Win480_Wout1.csv"
)
LABELS <- c("240","360","480")  # same length/order as FILES
FIG_NAME <- "combined_oos_sr_lasso_mebeme_monthly_ff3_mkt_N103_Wout1.png"

# ---- No edits below -------------------------------------------------------

if (length(FILES) != length(LABELS)) {
  stop("FILES and LABELS must have the same length.")
}

results_dir <- file.path("inst", "empirics", "results", RESULTS_SUBDIR)
fig_dir <- file.path("inst", "empirics", "figures", RESULTS_SUBDIR)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

plot_df <- data.frame()
max_points <- data.frame()
turn_df <- data.frame()
instab1_df <- data.frame()
instab2_df <- data.frame()
selinst_df <- data.frame()
for (i in seq_along(FILES)) {
  f <- FILES[i]
  lab <- LABELS[i]
  path <- file.path(results_dir, f)
  if (!file.exists(path)) {
    warning("Skipping missing file: ", path)
    next
  }
  dat <- readr::read_csv(path, show_col_types = FALSE)
  # Accept either columns named "k" and "SharpeRatio" or "k" and one method column
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
  if (!"k" %in% names(dat)) {
    warning("Skipping file missing k column: ", path)
    next
  }
  cur_df <- data.frame(k = dat$k, SR = dat[[sr_col]], label = lab)
  plot_df <- rbind(plot_df, cur_df)
  # max point for vertical line
  idx_max <- which.max(cur_df$SR)
  max_points <- rbind(max_points, data.frame(k = cur_df$k[idx_max], SR = cur_df$SR[idx_max], label = lab))

  # optional extra metrics
  if ("turnover" %in% names(dat)) {
    turn_df <- rbind(turn_df, data.frame(k = dat$k, val = dat$turnover, label = lab))
  }
  if ("weight_instability_l1" %in% names(dat)) {
    instab1_df <- rbind(instab1_df, data.frame(k = dat$k, val = dat$weight_instability_l1, label = lab))
  }
  if ("weight_instability_l2" %in% names(dat)) {
    instab2_df <- rbind(instab2_df, data.frame(k = dat$k, val = dat$weight_instability_l2, label = lab))
  }
  if ("selection_instability" %in% names(dat)) {
    selinst_df <- rbind(selinst_df, data.frame(k = dat$k, val = dat$selection_instability, label = lab))
  }
}

if (nrow(plot_df) == 0) stop("No data to plot.")

p <- ggplot(plot_df, aes(x = k, y = SR, color = label)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_vline(data = max_points, aes(xintercept = k, color = label), linetype = "dashed", alpha = 0.8) +
  labs(x = "Number of holdings k", y = "OOS Sharpe ratio", color = "In-sample window") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

fig_path <- file.path(fig_dir, FIG_NAME)
ggsave(fig_path, p, width = 8, height = 5, dpi = 150)
message("Saved plot to: ", fig_path)

# Additional plots if data available
plot_metric <- function(df, ylab, suffix) {
  if (nrow(df) == 0) return(NULL)
  p <- ggplot(df, aes(x = k, y = val, color = label)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(x = "Number of holdings k", y = ylab, color = "In-sample window") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14))
  fp <- file.path(fig_dir, sub("\\.png$", paste0("_", suffix, ".png"), FIG_NAME))
  ggsave(fp, p, width = 8, height = 5, dpi = 150)
  message("Saved plot to: ", fp)
}

plot_metric(turn_df, "Median turnover", "turnover")
plot_metric(instab1_df, "Median weight instability (ℓ1)", "weight_instability_l1")
plot_metric(instab2_df, "Median weight instability (ℓ2)", "weight_instability_l2")
plot_metric(selinst_df, "Median selection instability", "selection_instability")
