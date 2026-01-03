# Illustrative Sharpe ratio loss tradeoff plot (theoretical rates)
# Curves only (no dots/squares on the curves), with a smooth-looking line.

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}

# -----------------------------
# User-set parameters (edit)
# -----------------------------
N     <- 384   # number of assets
T_in  <- 500   # sample size (e.g., months)
k_max <- 20    # x-axis max

# -----------------------------
# Smooth grid for "curvier" lines
# -----------------------------
k_grid <- seq(1, k_max, length.out = 400)

efficiency_loss <- 1 / k_grid
estimation_risk <- sqrt(k_grid * log(N) / T_in)   # (k log N / T)^{1/2}

df_long <- data.frame(
  k = rep(k_grid, times = 2),
  component = factor(
    rep(c("Efficiency", "Estimation"), each = length(k_grid)),
    levels = c("Efficiency", "Estimation")
  ),
  value = c(efficiency_loss, estimation_risk)
)

# Intersection (continuous) as an "optimal" illustration:
# 1/k = sqrt(k log N / T)  =>  k^3 = T / log N
k_star <- (T_in / log(N))^(1/3)
y_star <- 1 / k_star

# -----------------------------
# Style constants (from your script)
# -----------------------------
base_size <- 14
line_w    <- 1.2

loss_cols <- c(
  "Efficiency" = "#006400",  # dark green
  "Estimation" = "#000000"   # black
)

# -----------------------------
# Plot
# -----------------------------
p <- ggplot2::ggplot(df_long, ggplot2::aes(x = k, y = value, color = component)) +
  ggplot2::geom_line(linewidth = line_w, lineend = "round", show.legend = FALSE) +
  ggplot2::scale_color_manual(values = loss_cols) +
  ggplot2::coord_cartesian(xlim = c(1, k_max), ylim = c(0, 1.1), expand = FALSE) +

  # k* marker + vertical dashed segment
  ggplot2::annotate(
    "segment",
    x = k_star, xend = k_star,
    y = 0, yend = y_star,
    linetype = "dashed",
    linewidth = 0.6,
    color = "gray50"
  ) +
  ggplot2::annotate(
    "point",
    x = k_star, y = y_star,
    size = 2.0,
    color = "black"
  ) +
  ggplot2::annotate(
    "text",
    x = k_star + 0.25,
    y = 0.5 * y_star,          # middle of the dashed line
    label = "k^\"*\"",
    parse = TRUE,
    hjust = 0, vjust = 0.5
  ) +

  # In-plot labels
  ggplot2::annotate(
    "text",
    x = 4.5, y = 0.6,
    label = "Efficiency loss",
    color = loss_cols["Efficiency"]
  ) +
  ggplot2::annotate(
    "text",
    x = 14.5, y = 0.5,
    label = "Estimation loss",
    color = loss_cols["Estimation"]
  ) +

  ggplot2::labs(
    x = "Number of holdings k",
    y = "Sharpe ratio loss"
  ) +

  # No ticks/labels (illustration) + keep left/bottom axis lines
  ggplot2::scale_x_continuous(breaks = NULL) +
  ggplot2::scale_y_continuous(breaks = NULL) +
  ggplot2::theme_minimal(base_size = base_size) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(linewidth = 0.4, color = "black")
  )

print(p)

ggplot2::ggsave(
  filename = file.path("inst", "theory_sharpe_ratio_loss_tradeoff.png"),
  plot = p,
  width = 12 / 2.54,
  height = 8 / 2.54,
  units = "in",
  dpi = 600
)


