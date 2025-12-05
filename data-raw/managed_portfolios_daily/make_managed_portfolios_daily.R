# data-raw/managed_portfolios_daily/make_managed_portfolios_daily.R
# -----------------------------------------------------------------------------
# Read daily FF5 factors + daily managed portfolios, clean sentinel missings,
# convert % → decimals, subtract RF, align on dates, and save RDS files under
# data/managed_portfolios_daily/.

START_DATE <- 19910101L  # set to NULL for full range
STOP_DATE  <- 20241231L  # set to NULL for full range

guess_dir <- function() {
  d <- tryCatch(normalizePath(dirname(sys.frame(1)$ofile)), error = function(...) NA_character_)
  if (!is.na(d) && file.exists(file.path(d, "F-F_Research_Data_5_Factors_2x3_daily.csv"))) return(d)
  alt <- normalizePath(file.path(getwd(), "data-raw", "managed_portfolios_daily"), mustWork = FALSE)
  if (file.exists(file.path(alt, "F-F_Research_Data_5_Factors_2x3_daily.csv"))) return(alt)
  stop("Cannot locate data-raw/managed_portfolios_daily directory.")
}

scriptdir <- guess_dir()
rawdir <- scriptdir
datadir <- normalizePath(file.path(scriptdir, "..", "..", "data", "managed_portfolios_daily"), mustWork = FALSE)
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)

read_raw_csv <- function(path) {
  df <- utils::read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  df[] <- lapply(df, as.numeric)
  df
}

clean_sentinels <- function(df, from_col = 2, threshold = -90) {
  num_cols <- seq(from_col, ncol(df))
  bad <- df[num_cols] < threshold
  df[num_cols][bad] <- NA_real_
  df
}

filter_window <- function(df) {
  if (!is.null(START_DATE)) df <- df[df[[1]] >= START_DATE, , drop = FALSE]
  if (!is.null(STOP_DATE)) df <- df[df[[1]] <= STOP_DATE, , drop = FALSE]
  df
}

percent_to_decimal <- function(df, from_col = 2) {
  num_cols <- seq(from_col, ncol(df))
  df[num_cols] <- df[num_cols] / 100
  df
}

rf_from_factors <- function(df) {
  rf_idx <- which(trimws(names(df)) == "RF")
  if (length(rf_idx) == 0) stop("RF column not found in factors header")
  df[, c(1, rf_idx), drop = FALSE]
}

align_and_excess <- function(df_ret, rf_tbl) {
  rf_map <- setNames(rf_tbl[[2]], rf_tbl[[1]])
  df_ret <- filter_window(df_ret)
  keep <- df_ret[[1]] %in% rf_tbl[[1]]
  df_ret <- df_ret[keep, , drop = FALSE]
  df_ret <- percent_to_decimal(df_ret)
  out <- df_ret
  rf_vec <- rf_map[as.character(df_ret[[1]])]
  out[-1] <- sweep(df_ret[-1], 1, rf_vec, "-")
  out
}

# FF5 factors (includes RF)
ff5_file <- file.path(rawdir, "F-F_Research_Data_5_Factors_2x3_daily.csv")
f_factors <- read_raw_csv(ff5_file)
f_factors <- clean_sentinels(f_factors, from_col = 2)
f_factors <- filter_window(f_factors)
f_factors <- percent_to_decimal(f_factors, from_col = 2)
rf <- rf_from_factors(f_factors)

saveRDS(f_factors, file.path(datadir, "factors_ff5_daily.rds"))
saveRDS(rf, file.path(datadir, "rf_daily.rds"))

portfolio_files <- list(
  "49_Industry_Portfolios_Daily.csv"          = "returns_ind49_daily",
  "25_Portfolios_5x5_Daily.csv"               = "returns_mebeme25_daily",
  "25_Portfolios_BEME_INV_5x5_daily.csv"      = "returns_bemeinv25_daily",
  "25_Portfolios_BEME_OP_5x5_daily.csv"       = "returns_bemeop25_daily",
  "25_Portfolios_ME_INV_5x5_daily.csv"        = "returns_meinv25_daily",
  "25_Portfolios_ME_OP_5x5_Daily.csv"         = "returns_meop25_daily",
  "25_Portfolios_ME_Prior_1_0_Daily.csv"      = "returns_meprior10_daily",
  "25_Portfolios_ME_Prior_12_2_Daily.csv"     = "returns_meprior122_daily",
  "25_Portfolios_ME_Prior_60_13_Daily.csv"    = "returns_meprior6013_daily",
  "25_Portfolios_OP_INV_5x5_daily.csv"        = "returns_opinv25_daily"
)

for (infile in names(portfolio_files)) {
  path <- file.path(rawdir, infile)
  if (!file.exists(path)) {
    warning("Skipping missing file: ", infile)
    next
  }
  df <- read_raw_csv(path)
  df <- clean_sentinels(df, from_col = 2)
  excess <- align_and_excess(df, rf)
  outpath <- file.path(datadir, paste0(portfolio_files[[infile]], ".rds"))
  saveRDS(excess, outpath)
  message(sprintf("✓ wrote %s (rows=%d, cols=%d)", basename(outpath), nrow(excess), ncol(excess)))
}

message("\n→ All daily datasets written under data/managed_portfolios_daily/ as .rds")
