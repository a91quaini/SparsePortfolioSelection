# data-raw/mebeme_ind_monthly/make_mebeme_ind_monthly.R
# -----------------------------------------------------------------------------
# Clean ME×BE/ME and industry monthly portfolios (value-weighted only), replace
# sentinel missings with NA, compute excess returns in percentages, and save to
# data/mebeme_ind_monthly/.

START_DATE <- 197101L  # start at Jan 1971
STOP_DATE  <- NULL  # e.g., 202512L; set to NULL for full range

guess_dir <- function() {
  d <- tryCatch(normalizePath(dirname(sys.frame(1)$ofile)), error = function(...) NA_character_)
  if (!is.na(d) && length(list.files(d, pattern = "\\.csv$")) > 0) return(d)
  alt <- normalizePath(file.path(getwd(), "data-raw", "mebeme_ind_monthly"), mustWork = FALSE)
  if (length(list.files(alt, pattern = "\\.csv$")) > 0) return(alt)
  stop("Cannot locate data-raw/mebeme_ind_monthly directory.")
}

scriptdir <- guess_dir()
rawdir <- scriptdir
datadir <- normalizePath(file.path(scriptdir, "..", "..", "data", "mebeme_ind_monthly"), mustWork = FALSE)
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)

is_yyyymm <- function(x) {
  is.character(x) && grepl("^\\s*\\d{6}\\s*$", x)
}

read_monthly_table <- function(path) {
  df_raw <- utils::read.csv(path,
                            header = TRUE,
                            stringsAsFactors = FALSE,
                            check.names = FALSE,
                            fill = TRUE,
                            blank.lines.skip = TRUE)
  dates <- integer(nrow(df_raw))
  keep <- logical(nrow(df_raw))
  for (i in seq_len(nrow(df_raw))) {
    v <- df_raw[i, 1]
    if (is.numeric(v) && length(v) == 1L && is.finite(v)) {
      d <- as.integer(round(v))
      if (d >= 101L && d <= 999912L) {
        dates[i] <- d
        keep[i] <- TRUE
      }
    } else if (is_yyyymm(v)) {
      dates[i] <- as.integer(trimws(v))
      keep[i] <- TRUE
    }
  }
  df <- df_raw[keep, , drop = FALSE]
  if (nrow(df) == 0) stop("No YYYYMM rows found in ", path)
  df[[1]] <- dates[keep]
  if (ncol(df) >= 2) {
    for (j in 2:ncol(df)) {
      x <- df[[j]]
      if (is.character(x)) x <- gsub("%", " ", x, fixed = TRUE)
      df[[j]] <- suppressWarnings(as.numeric(x))
    }
  }
  df
}

clean_sentinels <- function(df, from_col = 2, threshold = -90) {
  num_cols <- seq(from_col, ncol(df))
  bad <- df[num_cols] < threshold | df[num_cols] <= -999
  df[num_cols][bad] <- NA_real_
  df
}

filter_window <- function(df) {
  if (!is.null(START_DATE)) df <- df[df[[1]] >= START_DATE, , drop = FALSE]
  if (!is.null(STOP_DATE)) df <- df[df[[1]] <= STOP_DATE, , drop = FALSE]
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
  out <- df_ret
  rf_vec <- rf_map[as.character(df_ret[[1]])]
  if (ncol(out) >= 2) {
    out[-1] <- sweep(df_ret[-1], 1, rf_vec, "-")
  }
  out
}

# Factors (percent)
ff5_file <- file.path(rawdir, "F-F_Research_Data_5_Factors_2x3.csv")
f_factors <- read_monthly_table(ff5_file)
f_factors <- clean_sentinels(f_factors, from_col = 2)
f_factors <- filter_window(f_factors)
rf <- rf_from_factors(f_factors)
saveRDS(rf, file.path(datadir, "rf_mebeme_ind_monthly.rds"))

portfolio_files <- c(
  "100_Portfolios_10x10.csv" = "returns_mebeme100_monthly",
  "49_Industry_Portfolios.csv" = "returns_ind49_monthly"
)

for (infile in names(portfolio_files)) {
  path <- file.path(rawdir, infile)
  if (!file.exists(path)) {
    warning("Skipping missing file: ", infile)
    next
  }
  df <- read_monthly_table(path)
  df <- clean_sentinels(df, from_col = 2)
  excess <- align_and_excess(df, rf)
  outpath <- file.path(datadir, paste0(portfolio_files[[infile]], ".rds"))
  saveRDS(excess, outpath)
  message(sprintf("✓ wrote %s (rows=%d, cols=%d)", basename(outpath), nrow(excess), ncol(excess)))
}

message("\n→ Wrote ME×BE/ME and industry monthly excess returns under data/mebeme_ind_monthly/ as .rds (percent units).")
