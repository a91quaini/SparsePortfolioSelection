# data-raw/managed_portfolios_international_monthly/make_managed_portfolios_international_monthly.R
# -----------------------------------------------------------------------------
# Read international monthly managed-portfolio CSVs, clean sentinel missings,
# convert % → decimals, align panels to common dates, and save RDS under
# data/managed_portfolios_international_monthly/.

START_DATE <- 199007L  # set to NULL for full range
STOP_DATE  <- 202510L  # set to NULL for full range

guess_dir <- function() {
  d <- tryCatch(normalizePath(dirname(sys.frame(1)$ofile)), error = function(...) NA_character_)
  if (!is.na(d) && length(list.files(d, pattern = "\\.csv$")) > 0) return(d)
  alt <- normalizePath(file.path(getwd(), "data-raw", "managed_portfolios_international_monthly"), mustWork = FALSE)
  if (length(list.files(alt, pattern = "\\.csv$")) > 0) return(alt)
  stop("Cannot locate data-raw/managed_portfolios_international_monthly directory.")
}

scriptdir <- guess_dir()
rawdir <- scriptdir
datadir <- normalizePath(file.path(scriptdir, "..", "..", "data", "managed_portfolios_international_monthly"), mustWork = FALSE)
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)

is_yyyymm <- function(x) {
  is.character(x) && grepl("^\\s*\\d{6}\\s*$", x)
}

read_monthly_table <- function(path) {
  df_raw <- utils::read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
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

extract_dates <- function(df) as.integer(df[[1]])

align_to_dates <- function(df, dates) {
  idx <- match(dates, df[[1]])
  if (any(is.na(idx))) stop("Alignment failed: missing dates.")
  df[idx, , drop = FALSE]
}

portfolio_files <- c(
  # Asia Pacific ex Japan
  "Asia_Pacific_ex_Japan_25_Portfolios_ME_BE-ME.csv"  = "returns_apxj_mebeme25_int_monthly",
  "Asia_Pacific_ex_Japan_25_Portfolios_ME_INV.csv"    = "returns_apxj_meinv25_int_monthly",
  "Asia_Pacific_ex_Japan_25_Portfolios_ME_OP.csv"     = "returns_apxj_meop25_int_monthly",
  "Asia_Pacific_ex_Japan_25_Portfolios_ME_Prior_12_2.csv" = "returns_apxj_meprior122_int_monthly",
  "Asia_Pacific_ex_Japan_6_Portfolios_ME_BE-ME.csv"   = "returns_apxj_mebeme6_int_monthly",
  "Asia_Pacific_ex_Japan_6_Portfolios_ME_INV.csv"     = "returns_apxj_meinv6_int_monthly",
  "Asia_Pacific_ex_Japan_6_Portfolios_ME_OP.csv"      = "returns_apxj_meop6_int_monthly",
  "Asia_Pacific_ex_Japan_6_Portfolios_ME_Prior_12_2.csv" = "returns_apxj_meprior6_int_monthly",
  "Asia_Pacific_ex_Japan_32_Portfolios_ME_BE-ME_OP_2x4x4.csv" = "returns_apxj_mebemeop32_int_monthly",
  # Europe
  "Europe_25_Portfolios_ME_BE-ME.csv"                 = "returns_eu_mebeme25_int_monthly",
  "Europe_25_Portfolios_ME_INV.csv"                   = "returns_eu_meinv25_int_monthly",
  "Europe_25_Portfolios_ME_OP.csv"                    = "returns_eu_meop25_int_monthly",
  "Europe_25_Portfolios_ME_Prior_12_2.csv"            = "returns_eu_meprior122_int_monthly",
  "Europe_6_Portfolios_ME_BE-ME.csv"                  = "returns_eu_mebeme6_int_monthly",
  "Europe_6_Portfolios_ME_INV.csv"                    = "returns_eu_meinv6_int_monthly",
  "Europe_6_Portfolios_ME_OP.csv"                     = "returns_eu_meop6_int_monthly",
  "Europe_6_Portfolios_ME_Prior_12_2.csv"             = "returns_eu_meprior6_int_monthly",
  "Europe_32_Portfolios_ME_BE-ME_OP_2x4x4.csv"         = "returns_eu_mebemeop32_int_monthly",
  # Japan
  "Japan_25_Portfolios_ME_BE-ME.csv"                  = "returns_jp_mebeme25_int_monthly",
  "Japan_25_Portfolios_ME_INV.csv"                    = "returns_jp_meinv25_int_monthly",
  "Japan_25_Portfolios_ME_OP.csv"                     = "returns_jp_meop25_int_monthly",
  "Japan_25_Portfolios_ME_Prior_12_2.csv"             = "returns_jp_meprior122_int_monthly",
  "Japan_6_Portfolios_ME_BE-ME.csv"                   = "returns_jp_mebeme6_int_monthly",
  "Japan_6_Portfolios_ME_INV.csv"                     = "returns_jp_meinv6_int_monthly",
  "Japan_6_Portfolios_ME_OP.csv"                      = "returns_jp_meop6_int_monthly",
  "Japan_6_Portfolios_ME_Prior_12_2.csv"              = "returns_jp_meprior6_int_monthly",
  "Japan_32_Portfolios_ME_BE-ME_OP_2x4x4.csv"          = "returns_jp_mebemeop32_int_monthly",
  # North America
  "North_America_25_Portfolios_ME_BE-ME.csv"          = "returns_na_mebeme25_int_monthly",
  "North_America_25_Portfolios_ME_INV.csv"            = "returns_na_meinv25_int_monthly",
  "North_America_25_Portfolios_ME_OP.csv"             = "returns_na_meop25_int_monthly",
  "North_America_25_Portfolios_ME_Prior_12_2.csv"     = "returns_na_meprior122_int_monthly",
  "North_America_6_Portfolios_ME_BE-ME.csv"           = "returns_na_mebeme6_int_monthly",
  "North_America_6_Portfolios_ME_INV.csv"             = "returns_na_meinv6_int_monthly",
  "North_America_6_Portfolios_ME_OP.csv"              = "returns_na_meop6_int_monthly",
  "North_America_6_Portfolios_ME_Prior_12_2.csv"       = "returns_na_meprior6_int_monthly",
  "North_America_32_Portfolios_ME_BE-ME_OP_2x4x4.csv"   = "returns_na_mebemeop32_int_monthly"
)

factor_files <- c(
  "Asia_Pacific_ex_Japan_3_Factors.csv" = "ff3_factors_monthly_asia_pacific_ex_japan.rds",
  "Europe_3_Factors.csv"                = "ff3_factors_monthly_europe.rds",
  "Japan_3_Factors.csv"                 = "ff3_factors_monthly_japan.rds",
  "North_America_3_Factors.csv"         = "ff3_factors_monthly_north_america.rds"
)

ff_file <- file.path(rawdir, "F-F_Research_Data_5_Factors_2x3.csv")
ff_factors <- NULL
rf_tbl <- NULL

raw_tables <- list()
dates_list <- list()
factor_tables <- list()

for (infile in names(portfolio_files)) {
  path <- file.path(rawdir, infile)
  if (!file.exists(path)) {
    warning("Skipping missing file: ", infile)
    next
  }
  df <- read_monthly_table(path)
  df <- filter_window(df)
  df <- clean_sentinels(df, from_col = 2)
  raw_tables[[infile]] <- df
  dates_list[[infile]] <- extract_dates(df)
}

if (length(raw_tables) == 0) stop("No input CSVs found.")

common_dates <- Reduce(intersect, dates_list)
common_dates <- sort(common_dates)

# Save region-specific FF3 factors and build RF from the first available set
rf_built <- FALSE
for (f_in in names(factor_files)) {
  f_path <- file.path(rawdir, f_in)
  if (!file.exists(f_path)) {
    warning("Skipping missing factor file: ", f_in)
    next
  }
  f_df <- read_monthly_table(f_path)
  f_df <- clean_sentinels(f_df, from_col = 2)
  f_df <- filter_window(f_df)
  f_df <- percent_to_decimal(f_df, from_col = 2)
  factor_tables[[f_in]] <- f_df
  saveRDS(f_df, file.path(datadir, factor_files[[f_in]]))
  message(sprintf("✓ wrote %s (rows=%d, cols=%d)", factor_files[[f_in]], nrow(f_df), ncol(f_df)))
  if (!rf_built && "RF" %in% names(f_df)) {
    rf_tbl <- rf_from_factors(f_df)
    rf_tbl <- rf_tbl[rf_tbl[[1]] %in% common_dates, , drop = FALSE]
    saveRDS(rf_tbl, file.path(datadir, "rf_monthly.rds"))
    rf_built <- TRUE
  }
}
if (!rf_built) stop("No factor file provided RF; cannot build excess returns.")

for (infile in names(portfolio_files)) {
  if (!infile %in% names(raw_tables)) next
  df <- raw_tables[[infile]]
  df <- align_to_dates(df, common_dates)
  df <- percent_to_decimal(df, from_col = 2)
  # compute excess returns using common RF
  rf_map <- setNames(rf_tbl[[2]], rf_tbl[[1]])
  rf_vec <- rf_map[as.character(df[[1]])]
  if (ncol(df) >= 2) {
    df[-1] <- sweep(df[-1], 1, rf_vec, "-")
  }
  outpath <- file.path(datadir, paste0(portfolio_files[[infile]], ".rds"))
  saveRDS(df, outpath)
  message(sprintf("✓ wrote %s (rows=%d, cols=%d)", basename(outpath), nrow(df), ncol(df)))
}

message("\n→ Wrote international monthly datasets under data/managed_portfolios_international_monthly/ as .rds")
