# -------------------------------
# Create Data on French Factors, Portfolio Excess Returns, and CRSP Excess Returns
# Dates: >= 200805 and <= 202212
# -------------------------------

### Section A: FF5 and Momentum Factors
# define file path
path <- "data-raw/"
# Read the FF5 and Momentum CSV files (assumes first column "Date" in YYYYMM format)
factors_ff5 <- read.csv(file.path(path, "F-F_Research_Data_5_Factors_2x3.csv"), stringsAsFactors = FALSE)
factor_mom   <- read.csv(file.path(path, "F-F_Momentum_Factor.csv"), stringsAsFactors = FALSE)

# Filter rows for the desired date range
factors_ff5 <- subset(factors_ff5, as.numeric(Date) >= 200805 & as.numeric(Date) <= 202212)
factor_mom  <- subset(factor_mom,  as.numeric(Date) >= 200805 & as.numeric(Date) <= 202212)

# Convert return columns (all except Date) from percentages to decimals
factors_ff5[,-1] <- factors_ff5[,-1] / 100
factor_mom[,-1]  <- factor_mom[,-1] / 100

# Convert data.frames to matrices (so that all columns are numeric)
factors_ff5 <- as.matrix(factors_ff5)
factor_mom  <- as.matrix(factor_mom)

# Extract the risk-free factor from FF5 (assumes columns: Date, Mkt-RF, SMB, HML, RMW, CMA, RF)
rf <- factors_ff5[, c(1, 7)]

# Save the FF5 factors, Momentum factors, and RF from FF5
usethis::use_data(factors_ff5, factor_mom, rf, overwrite = TRUE)



################################################################################
### Section B: CRSP Returns Excess Returns
# Read CRSP returns CSV using data.table::fread and convert to a tibble for dplyr piping
returns_crsp <- read.csv(file.path(path, "CRSP_returns.csv"), stringsAsFactors = FALSE)

# Compute excess returns:
# Convert asset returns (in percent) to decimals and subtract the risk-free rate (row-wise)
returns_crsp[,-1] <- as.matrix(returns_crsp[,-1]) / 100 - rf[,-1]

# Convert the result to a matrix
returns_crsp <- as.matrix(returns_crsp)

# Save the CRSP returns matrix
usethis::use_data(returns_crsp, overwrite = TRUE)


## --- Section C: Additional Portfolio Files ---

# 1. returns_ind17 from "17_Industry_Portfolios.CSV"
returns_ind17 <- read.csv(file.path(path, "17_Industry_Portfolios.CSV"), stringsAsFactors = FALSE)
returns_ind17 <- returns_ind17[as.numeric(returns_ind17$Date) >= 200805 & as.numeric(returns_ind17$Date) <= 202212, ]
returns_ind17[,-1] <- (returns_ind17[,-1] / 100) - rf[,-1]
returns_ind17 <- as.matrix(returns_ind17)
usethis::use_data(returns_ind17, overwrite = TRUE)

# 2. returns_bemeinv25 from "25_Portfolios_BEME_INV_5x5.CSV"
returns_bemeinv25 <- read.csv(file.path(path, "25_Portfolios_BEME_INV_5x5.CSV"), stringsAsFactors = FALSE)
returns_bemeinv25 <- returns_bemeinv25[as.numeric(returns_bemeinv25$Date) >= 200805 & as.numeric(returns_bemeinv25$Date) <= 202212, ]
returns_bemeinv25[,-1] <- (returns_bemeinv25[,-1] / 100) - rf[,-1]
returns_bemeinv25 <- as.matrix(returns_bemeinv25)
usethis::use_data(returns_bemeinv25, overwrite = TRUE)

# 3. returns_bemeop25 from "25_Portfolios_BEME_OP_5x5.CSV"
returns_bemeop25 <- read.csv(file.path(path, "25_Portfolios_BEME_OP_5x5.CSV"), stringsAsFactors = FALSE)
returns_bemeop25 <- returns_bemeop25[as.numeric(returns_bemeop25$Date) >= 200805 & as.numeric(returns_bemeop25$Date) <= 202212, ]
returns_bemeop25[,-1] <- (returns_bemeop25[,-1] / 100) - rf[,-1]
returns_bemeop25 <- as.matrix(returns_bemeop25)
usethis::use_data(returns_bemeop25, overwrite = TRUE)

# 4. returns_meac25 from "25_Portfolios_ME_AC_5x5.csv"
returns_meac25 <- read.csv(file.path(path, "25_Portfolios_ME_AC_5x5.csv"), stringsAsFactors = FALSE)
returns_meac25 <- returns_meac25[as.numeric(returns_meac25$Date) >= 200805 & as.numeric(returns_meac25$Date) <= 202212, ]
returns_meac25[,-1] <- (returns_meac25[,-1] / 100) - rf[,-1]
returns_meac25 <- as.matrix(returns_meac25)
usethis::use_data(returns_meac25, overwrite = TRUE)

# 5. returns_mebeta25 from "25_Portfolios_ME_BETA_5x5.csv"
returns_mebeta25 <- read.csv(file.path(path, "25_Portfolios_ME_BETA_5x5.csv"), stringsAsFactors = FALSE)
returns_mebeta25 <- returns_mebeta25[as.numeric(returns_mebeta25$Date) >= 200805 & as.numeric(returns_mebeta25$Date) <= 202212, ]
returns_mebeta25[,-1] <- (returns_mebeta25[,-1] / 100) - rf[,-1]
returns_mebeta25 <- as.matrix(returns_mebeta25)
usethis::use_data(returns_mebeta25, overwrite = TRUE)

# 6. returns_meinv25 from "25_Portfolios_ME_INV_5x5.CSV"
returns_meinv25 <- read.csv(file.path(path, "25_Portfolios_ME_INV_5x5.CSV"), stringsAsFactors = FALSE)
returns_meinv25 <- returns_meinv25[as.numeric(returns_meinv25$Date) >= 200805 & as.numeric(returns_meinv25$Date) <= 202212, ]
returns_meinv25[,-1] <- (returns_meinv25[,-1] / 100) - rf[,-1]
returns_meinv25 <- as.matrix(returns_meinv25)
usethis::use_data(returns_meinv25, overwrite = TRUE)

# 7. returns_meni25 from "25_Portfolios_ME_NI_5x5.csv"
returns_meni25 <- read.csv(file.path(path, "25_Portfolios_ME_NI_5x5.csv"), stringsAsFactors = FALSE)
returns_meni25 <- returns_meni25[as.numeric(returns_meni25$Date) >= 200805 & as.numeric(returns_meni25$Date) <= 202212, ]
returns_meni25[,-1] <- (returns_meni25[,-1] / 100) - rf[,-1]
returns_meni25 <- as.matrix(returns_meni25)
usethis::use_data(returns_meni25, overwrite = TRUE)

# 8. returns_meop25 from "25_Portfolios_ME_OP_5x5.CSV"
returns_meop25 <- read.csv(file.path(path, "25_Portfolios_ME_OP_5x5.CSV"), stringsAsFactors = FALSE)
returns_meop25 <- returns_meop25[as.numeric(returns_meop25$Date) >= 200805 & as.numeric(returns_meop25$Date) <= 202212, ]
returns_meop25[,-1] <- (returns_meop25[,-1] / 100) - rf[,-1]
returns_meop25 <- as.matrix(returns_meop25)
usethis::use_data(returns_meop25, overwrite = TRUE)

# 9. returns_meprior10 from "25_Portfolios_ME_Prior_1_0.CSV"
returns_meprior10 <- read.csv(file.path(path, "25_Portfolios_ME_Prior_1_0.CSV"), stringsAsFactors = FALSE)
returns_meprior10 <- returns_meprior10[as.numeric(returns_meprior10$Date) >= 200805 & as.numeric(returns_meprior10$Date) <= 202212, ]
returns_meprior10[,-1] <- (returns_meprior10[,-1] / 100) - rf[,-1]
returns_meprior10 <- as.matrix(returns_meprior10)
usethis::use_data(returns_meprior10, overwrite = TRUE)

# 10. returns_meprior122 from "25_Portfolios_ME_Prior_12_2.CSV"
returns_meprior122 <- read.csv(file.path(path, "25_Portfolios_ME_Prior_12_2.CSV"), stringsAsFactors = FALSE)
returns_meprior122 <- returns_meprior122[as.numeric(returns_meprior122$Date) >= 200805 & as.numeric(returns_meprior122$Date) <= 202212, ]
returns_meprior122[,-1] <- (returns_meprior122[,-1] / 100) - rf[,-1]
returns_meprior122 <- as.matrix(returns_meprior122)
usethis::use_data(returns_meprior122, overwrite = TRUE)

# 11. returns_meprior6013 from "25_Portfolios_ME_Prior_60_13.CSV"
returns_meprior6013 <- read.csv(file.path(path, "25_Portfolios_ME_Prior_60_13.CSV"), stringsAsFactors = FALSE)
returns_meprior6013 <- returns_meprior6013[as.numeric(returns_meprior6013$Date) >= 200805 & as.numeric(returns_meprior6013$Date) <= 202212, ]
returns_meprior6013[,-1] <- (returns_meprior6013[,-1] / 100) - rf[,-1]
returns_meprior6013 <- as.matrix(returns_meprior6013)
usethis::use_data(returns_meprior6013, overwrite = TRUE)

# 12. returns_mevar25 from "25_Portfolios_ME_VAR_5x5.csv"
returns_mevar25 <- read.csv(file.path(path, "25_Portfolios_ME_VAR_5x5.csv"), stringsAsFactors = FALSE)
returns_mevar25 <- returns_mevar25[as.numeric(returns_mevar25$Date) >= 200805 & as.numeric(returns_mevar25$Date) <= 202212, ]
returns_mevar25[,-1] <- (returns_mevar25[,-1] / 100) - rf[,-1]
returns_mevar25 <- as.matrix(returns_mevar25)
usethis::use_data(returns_mevar25, overwrite = TRUE)

# 13. returns_opinv25 from "25_Portfolios_OP_INV_5x5.CSV"
returns_opinv25 <- read.csv(file.path(path, "25_Portfolios_OP_INV_5x5.CSV"), stringsAsFactors = FALSE)
returns_opinv25 <- returns_opinv25[as.numeric(returns_opinv25$Date) >= 200805 & as.numeric(returns_opinv25$Date) <= 202212, ]
returns_opinv25[,-1] <- (returns_opinv25[,-1] / 100) - rf[,-1]
returns_opinv25 <- as.matrix(returns_opinv25)
usethis::use_data(returns_opinv25, overwrite = TRUE)

# 2. returns_mebeme25 from "25_Portfolios_5x5.CSV"
returns_mebeme25 <- read.csv(file.path(path, "25_Portfolios_5x5.CSV"), stringsAsFactors = FALSE)
returns_mebeme25 <- returns_mebeme25[as.numeric(returns_mebeme25$Date) >= 200805 & as.numeric(returns_mebeme25$Date) <= 202212, ]
returns_mebeme25[,-1] <- (returns_mebeme25[,-1] / 100) - rf[,-1]
returns_mebeme25 <- as.matrix(returns_mebeme25)
usethis::use_data(returns_mebeme25, overwrite = TRUE)
