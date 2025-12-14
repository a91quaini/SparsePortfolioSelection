#' FF5 daily factors (with RF)
#'
#' Daily Fama–French 5 factors including RF (decimals). Not exported; used internally.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMMDD format}
#'   \item{Asset}{Factor columns}
#' }
factors_ff5_daily <- NULL

#' Daily risk-free rates (internal)
#'
#' DATE and RF (decimals) aligned to the factors. Not exported.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMMDD format}
#'   \item{Asset}{RF}
#' }
rf_daily <- NULL

#' FF5 monthly factors (with RF)
#'
#' Monthly Fama–French 5 factors including RF (decimals). Not exported; used internally.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMM format}
#'   \item{Asset}{Factor columns}
#' }
factors_ff5_monthly <- NULL

#' Monthly risk-free rates (internal)
#'
#' DATE and RF (decimals) aligned to the factors. Not exported.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMM format}
#'   \item{Asset}{RF}
#' }
rf_monthly <- NULL

#' Daily excess returns: 49 industry portfolios (US)
#'
#' DATE plus excess returns (decimals) for 49 industry portfolios.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMMDD format}
#'   \item{Asset}{Excess returns}
#' }
returns_ind49_daily <- NULL

#' Daily excess returns: ME×BE/ME 5x5 (US)
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMMDD format}
#'   \item{Asset}{Excess returns}
#' }
returns_mebeme25_daily <- NULL

#' Daily excess returns: BE/ME×INV 5x5 (US)
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMMDD format}
#'   \item{Asset}{Excess returns}
#' }
returns_bemeinv25_daily <- NULL

#' Daily excess returns: BE/ME×OP 5x5 (US)
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMMDD format}
#'   \item{Asset}{Excess returns}
#' }
returns_bemeop25_daily <- NULL

#' Daily excess returns: ME×INV 5x5 (US)
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMMDD format}
#'   \item{Asset}{Excess returns}
#' }
returns_meinv25_daily <- NULL

#' Daily excess returns: ME×OP 5x5 (US)
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMMDD format}
#'   \item{Asset}{Excess returns}
#' }
returns_meop25_daily <- NULL

#' Daily excess returns: ME×Prior 1–0 (US)
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMMDD format}
#'   \item{Asset}{Excess returns}
#' }
returns_meprior10_daily <- NULL

#' Daily excess returns: ME×Prior 12–2 (US)
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMMDD format}
#'   \item{Asset}{Excess returns}
#' }
returns_meprior122_daily <- NULL

#' Daily excess returns: ME×Prior 60–13 (US)
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMMDD format}
#'   \item{Asset}{Excess returns}
#' }
returns_meprior6013_daily <- NULL

#' Daily excess returns: OP×INV 5x5 (US)
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMMDD format}
#'   \item{Asset}{Excess returns}
#' }
returns_opinv25_daily <- NULL

# International panels (decimals, with DATE column)
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMMDD format}
#'   \item{Asset}{Returns}
#' }
returns_apxj_mebeme25_int_daily <- NULL
returns_apxj_meinv25_int_daily <- NULL
returns_apxj_meop25_int_daily <- NULL
returns_apxj_meprior25020_int_daily <- NULL
returns_eu_mebeme25_int_daily <- NULL
returns_eu_meinv25_int_daily <- NULL
returns_eu_meop25_int_daily <- NULL
returns_eu_meprior25020_int_daily <- NULL
returns_jp_mebeme25_int_daily <- NULL
returns_jp_meinv25_int_daily <- NULL
returns_jp_meop25_int_daily <- NULL
returns_jp_meprior25020_int_daily <- NULL
returns_na_mebeme25_int_daily <- NULL
returns_na_meinv25_int_daily <- NULL
returns_na_meop25_int_daily <- NULL
returns_na_meprior25020_int_daily <- NULL

# International panels (monthly, decimals, with DATE column)
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMM format}
#'   \item{Asset}{Returns}
#' }
returns_apxj_mebeme25_int_monthly <- NULL
returns_apxj_meinv25_int_monthly <- NULL
returns_apxj_meop25_int_monthly <- NULL
returns_apxj_meprior122_int_monthly <- NULL
returns_eu_mebeme25_int_monthly <- NULL
returns_eu_meinv25_int_monthly <- NULL
returns_eu_meop25_int_monthly <- NULL
returns_eu_meprior122_int_monthly <- NULL
returns_jp_mebeme25_int_monthly <- NULL
returns_jp_meinv25_int_monthly <- NULL
returns_jp_meop25_int_monthly <- NULL
returns_jp_meprior122_int_monthly <- NULL
returns_na_mebeme25_int_monthly <- NULL
returns_na_meinv25_int_monthly <- NULL
returns_na_meop25_int_monthly <- NULL
returns_na_meprior122_int_monthly <- NULL

# Monthly excess returns (US, decimals)
#' @format A data frame with columns:
#' \describe{
#'   \item{DATE}{Integer date in YYYYMM format}
#'   \item{Asset}{Excess returns}
#' }
returns_ind17_monthly <- NULL
returns_bemeinv25_monthly <- NULL
returns_bemeop25_monthly <- NULL
returns_meac25_monthly <- NULL
returns_mebeta25_monthly <- NULL
returns_meinv25_monthly <- NULL
returns_meni25_monthly <- NULL
returns_meop25_monthly <- NULL
returns_meprior10_monthly <- NULL
returns_meprior122_monthly <- NULL
returns_meprior6013_monthly <- NULL
returns_mevar25_monthly <- NULL
returns_opinv25_monthly <- NULL
returns_mebeme25_monthly <- NULL
