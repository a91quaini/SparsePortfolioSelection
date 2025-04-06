#' Load Returns or Factors Data
#'
#' This function loads and assembles data based on the specified type. The input
#' \code{type} can be one of the following:
#'
#' "p" loads all double-sorted and industry portfolio objects from the SparsePortfolioSelection package.
#' The function removes the Date column from each portfolio object and aggregates their returns (columns represent portfolios)
#' into a single numeric matrix.
#'
#' "i" loads the CRSP individual asset returns data from the SparsePortfolioSelection package.
#' It removes the Date column and returns a numeric matrix of individual asset returns.
#'
#' "f" loads the Fama–French 5 factors (tradable version) data from the SparsePortfolioSelection package.
#' It removes the Date column and returns a numeric matrix of factors.
#'
#' All data are assumed to have been filtered to include only observations with dates from May 2008 to December 2022.
#'
#' @param type A character string indicating the type of data to load. Valid options are "p" for portfolios, "i" for individual assets, and "f" for factors.
#' @param do_checks Logical. If \code{TRUE} (the default), the function checks that \code{type} is one of "p", "i", or "f".
#'
#' @return A numeric matrix named \code{data} containing the data with the Date column removed. For \code{type = "p"}, each column corresponds to a portfolio's excess returns; for \code{type = "i"}, the matrix contains individual asset excess returns; for \code{type = "f"}, the matrix contains the Fama–French 5 factors (tradable).
#'
#' @examples
#' \dontrun{
#'   # Load portfolio returns
#'   data <- load_data(type = "p")
#'
#'   # Load individual asset returns (CRSP)
#'   data <- load_data(type = "i")
#'
#'   # Load factor data (Fama–French 5)
#'   data <- load_data(type = "f")
#' }
#'
#' @export
load_data <- function(type = "p", do_checks = TRUE) {
  # Check that type is valid
  if (do_checks) {
    if (!is.character(type) || !(type %in% c("p", "i", "f"))) {
      stop("Invalid 'type'. It must be 'p' (for portfolios), 'i' (for individual assets), or 'f' (for factors).")
    }
  }

  if (type == "p") {
    # Load portfolio objects from the SparsePortfolioSelection package and combine them.
    data <- cbind(
      SparsePortfolioSelection::returns_mebeme25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_meop25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_opinv25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_bemeinv25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_bemeop25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_meac25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_mebeta25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_meinv25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_meprior10[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_meprior122[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_meprior6013[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_mevar25[,-1, drop = FALSE],
      SparsePortfolioSelection::returns_ind17[,-1, drop = FALSE]
    )
    return(data)

  } else if (type == "i") {
    data <- SparsePortfolioSelection::returns_crsp[,-1, drop = FALSE]
    return(data)

  } else if (type == "f") {
    data <- SparsePortfolioSelection::factors_ff5[,-1, drop = FALSE]
    return(data)
  }
}
