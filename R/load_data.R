#' Load Returns Data
#'
#' This function loads and assembles returns data based on the specified type.
#' If \code{type = 0}, it loads all double-sorted and industry portfolios,
#' and combines their returns into a single numeric matrix (columns represent portfolios).
#' If \code{type = 1}, it loads the CRSP returns data, and returns
#' a numeric matrix of individual asset returns.
#' All data are filtered to include only observations with dates from 05/2008 to 12/2022.
#'
#' @param type An integer indicating the type of data to load. Use \code{0} (default) for portfolio data and
#' \code{1} for individual asset data (CRSP returns).
#' @param do_checks Logical. If \code{TRUE} (default), the function will check that \code{type} is valid.
#'
#' @return A numeric matrix containing the returns data.
#'   For type = 0: Each column corresponds to a portfolio's excess returns, with the Date column removed.
#'   For type = 1: The matrix contains individual asset excess returns from CRSP, with the Date column removed.
#'
#' @details
#' When \code{type = 0}, the function expects the following objects to be present in the current environment:
#' \code{returns_mebeme25}, \code{returns_opinv25}, \code{returns_bemeinv25}, \code{returns_bemeop25},
#' \code{returns_meac25}, \code{returns_mebeta25}, \code{returns_meinv25}, \code{returns_meni25},
#' \code{returns_meprior10}, \code{returns_meprior122}, \code{returns_meprior6013},
#' \code{returns_mevar25}, and \code{returns_ind17}. For \code{type = 1}, the object
#' \code{CRSP_returns} must be available.
#'
#'
#' @examples
#' \dontrun{
#' # Load portfolio returns (double-sorted and industry portfolios)
#' portfolio_returns <- load_data(type = 0)
#'
#' # Load individual asset returns from CRSP
#' crsp_returns <- load_data(type = 1)
#' }
#'
#' @export
load_data <- function(type = 0, do_checks = TRUE) {
  # Check that type is valid
  if (do_checks) {
    if (!is.numeric(type) || !(type %in% c(0, 1))) {
      stop("Invalid 'type'. It must be 0 (for portfolios) or 1 (for individual assets).")
    }
  }

  if (type == 0) {
    # Define the names of portfolio objects to load
    portfolio_names <- c("returns_mebeme25", "returns_opinv25", "returns_bemeinv25",
                         "returns_bemeop25", "returns_meac25", "returns_mebeta25",
                         "returns_meinv25", "returns_meprior10", "returns_meprior122",
                         "returns_meprior6013", "returns_mevar25",
                         "returns_ind17")

    returns <- cbind(
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

    return(returns)

  } else if (type == 1) {

    returns <- SparsePortfolioSelection::returns_crsp[,-1, drop = FALSE]

    return(returns)
  }
}
