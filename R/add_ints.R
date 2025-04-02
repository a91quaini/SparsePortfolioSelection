#' Add two integers using C++
#'
#' @param x integer
#' @param y integer
#' @return integer sum
#' @export
add_ints <- function(x, y) {
  .Call(`_SparsePortfolioSelection_add_ints`, x, y)
}
