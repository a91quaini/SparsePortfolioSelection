library(testthat)

# --- Simulate the SparsePortfolioSelection namespace for testing ---
if (!exists("SparsePortfolioSelection")) {
  SparsePortfolioSelection <- new.env()
}

# Dummy data dimensions:
# We'll use 5 rows for all dummy matrices.
# For portfolio objects, each dummy is a 5 x 2 matrix (first column: Date, second: return).
dummy_portfolios <- list(
  returns_mebeme25    = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_meop25      = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_opinv25     = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_bemeinv25   = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_bemeop25    = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_meac25      = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_mebeta25    = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_meinv25     = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_meprior10   = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_meprior122  = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_meprior6013 = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_mevar25     = matrix(c(rep(200805, 5), runif(5)), ncol = 2),
  returns_ind17       = matrix(c(rep(200805, 5), runif(5)), ncol = 2)
)

# Assign dummy portfolio objects to the simulated namespace
for (nm in names(dummy_portfolios)) {
  assign(nm, dummy_portfolios[[nm]], envir = SparsePortfolioSelection)
}

# Dummy CRSP_returns: a 5 x 6 matrix (first column is Date, 5 return columns)
CRSP_returns <- matrix(c(rep(200805, 5), runif(5 * 5)), ncol = 6)
assign("returns_crsp", CRSP_returns, envir = SparsePortfolioSelection)

# Dummy factors_ff5: a 5 x 7 matrix (first column is Date, 6 factor columns)
factors_ff5 <- matrix(c(rep(200805, 5), runif(5 * 6)), ncol = 7)
assign("factors_ff5", factors_ff5, envir = SparsePortfolioSelection)

# --- Begin tests ---

test_that("load_data rejects invalid type", {
  expect_error(load_data(type = "x"), "Invalid 'type'")
})

test_that("load_data('p') returns a matrix with correct dimensions", {
  # For portfolios: each dummy portfolio is 5 x 2 (Date + 1 return column), so after dropping Date,
  # each contributes 1 column. There are 13 portfolio objects.
  data <- load_data(type = "p")
  expect_true(is.matrix(data))
  expect_equal(nrow(data), 176)
  expect_equal(ncol(data), 317)
})

test_that("load_data('i') returns a matrix with correct dimensions", {
  # For individual assets (CRSP): dummy CRSP_returns is 5 x 6; dropping Date leaves 5 columns.
  data <- load_data(type = "i")
  expect_true(is.matrix(data))
  expect_equal(nrow(data), 176)
  expect_equal(ncol(data), 200)
})

test_that("load_data('f') returns a matrix with correct dimensions", {
  # For factors: dummy factors_ff5 is 5 x 7; dropping Date leaves 6 columns.
  data <- load_data(type = "f")
  expect_true(is.matrix(data))
  expect_equal(nrow(data), 176)
  expect_equal(ncol(data), 6)
})
